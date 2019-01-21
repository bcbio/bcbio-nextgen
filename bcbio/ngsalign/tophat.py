"""Next-gen alignments with TopHat a spliced read mapper for RNA-seq experiments.

http://tophat.cbcb.umd.edu
"""
import os
import shutil
import glob
import subprocess

import numpy
import pysam

from bcbio.pipeline import config_utils
from bcbio.ngsalign import bowtie, bowtie2
from bcbio.utils import safe_makedir, file_exists, get_in, symlink_plus
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import bam, broad, utils
import bcbio.pipeline.datadict as dd


def _set_quality_flag(options, data):
    qual_format = dd.get_quality_format(data)
    if qual_format.lower() == "illumina":
        options["solexa1.3-quals"] = True
    elif qual_format.lower() == "solexa":
        options["solexa-quals"] = True
    return options

def _set_transcriptome_option(options, data, ref_file):
    # prefer transcriptome-index vs a GTF file if available
    transcriptome_index = get_in(data, ("genome_resources", "rnaseq",
                                        "transcriptome_index", "tophat"))
    fusion_mode = _should_run_fusion(data)
    if transcriptome_index and file_exists(transcriptome_index) and not fusion_mode:
        options["transcriptome-index"] = os.path.splitext(transcriptome_index)[0]
        return options

    gtf_file = dd.get_gtf_file(data)
    if gtf_file:
        options["GTF"] = gtf_file
        return options

    return options

def _set_cores(options, config):
    num_cores = config["algorithm"].get("num_cores", 0)
    if num_cores > 1 and "num-threads" not in options:
        options["num-threads"] = num_cores
    return options

def _set_rg_options(options, names):
    if not names:
        return options
    options["rg-id"] = names["rg"]
    options["rg-sample"] = names["sample"]
    options["rg-library"] = names["lb"] or names["pl"]
    options["rg-platform-unit"] = names["pu"]
    return options

def _set_stranded_flag(options, config):
    strand_flag = {"unstranded": "fr-unstranded",
                   "firststrand": "fr-firststrand",
                   "secondstrand": "fr-secondstrand"}
    stranded = get_in(config, ("algorithm", "strandedness"), "unstranded").lower()
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
                                     "Valid values are 'firststrand', "
                                     "'secondstrand' and 'unstranded" % (stranded))
    flag = strand_flag[stranded]
    options["library-type"] = flag
    return options

def _set_fusion_mode(options, data):
    if _should_run_fusion(data):
        options["fusion-search"] = True
    return options

def tophat_align(fastq_file, pair_file, ref_file, out_base, align_dir, data,
                 names=None):
    """
    run alignment using Tophat v2
    """
    config = data["config"]
    options = get_in(config, ("resources", "tophat", "options"), {})
    options = _set_fusion_mode(options, data)
    options = _set_quality_flag(options, data)
    options = _set_transcriptome_option(options, data, ref_file)
    options = _set_cores(options, config)
    options = _set_rg_options(options, names)
    options = _set_stranded_flag(options, config)

    ref_file, runner = _determine_aligner_and_reference(ref_file, data)

    # fusion search does not work properly with Bowtie2
    if options.get("fusion-search", False):
        ref_file = ref_file.replace("/bowtie2", "/bowtie")

    if _tophat_major_version(config) == 1:
        raise NotImplementedError("Tophat versions < 2.0 are not supported, please "
                                  "download the newest version of Tophat here: "
                                  "http://tophat.cbcb.umd.edu")

    if _ref_version(ref_file) == 1 or options.get("fusion-search", False):
        options["bowtie1"] = True

    out_dir = os.path.join(align_dir, "%s_tophat" % out_base)
    final_out = os.path.join(out_dir, "{0}.bam".format(names["sample"]))
    if file_exists(final_out):
        return final_out

    out_file = os.path.join(out_dir, "accepted_hits.bam")
    unmapped = os.path.join(out_dir, "unmapped.bam")
    files = [ref_file, fastq_file]
    if not file_exists(out_file):
        with file_transaction(config, out_dir) as tx_out_dir:
            safe_makedir(tx_out_dir)
            if pair_file and not options.get("mate-inner-dist", None):
                d, d_stdev = _estimate_paired_innerdist(fastq_file, pair_file,
                                                        ref_file, out_base,
                                                        tx_out_dir, data)
                options["mate-inner-dist"] = d
                options["mate-std-dev"] = d_stdev
                files.append(pair_file)
            options["output-dir"] = tx_out_dir
            options["no-coverage-search"] = True
            options["no-mixed"] = True
            cmd = [utils.get_program_python("tophat"), config_utils.get_program("tophat", config)]
            for k, v in options.items():
                if v is True:
                    cmd.append("--%s" % k)
                else:
                    assert not isinstance(v, bool)
                    cmd.append("--%s=%s" % (k, v))
            # tophat requires options before arguments, otherwise it silently ignores them
            cmd += files
            do.run(cmd, "Running Tophat on %s and %s." % (fastq_file, pair_file))
    if pair_file and _has_alignments(out_file):
        fixed = _fix_mates(out_file, os.path.join(out_dir, "%s-align.bam" % out_base),
                           ref_file, config)
    else:
        fixed = out_file
    fixed_unmapped = _fix_unmapped(fixed, unmapped, data)
    fixed = merge_unmapped(fixed, fixed_unmapped, config)
    fixed = _add_rg(fixed, config, names)
    fixed = bam.sort(fixed, config)
    picard = broad.runner_from_path("picard", config)
    # set the contig order to match the reference file so GATK works
    fixed = picard.run_fn("picard_reorder", fixed, data["sam_ref"],
                          os.path.splitext(fixed)[0] + ".picard.bam")
    fixed = fix_insert_size(fixed, config)
    if not file_exists(final_out):
        symlink_plus(fixed, final_out)
    return final_out

def merge_unmapped(bam_file, unmapped_bam, config):
    merged_bam = os.path.join(os.path.dirname(bam_file), "merged.bam")
    if not file_exists(merged_bam):
        merged_bam = bam.merge([bam_file, unmapped_bam], merged_bam, config)
    return merged_bam

def _has_alignments(sam_file):
    with open(sam_file) as in_handle:
        try:
            for line in in_handle:
                if line.startswith("File removed to save disk space"):
                    return False
                elif not line.startswith("@"):
                    return True
        except UnicodeDecodeError:
            return not bam.is_empty(sam_file)
    return False

def _fix_mates(orig_file, out_file, ref_file, config):
    """Fix problematic unmapped mate pairs in TopHat output.

    TopHat 2.0.9 appears to have issues with secondary reads:
    https://groups.google.com/forum/#!topic/tuxedo-tools-users/puLfDNbN9bo
    This cleans the input file to only keep properly mapped pairs,
    providing a general fix that will handle correctly mapped secondary
    reads as well.
    """
    if not file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            samtools = config_utils.get_program("samtools", config)
            cmd = "{samtools} view -bS -h -t {ref_file}.fai -F 8 {orig_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "Fix mate pairs in TopHat output", {})
    return out_file

def _add_rg(unmapped_file, config, names):
    """Add the missing RG header."""
    picard = broad.runner_from_path("picard", config)
    rg_fixed = picard.run_fn("picard_fix_rgs", unmapped_file, names)
    return rg_fixed


def _fix_unmapped(mapped_file, unmapped_file, data):
    """
    The unmapped.bam file up until at least Tophat 2.1.1 is broken in various
    ways, see https://github.com/cbrueffer/tophat-recondition for details.
    Run TopHat-Recondition to fix these issues.
    """
    out_file = os.path.splitext(unmapped_file)[0] + "_fixup.bam"
    if file_exists(out_file):
        return out_file

    assert os.path.dirname(mapped_file) == os.path.dirname(unmapped_file)

    cmd = config_utils.get_program("tophat-recondition", data)
    cmd += " -q"
    tophat_out_dir = os.path.dirname(mapped_file)
    tophat_logfile = os.path.join(tophat_out_dir, 'tophat-recondition.log')

    with file_transaction(data, tophat_logfile) as tx_logfile:
        cmd += ' --logfile %s' % tx_logfile
        cmd += " -m %s" % mapped_file
        cmd += " -u %s" % unmapped_file
        cmd += " %s" % tophat_out_dir
        do.run(cmd, "Fixing unmapped reads with Tophat-Recondition.", None)

    return out_file


def align(fastq_file, pair_file, ref_file, names, align_dir, data,):
    out_files = tophat_align(fastq_file, pair_file, ref_file, names["lane"],
                             align_dir, data, names)

    return out_files


def _estimate_paired_innerdist(fastq_file, pair_file, ref_file, out_base,
                               out_dir, data):
    """Use Bowtie to estimate the inner distance of paired reads.
    """
    mean, stdev = _bowtie_for_innerdist("100000", fastq_file, pair_file, ref_file,
                                        out_base, out_dir, data, True)
    if not mean or not stdev:
        mean, stdev = _bowtie_for_innerdist("1", fastq_file, pair_file, ref_file,
                                            out_base, out_dir, data, True)
    # No reads aligning so no data to process, set some default values
    if not mean or not stdev:
        mean, stdev = 200, 50

    return mean, stdev


def _bowtie_for_innerdist(start, fastq_file, pair_file, ref_file, out_base,
                          out_dir, data, remove_workdir=False):
    work_dir = os.path.join(out_dir, "innerdist_estimate")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    safe_makedir(work_dir)
    extra_args = ["-s", str(start), "-u", "250000"]
    ref_file, bowtie_runner = _determine_aligner_and_reference(ref_file, data)
    out_sam = bowtie_runner.align(fastq_file, pair_file, ref_file, {"lane": out_base},
                                  work_dir, data, extra_args)
    dists = []
    with pysam.Samfile(out_sam) as work_sam:
        for read in work_sam:
            if read.is_proper_pair and read.is_read1:
                dists.append(abs(read.isize) - 2 * read.rlen)
    if dists:
        median = float(numpy.median(dists))
        deviations = []
        for d in dists:
            deviations.append(abs(d - median))
        # this is the median absolute deviation estimator of the
        # standard deviation
        mad = 1.4826 * float(numpy.median(deviations))
        return int(median), int(mad)
    else:
        return None, None

def _calculate_average_read_length(sam_file):
    with pysam.Samfile(sam_file) as work_sam:
        count = 0
        read_lengths = []
        for read in work_sam:
            count = count + 1
            read_lengths.append(read.rlen)
    avg_read_length = int(float(sum(read_lengths)) / float(count))
    return avg_read_length


def _bowtie_major_version(stdout):
    """
    bowtie --version returns strings like this:
    bowtie version 0.12.7
    32-bit
    Built on Franklin.local
    Tue Sep  7 14:25:02 PDT 2010
    """
    version_line = stdout.split("\n")[0]
    version_string = version_line.strip().split()[2]
    major_version = int(version_string.split(".")[0])
    # bowtie version 1 has a leading character of 0 or 1
    if major_version == 0 or major_version == 1:
        major_version = 1
    return major_version


def _should_run_fusion(data):
    return dd.get_fusion_caller(data)

def _determine_aligner_and_reference(ref_file, data):
    fusion_mode = _should_run_fusion(data)
    # fusion_mode only works with bowtie1
    if fusion_mode:
        return _get_bowtie_with_reference(ref_file, 1)
    else:
        return _get_bowtie_with_reference(ref_file, 2)

def _get_bowtie_with_reference(ref_file, version):
    if version == 1:
        ref_file = ref_file.replace("/bowtie2/", "/bowtie/")
        return ref_file, bowtie
    else:
        ref_file = ref_file.replace("/bowtie/", "/bowtie2/")
        return ref_file, bowtie2


def _tophat_major_version(config):
    cmd =  [
        utils.get_program_python("tophat"),
        config_utils.get_program("tophat", config, default="tophat"),
        "--version"
    ]

    # tophat --version returns strings like this: Tophat v2.0.4
    version_string = str(subprocess.check_output(cmd)).strip().split()[1]
    major_version = int(version_string.split(".")[0][1:])
    return major_version


def _ref_version(ref_file):
    for ext in [os.path.splitext(x)[1] for x in glob.glob(ref_file + "*")]:
        if ext == ".ebwt":
            return 1
        elif ext == ".bt2":
            return 2
    raise ValueError("Cannot detect which reference version %s is. "
                     "Should end in either .ebwt (bowtie) or .bt2 "
                     "(bowtie2)." % (ref_file))

def fix_insert_size(in_bam, config):
    """
    Tophat sets PI in the RG to be the inner distance size, but the SAM spec
    states should be the insert size. This fixes the RG in the alignment
    file generated by Tophat header to match the spec
    """
    fixed_file = os.path.splitext(in_bam)[0] + ".pi_fixed.bam"
    if file_exists(fixed_file):
        return fixed_file
    header_file = os.path.splitext(in_bam)[0] + ".header.sam"
    read_length = bam.estimate_read_length(in_bam)
    bam_handle= bam.open_samfile(in_bam)
    header = bam_handle.header.copy()
    rg_dict = header['RG'][0]
    if 'PI' not in rg_dict:
        return in_bam
    PI = int(rg_dict.get('PI'))
    PI = PI + 2*read_length
    rg_dict['PI'] = PI
    header['RG'][0] = rg_dict
    with pysam.Samfile(header_file, "wb", header=header) as out_handle:
        with bam.open_samfile(in_bam) as in_handle:
            for record in in_handle:
                out_handle.write(record)
    shutil.move(header_file, fixed_file)
    return fixed_file
