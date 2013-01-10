"""Next-gen alignments with TopHat a spliced read mapper for RNA-seq experiments.

http://tophat.cbcb.umd.edu
"""
import sh
import os
import shutil
import subprocess
from contextlib import closing
import glob

import pysam
import numpy

from bcbio.pipeline import config_utils
from bcbio.ngsalign import bowtie, bowtie2
from bcbio.utils import safe_makedir, file_exists
from bcbio.distributed.transaction import file_transaction


galaxy_location_file = "bowtie_indices.loc"

_out_fnames = ["accepted_hits.sam", "junctions.bed",
               "insertions.bed", "deletions.bed"]


def tophat_align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
                 rg_name=None):
    """
    run alignment using Tophat v1
    """
    if _ref_version(ref_file) != 1:
        raise ValueError("Tophat v1 needs Bowtie version 1 libraries "
                         "(the extension should be .ebwt). Download "
                         "v1 libraries or build new ones with bowtie-build.")

    qual_format = config["algorithm"].get("quality_format", None)

    if qual_format is None or qual_format.lower() == "illumina":
        qual_flags = ["--solexa1.3-quals"]
    elif qual_format == "solexa":
        qual_flags = ["--solexa-quals"]
    else:
        qual_flags = []

    gtf_file = config.get("gtf", None)
    if gtf_file is not None:
        gtf_flags = ["--GTF", gtf_file]
    else:
        gtf_flags = []

    cores = config.get("resources", {}).get("tophat", {}).get("cores", None)
    core_flags = ["-p", str(cores)] if cores else []
    out_dir = os.path.join(align_dir, "%s_tophat" % out_base)
    out_file = os.path.join(out_dir, _out_fnames[0])
    files = [ref_file, fastq_file]
    if not file_exists(out_file):
        with file_transaction(out_dir) as tx_out_dir:
            safe_makedir(tx_out_dir)
            cl = [config_utils.get_program("tophat", config)]
            cl += core_flags
            cl += qual_flags
            cl += gtf_flags
            cl += ["-m", str(config["algorithm"].get("max_errors", 0)),
                   "--output-dir", tx_out_dir,
                   "--no-convert-bam"]
            if pair_file:
                d, d_stdev = _estimate_paired_innerdist(fastq_file, pair_file, ref_file,
                                                        out_base, tx_out_dir, config)
                cl += ["--mate-inner-dist", str(d),
                       "--mate-std-dev", str(d_stdev)]
                files.append(pair_file)
            cl += files
            child = subprocess.check_call(cl)
    out_file_final = os.path.join(out_dir, "%s.sam" % out_base)
    if not os.path.exists(out_file_final):
        os.symlink(out_file, out_file_final)
    return out_file_final


def tophat2_align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
                 rg_name=None):
    """
    run alignment using Tophat v2
    """

    """
    select the correct bowtie to run based on the reference version
    """

    qual_format = config["algorithm"].get("quality_format", None)

    if qual_format is None or qual_format.lower() == "illumina":
        qual_flags = ["--solexa1.3-quals"]
    elif qual_format == "solexa":
        qual_flags = ["--solexa-quals"]
    else:
        qual_flags = []

    gtf_file = config.get("gtf", None)
    if gtf_file is not None:
        gtf_flags = ["--GTF", gtf_file]
    else:
        gtf_flags = []

    if _ref_version(ref_file) == 1:
        ref_flags = ["--bowtie1"]
    else:
        ref_flags = []

    cores = config.get("resources", {}).get("tophat", {}).get("cores", None)
    core_flags = ["-p", str(cores)] if cores else []
    out_dir = os.path.join(align_dir, "%s_tophat" % out_base)
    out_file = os.path.join(out_dir, _out_fnames[0])
    files = [ref_file, fastq_file]
    if not file_exists(out_file):
        with file_transaction(out_dir) as tx_out_dir:
            safe_makedir(tx_out_dir)
            cl = [config_utils.get_program("tophat", config)]
            cl += core_flags
            cl += qual_flags
            cl += gtf_flags
            cl += ref_flags
            cl += ["-m", str(config["algorithm"].get("max_errors", 0)),
                   "--output-dir", tx_out_dir,
                   "--no-convert-bam"]
            if pair_file:
                d, d_stdev = _estimate_paired_innerdist(fastq_file, pair_file, ref_file,
                                                        out_base, tx_out_dir, config)
                cl += ["--mate-inner-dist", str(d),
                       "--mate-std-dev", str(d_stdev)]
                files.append(pair_file)
            cl += files
            child = subprocess.check_call(cl)
    out_file_final = os.path.join(out_dir, "%s.sam" % out_base)
    if not os.path.exists(out_file_final):
        os.symlink(out_file, out_file_final)
    return out_file_final


def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          rg_name=None):

    out_dir = os.path.join(align_dir, "%s_tophat" % out_base)
    out_file = os.path.join(out_dir, _out_fnames[0])

    if file_exists(out_file):
        return out_file

    if not _bowtie_ref_match(ref_file, config):
        raise ValueError("Bowtie version %d was detected but the reference "
                         "file %s is built for version %d. Download version "
                         "%d or build it with bowtie-build."
                         % (_bowtie_major_version(config), ref_file,
                            _ref_version(ref_file),
                            _bowtie_major_version(config)))

    tophat_version = _tophat_major_version(config)
    if tophat_version == 1:
        out_files = tophat_align(fastq_file, pair_file, ref_file, out_base,
                                 align_dir, config, rg_name=None)
    elif tophat_version == 2:
        out_files = tophat2_align(fastq_file, pair_file, ref_file, out_base,
                                  align_dir, config, rg_name=None)
    else:
        raise ValueError("Tophat major version %d is not supported yet."
                         % (tophat_version))

    return out_files


def _estimate_paired_innerdist(fastq_file, pair_file, ref_file, out_base,
                               out_dir, config):
    """Use Bowtie to estimate the inner distance of paired reads.
    """
    # skip initial reads for large file, but not for smaller
    dists = _bowtie_for_innerdist("1000000", fastq_file, pair_file, ref_file,
                                  out_base, out_dir, config)
    if len(dists) == 0:
        dists = _bowtie_for_innerdist("1", fastq_file, pair_file, ref_file,
                                      out_base, out_dir, config, True)
    return int(round(numpy.mean(dists))), int(round(numpy.std(dists)))


def _bowtie_for_innerdist(start, fastq_file, pair_file, ref_file, out_base,
                          out_dir, config, remove_workdir=False):
    work_dir = os.path.join(out_dir, "innerdist_estimate")
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    safe_makedir(work_dir)
    extra_args = ["-s", str(start), "-u", "250000"]
    bowtie_runner = _select_bowtie_version(config)
    out_sam = bowtie_runner.align(fastq_file, pair_file, ref_file, out_base,
                                  work_dir, config, extra_args)
    dists = []
    with closing(pysam.Samfile(out_sam)) as work_sam:
        for read in work_sam:
            if not read.is_unmapped and read.is_read1:
                dists.append(abs(read.isize) - 2 * read.rlen)
    return dists


def _bowtie_major_version(config):
    bowtie_runner = sh.Command(config_utils.get_program("bowtie", config,
                                                        default="tophat"))
    """
    bowtie --version returns strings like this:
    bowtie version 0.12.7
    32-bit
    Built on Franklin.local
    Tue Sep  7 14:25:02 PDT 2010
    """
    version_line = str(bowtie_runner(version=True)).split("\n")[0]
    version_string = version_line.strip().split()[2]
    major_version = int(version_string.split(".")[0])
    # bowtie version 1 has a leading character of 0
    if major_version == 0:
        major_version += 1
    return major_version


def _tophat_major_version(config):
    tophat_runner = sh.Command(config_utils.get_program("tophat", config,
                                                        default="tophat"))

    # tophat --version returns strings like this: Tophat v2.0.4
    version_string = str(tophat_runner(version=True)).strip().split()[1]
    major_version = int(version_string.split(".")[0][1:])
    return major_version


def _bowtie_ref_match(ref_file, config):
    return _ref_version(ref_file) == _bowtie_major_version(config)


def _select_bowtie_version(config):
    if _bowtie_major_version(config) == 1:
        return bowtie
    else:
        return bowtie2


def _ref_version(ref_file):
    _, ext = os.path.splitext(glob.glob(ref_file + "*")[0])
    if ext == ".ebwt":
        return 1
    elif ext == ".bt2":
        return 2
    else:
        raise ValueError("Cannot detect which reference version %s is. "
                         "Should end in either .ebwt (bowtie) or .bt2 "
                         "(bowtie2)." % (ref_file))
