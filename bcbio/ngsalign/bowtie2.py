"""Next gen sequence alignments with Bowtie2.

http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
"""
import os

from bcbio.pipeline import config_utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio import bam, utils
from bcbio.pipeline import datadict as dd
from bcbio.rnaseq import gtf
from bcbio.ngsalign import alignprep, postalign

def _bowtie2_args_from_config(config, curcl):
    """Configurable high level options for bowtie2.
    """
    qual_format = config["algorithm"].get("quality_format", "")
    if qual_format.lower() == "illumina":
        qual_flags = ["--phred64-quals"]
    else:
        qual_flags = []
    num_cores = config["algorithm"].get("num_cores", 1)
    core_flags = ["-p", str(num_cores)] if num_cores > 1 else []
    user_opts = config_utils.get_resources("bowtie2", config).get("options", [])
    for flag_opt in (o for o in user_opts if str(o).startswith("-")):
        if flag_opt in curcl:
            raise ValueError("Duplicate option %s in resources and bcbio commandline: %s %s" %
                             flag_opt, user_opts, curcl)
    return core_flags + qual_flags + user_opts

def align(fastq_file, pair_file, ref_file, names, align_dir, data,
          extra_args=None):
    """Alignment with bowtie2.
    """
    config = data["config"]
    analysis_config = ANALYSIS.get(data["analysis"].lower())
    assert analysis_config, "Analysis %s is not supported by bowtie2" % (data["analysis"])
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(dd.get_sample_name(data)))
    if data.get("align_split"):
        final_file = out_file
        out_file, data = alignprep.setup_combine(final_file, data)
        fastq_file, pair_file = alignprep.split_namedpipe_cls(fastq_file, pair_file, data)
    else:
        final_file = None
    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        with postalign.tobam_cl(data, out_file, pair_file is not None) as (tobam_cl, tx_out_file):
            cl = [config_utils.get_program("bowtie2", config)]
            cl += extra_args if extra_args is not None else []
            cl += ["-q",
                   "-x", ref_file]
            cl += analysis_config.get("params", [])
            if pair_file:
                cl += ["-1", fastq_file, "-2", pair_file]
            else:
                cl += ["-U", fastq_file]
            if names and "rg" in names:
                cl += ["--rg-id", names["rg"]]
                for key, tag in [("sample", "SM"), ("pl", "PL"), ("pu", "PU"), ("lb", "LB")]:
                    if names.get(key):
                        cl += ["--rg", "%s:%s" % (tag, names[key])]
            cl += _bowtie2_args_from_config(config, cl)
            cl = [str(i) for i in cl]
            cmd = "unset JAVA_HOME && " + " ".join(cl) + " | " + tobam_cl
            do.run(cmd, "Aligning %s and %s with Bowtie2." % (fastq_file, pair_file))
    return out_file

# Optional galaxy location file. Falls back on remap_index_fn if not found
galaxy_location_file = "bowtie2_indices.loc"

def remap_index_fn(ref_file):
    """Map sequence references to equivalent bowtie2 indexes.
    """
    return os.path.splitext(ref_file)[0].replace("/seq/", "/bowtie2/")

def filter_multimappers(align_file, data):
    """
    It does not seem like bowtie2 has a corollary to the -m 1 flag in bowtie,
    there are some options that are close but don't do the same thing. Bowtie2
    sets the XS flag for reads mapping in more than one place, so we can just
    filter on that. This will not work for other aligners.
    """
    config = dd.get_config(data)
    type_flag = "" if bam.is_bam(align_file) else "S"
    base, ext = os.path.splitext(align_file)
    out_file = base + ".unique" + ext
    bed_file = dd.get_variant_regions(data) or dd.get_sample_callable(data)
    bed_cmd = '-L {0}'.format(bed_file) if bed_file else " "
    if utils.file_exists(out_file):
        return out_file
    base_filter = '-F "[XS] == null and not unmapped {paired_filter} and not duplicate" '
    if bam.is_paired(align_file):
        paired_filter = "and paired and proper_pair"
    else:
        paired_filter = ""
    filter_string = base_filter.format(paired_filter=paired_filter)
    sambamba = config_utils.get_program("sambamba", config)
    num_cores = dd.get_num_cores(data)
    with file_transaction(out_file) as tx_out_file:
        cmd = ('{sambamba} view -h{type_flag} '
               '--nthreads {num_cores} '
               '-f bam {bed_cmd} '
               '{filter_string} '
               '{align_file} '
               '> {tx_out_file}')
        message = "Removing multimapped reads from %s." % align_file
        do.run(cmd.format(**locals()), message)
    bam.index(out_file, config)
    return out_file

ANALYSIS = {"chip-seq": {"params": ["-X", 2000]},
            "variant2": {"params": ["-X", 2000]},
            "standard": {"params": ["-X", 2000]},
            "rna-seq": {"params": ["--sensitive", "-X", 2000]},
            "smallrna-seq": {"params": ["-N", 1, "-k", 1000, "--sensitive", "-X", 200]}}

def index_transcriptome(gtf_file, ref_file, data):
    """
    use a GTF file and a reference FASTA file to index the transcriptome
    """
    gtf_fasta = gtf.gtf_to_fasta(gtf_file, ref_file)
    bowtie2_index = os.path.splitext(gtf_fasta)[0]
    bowtie2_build = config_utils.get_program("bowtie2", data["config"]) + "-build"
    cmd = "{bowtie2_build} --offrate 1 {gtf_fasta} {bowtie2_index}".format(**locals())
    message = "Creating transcriptome index of %s with bowtie2." % (gtf_fasta)
    do.run(cmd, message)
    return bowtie2_index


def align_transcriptome(fastq_file, pair_file, ref_file, data):
    """
    bowtie2 with settings for aligning to the transcriptome for eXpress/RSEM/etc
    """
    work_bam = dd.get_work_bam(data)
    base, ext = os.path.splitext(work_bam)
    out_file = base + ".transcriptome" + ext
    if utils.file_exists(out_file):
        data = dd.set_transcriptome_bam(data, out_file)
        return data
    bowtie2 = config_utils.get_program("bowtie2", data["config"])
    gtf_file = dd.get_gtf_file(data)
    gtf_index = index_transcriptome(gtf_file, ref_file, data)
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    fastq_cmd = "-1 %s" % fastq_file if pair_file else "-U %s" % fastq_file
    pair_cmd = "-2 %s " % pair_file if pair_file else ""
    cmd = ("{bowtie2} -p {num_cores} -a -X 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed -x {gtf_index} {fastq_cmd} {pair_cmd} ")
    with file_transaction(data, out_file) as tx_out_file:
        message = "Aligning %s and %s to the transcriptome." % (fastq_file, pair_file)
        cmd += "| " + postalign.sam_to_sortbam_cl(data, tx_out_file, name_sort=True)
        do.run(cmd.format(**locals()), message)
    data = dd.set_transcriptome_bam(data, out_file)
    return data
