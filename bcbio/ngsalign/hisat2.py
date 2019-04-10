import os
from bcbio.utils import file_exists, safe_makedir
from bcbio import bam
from bcbio import bed
import bcbio.pipeline.datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.ngsalign import alignprep, postalign
from bcbio.provenance import do

def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    paired = True if pair_file else False
    hisat2 = config_utils.get_program("hisat2", data)
    num_cores = dd.get_num_cores(data)
    quality_flag = _get_quality_flag(data)
    stranded_flag = _get_stranded_flag(data, paired)
    rg_flags = _get_rg_flags(names)
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(dd.get_sample_name(data)))
    if data.get("align_split"):
        final_file = out_file
        out_file, data = alignprep.setup_combine(final_file, data)
        fastq_file, pair_file = alignprep.split_namedpipe_cls(fastq_file, pair_file, data)
    else:
        final_file = None
    if not file_exists(out_file) and (final_file is None or not file_exists(final_file)):
        cmd = ("{hisat2} --new-summary -x {ref_file} -p {num_cores} {quality_flag} {stranded_flag} "
               "{rg_flags} ")
        if paired:
            cmd += "-1 {fastq_file} -2 {pair_file} "
        else:
            cmd += "-U {fastq_file} "
        if dd.get_analysis(data).lower() == "smallrna-seq":
            cmd += "-k 1000 "
        # if assembling transcripts, set flags that cufflinks/stringtie can use
        if dd.get_transcript_assembler(data):
            cmd += "--dta-cufflinks "
        if dd.get_analysis(data).lower() == "rna-seq":
            splicesites = get_known_splicesites_file(align_dir, data)
            if file_exists(splicesites):
                cmd += "--known-splicesite-infile {splicesites} "
        novel_splicesite_file = os.path.join(align_dir, "{0}-novelsplicesites.bed".format(dd.get_sample_name(data)))
        cmd += "--novel-splicesite-outfile {novel_splicesite_file} "
        # apply additional hisat2 options
        cmd += " ".join(_get_options_from_config(data))

        message = "Aligning %s and %s with hisat2." % (fastq_file, pair_file)
        with postalign.tobam_cl(data, out_file, pair_file is not None) as (tobam_cl, tx_out_file):
            cmd += " | " + tobam_cl
            do.run(cmd.format(**locals()), message)
    data = dd.set_work_bam(data, out_file)
    junctionbed = get_splicejunction_file(align_dir, data)
    data = dd.set_junction_bed(data, junctionbed)
    return data

def get_known_splicesites_file(align_dir, data):
    gtf_file = dd.get_gtf_file(data)
    splicesites = os.path.join(os.path.dirname(gtf_file),
                               "ref-transcripts-splicesites.txt")
    if not file_exists(splicesites):
        splicesites = create_splicesites_file(gtf_file, align_dir, data)
    return splicesites

def create_splicesites_file(gtf_file, align_dir, data):
    """
    if not pre-created, make a splicesites file to use with hisat2
    """
    out_file = os.path.join(align_dir, "ref-transcripts-splicesites.txt")
    if file_exists(out_file):
        return out_file
    safe_makedir(align_dir)
    hisat2_ss = config_utils.get_program("hisat2_extract_splice_sites.py", data)
    cmd = "{hisat2_ss} {gtf_file} > {tx_out_file}"
    message = "Creating hisat2 splicesites file from %s." % gtf_file
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    return out_file

def _get_quality_flag(data):
    qual_format = dd.get_quality_format(data)
    if qual_format.lower() == "illumina":
        return "--phred64"
    elif qual_format.lower() == "solexa":
        return "--solexa-quals"
    else:
        return "--phred33"

def _get_options_from_config(config):
    opts = []
    resources = config_utils.get_resources("hisat2", config)
    if resources.get("options"):
        opts += [str(x) for x in resources["options"]]
    return opts

def _get_stranded_flag(data, paired):
    strandedness = dd.get_strandedness(data)
    base = "--rna-strandness "
    if paired:
        if strandedness and strandedness.lower() in ["firstrand", "firststrand"]:
            return base + "RF"
        elif strandedness and strandedness.lower() == "secondstrand":
            return base + "FR"
        else:
            return ""
    else:
        if strandedness and strandedness.lower() in ["firstrand", "firststrand"]:
            return base + "R"
        elif strandedness and strandedness.lower() == "secondstrand":
            return base + "F"
        else:
            return ""

def _get_rg_flags(names):
    rg_id = names["rg"]
    rg_sample = names["sample"]
    rg_library = names["pl"]
    rg_platform_unit = names["pu"]
    rg_lb = ("--rg LB:%s " % names.get("lb")) if names.get("lb") else ""
    flags = ("--rg-id {rg_id} --rg PL:{rg_library} --rg PU:{rg_platform_unit} "
             "--rg SM:{rg_sample} {rg_lb}")
    return flags.format(**locals())

def remap_index_fn(ref_file):
    """Map sequence references to equivalent hisat2 indexes
    """
    return os.path.splitext(ref_file)[0].replace("/seq/", "/hisat2/")

def get_splicejunction_file(align_dir, data):
    """
    locate the splice junction file from hisat2. hisat2 outputs a novel
    splicesites file to go along with the provided file, if available.
    this combines the two together and outputs a combined file of all
    of the known and novel splice junctions
    """
    samplename = dd.get_sample_name(data)
    align_dir = os.path.dirname(dd.get_work_bam(data))
    knownfile = get_known_splicesites_file(align_dir, data)
    novelfile = os.path.join(align_dir, "%s-novelsplicesites.bed" % samplename)
    bed_files = [x for x in [knownfile, novelfile] if file_exists(x)]
    splicejunction = bed.concat(bed_files)
    splicejunctionfile = os.path.join(align_dir,
                                      "%s-splicejunctions.bed" % samplename)
    if splicejunction:
        splicejunction.saveas(splicejunctionfile)
        return splicejunctionfile
    else:
        return None
