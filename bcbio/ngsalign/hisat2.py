import os
from bcbio.utils import file_exists
import bcbio.pipeline.datadict as dd
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.ngsalign import postalign
from bcbio.provenance import do

def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    paired = True if pair_file else False
    hisat2 = config_utils.get_program("hisat2", data)
    num_cores = dd.get_num_cores(data)
    quality_flag = _get_quality_flag(data)
    stranded_flag = _get_stranded_flag(data, paired)
    rg_flags = _get_rg_flags(names)
    out_file = os.path.join(align_dir, dd.get_lane(data)) + ".bam"
    if file_exists(out_file):
        data = dd.set_work_bam(data, out_file)
        return data
    cmd = ("{hisat2} -x {ref_file} -p {num_cores} {quality_flag} {stranded_flag} "
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
        gtf_file = dd.get_gtf_file(data)
        splicesites = os.path.join(os.path.dirname(gtf_file),
                                   "ref-transcripts-splicesites.txt")
        cmd += "--known-splicesite-infile {splicesites} "
    message = "Aligning %s and %s with hisat2." %(fastq_file, pair_file)
    with file_transaction(out_file) as tx_out_file:
        cmd += " | " + postalign.sam_to_sortbam_cl(data, tx_out_file)
        do.run(cmd.format(**locals()), message)
    data = dd.set_work_bam(data, out_file)
    return data

def _get_quality_flag(data):
    qual_format = dd.get_quality_format(data)
    if qual_format.lower() == "illumina":
        return "--phred64"
    elif qual_format.lower() == "solexa":
        return "--solexa-quals"
    else:
        return "--phred33"

def _get_stranded_flag(data, paired):
    strandedness = dd.get_strandedness(data)
    base = "--rna-strandness "
    if paired:
        if strandedness == "firststrand":
            return base + "RF"
        elif strandedness == "secondstrand":
            return base + "FR"
        else:
            return ""
    else:
        if strandedness == "firstrand":
            return base + "R"
        elif strandedness == "secondstrand":
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
