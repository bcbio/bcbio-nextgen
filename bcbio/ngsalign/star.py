import os
import tempfile

from bcbio.pipeline import config_utils
from bcbio.utils import safe_makedir, file_exists, get_in, symlink_plus
from bcbio.provenance import do
from bcbio import bam

CLEANUP_FILES = ["Aligned.out.sam", "Log.out", "Log.progress.out"]

def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    config = data["config"]
    out_prefix = os.path.join(align_dir, names["lane"])
    out_file = out_prefix + "Aligned.out.sam"
    out_dir = os.path.join(align_dir, "%s_star" % names["lane"])

    final_out = os.path.join(out_dir, "{0}.bam".format(names["sample"]))
    if file_exists(final_out):
        return final_out
    star_path = config_utils.get_program("STAR", config)
    fastq = " ".join([fastq_file, pair_file]) if pair_file else fastq_file
    num_cores = config["algorithm"].get("num_cores", 1)

    safe_makedir(align_dir)
    cmd = ("{star_path} --genomeDir {ref_file} --readFilesIn {fastq} "
           "--runThreadN {num_cores} --outFileNamePrefix {out_prefix} "
           "--outReadsUnmapped Fastx --outFilterMultimapNmax 10 "
           "--outSAMunmapped Within")
    cmd += _read_group_option(names)
    fusion_mode = get_in(data, ("config", "algorithm", "fusion_mode"), False)
    if fusion_mode:
        cmd += " --chimSegmentMin 15 --chimJunctionOverhangMin 15"
    strandedness = get_in(data, ("config", "algorithm", "strandedness"),
                          "unstranded").lower()
    if strandedness == "unstranded":
        cmd += " --outSAMstrandField intronMotif"
    run_message = "Running STAR aligner on %s and %s." % (pair_file, ref_file)
    do.run(cmd.format(**locals()), run_message, None)
    out_file = bam.sam_to_bam(out_file, config)
    out_file = _fix_sam_header(out_file, config)
    if not file_exists(final_out):
        symlink_plus(out_file, final_out)
    return final_out

def _fix_sam_header(in_file, config):
    """
    STAR outputs a duplicate cl: line in the header which breaks some downstream
    tools like FastQC
    https://groups.google.com/d/msg/rna-star/xxE4cUnafJQ/EUsgYId-dB8J
    This can be safely removed whenever that bug gets fixed.
    """
    with bam.open_samfile(in_file) as in_handle:
        header = in_handle.header
    with tempfile.NamedTemporaryFile(delete=False) as header_handle:
        for key, line in header.items():
            line_key = "@" + str(key)
            for line_item in line:
                out_line = [line_key]
                out_line += [":".join([str(k), str(v)])
                             for k, v in line_item.items()
                             if k != "cl"]
                header_handle.write("\t".join(out_line) + "\n")
    header_name = header_handle.name
    header_handle.close()

    return bam.reheader(header_name, in_file, config)


def _read_group_option(names):
    rg_id = names["rg"]
    rg_sample = names["sample"]
    rg_library = names["pl"]
    rg_platform_unit = names["pu"]

    return (" --outSAMattrRGline ID:{rg_id} PL:{rg_library} "
            "PU:{rg_platform_unit} SM:{rg_sample} ").format(**locals())

def _get_quality_format(config):
    qual_format = config["algorithm"].get("quality_format", None)
    if qual_format.lower() == "illumina":
        return "fastq-illumina"
    elif qual_format.lower() == "solexa":
        return "fastq-solexa"
    else:
        return "fastq-sanger"

def remap_index_fn(ref_file):
    """Map sequence references to equivalent star indexes
    """
    return os.path.join(os.path.dirname(os.path.dirname(ref_file)), "star")
