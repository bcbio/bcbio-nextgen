"""Alignment with minimap2: https://github.com/lh3/minimap2
"""
import os

from bcbio import utils
from bcbio.pipeline import datadict as dd
from bcbio.ngsalign import alignprep, novoalign, postalign
from bcbio.provenance import do

def align(fastq_file, pair_file, index_dir, names, align_dir, data):
    """Perform piped alignment of fastq input files, generating sorted, deduplicated BAM.
    """
    umi_ext = "-cumi" if "umi_bam" in data else ""
    out_file = os.path.join(align_dir, "{0}-sort{1}.bam".format(dd.get_sample_name(data), umi_ext))
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    rg_info = novoalign.get_rg_info(names)
    preset = "sr"

    pair_file = pair_file if pair_file else ""
    if data.get("align_split"):
        final_file = out_file
        out_file, data = alignprep.setup_combine(final_file, data)
        fastq_file, pair_file = alignprep.split_namedpipe_cls(fastq_file, pair_file, data)
    else:
        final_file = None

    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        with postalign.tobam_cl(data, out_file, pair_file != "") as (tobam_cl, tx_out_file):
            index_file = None
            # Skip trying to use indices now as they provide only slight speed-ups
            # and give inconsitent outputs in BAM headers
            # If a single index present, index_dir points to that
            # if index_dir and os.path.isfile(index_dir):
            #     index_dir = os.path.dirname(index_dir)
            #     index_file = os.path.join(index_dir, "%s-%s.mmi" % (dd.get_genome_build(data), preset))
            if not index_file or not os.path.exists(index_file):
                index_file = dd.get_ref_file(data)
            cmd = ("minimap2 -a -x {preset} -R '{rg_info}' -t {num_cores} {index_file} "
                   "{fastq_file} {pair_file} | ")
            do.run(cmd.format(**locals()) + tobam_cl, "minimap2 alignment: %s" % dd.get_sample_name(data))
    data["work_bam"] = out_file
    return data

def remap_index_fn(ref_file):
    """minimap2 can build indexes on the fly but will also store commons ones.
    """
    index_dir = os.path.join(os.path.dirname(ref_file), os.pardir, "minimap2")
    if os.path.exists(index_dir) and os.path.isdir(index_dir):
        return index_dir
    else:
        return os.path.dirname(ref_file)
