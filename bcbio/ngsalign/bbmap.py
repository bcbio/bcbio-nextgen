"""Alignment with bbmap: https://sourceforge.net/projects/bbmap/
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
    rg_info = "rgid={rg} rgpl={pl} rgpu={pu} rgsm={sample}".format(**names)
    pair_file = pair_file if pair_file else ""
    final_file = None
    if data.get("align_split"):
        # BBMap does not accept input fastq streams
        raise ValueError("bbmap is not compatible with alignment splitting, set `align_split: false`")
    pair_arg = "in2=%s" % pair_file if pair_file else ""
    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        with postalign.tobam_cl(data, out_file, pair_file != "") as (tobam_cl, tx_out_file):
            if index_dir.endswith(("/ref", "/ref/")):
                index_dir = os.path.dirname(index_dir)
            # sam=1.3 required for compatibility with strelka2
            cmd = ("bbmap.sh sam=1.3 mdtag=t {rg_info} path={index_dir} in1={fastq_file} "
                   "{pair_arg} out=stdout.sam | ")
            do.run(cmd.format(**locals()) + tobam_cl, "bbmap alignment: %s" % dd.get_sample_name(data))
    data["work_bam"] = out_file
    return data

def remap_index_fn(ref_file):
    index_dir = os.path.join(os.path.dirname(ref_file), os.pardir, "bbmap")
    if os.path.exists(index_dir) and os.path.isdir(index_dir):
        return index_dir
    else:
        return os.path.dirname(ref_file)
