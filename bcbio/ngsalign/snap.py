"""Alignment with SNAP: http://snap.cs.berkeley.edu/
"""
import os

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.ngsalign import novoalign, postalign
from bcbio.provenance import do

def align(fastq_file, pair_file, index_dir, names, align_dir, data):
    """Perform piped alignment of fastq input files, generating sorted, deduplicated BAM.

    TODO: Use streaming with new development version of SNAP to feed into
    structural variation preparation de-duplication.
    """
    pair_file = pair_file if pair_file else ""
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(names["lane"]))
    assert not data.get("align_split"), "Split alignments not supported with SNAP"
    snap = config_utils.get_program("snap", data["config"])
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    resources = config_utils.get_resources("snap", data["config"])
    rg_info = novoalign.get_rg_info(names)
    is_paired = bam.is_paired(fastq_file) if fastq_file.endswith(".bam") else pair_file
    if not utils.file_exists(out_file):
        with postalign.tobam_cl(data, out_file, is_paired) as (tobam_cl, tx_out_file):
            cmd_name = "paired" if is_paired else "single"
            cmd = ("{snap} {cmd_name} {index_dir} {fastq_file} {pair_file} "
                   "-R '{rg_info}' -t {num_cores} -M -o -sam - | ")
            do.run(cmd.format(**locals()) + tobam_cl, "SNAP alignment: %s" % names["sample"])
    data["work_bam"] = out_file
    return data

def align_bam(bam_file, index_dir, names, align_dir, data):
    return align(bam_file, None, index_dir, names, align_dir, data)

# Optional galaxy location file. Falls back on remap_index_fn if not found
galaxy_location_file = "snap_indices.loc"

def remap_index_fn(ref_file):
    """Map sequence references to snap reference directory, using standard layout.
    """
    snap_dir = os.path.join(os.path.dirname(ref_file), os.pardir, "snap")
    assert os.path.exists(snap_dir) and os.path.isdir(snap_dir), snap_dir
    return snap_dir
