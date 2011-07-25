"""Pipeline support for barcode analysis and de-mulitplexing.
"""
import os
import copy
import subprocess

from Bio import SeqIO

from bcbio import utils
from bcbio.pipeline.fastq import get_fastq_files

def split_by_barcode(fastq1, fastq2, multiplex, base_name, dirs, config):
    """Split a fastq file into multiplex pieces using barcode details.
    """
    if not multiplex:
        return [("", "", fastq1, fastq2)]
    bc_dir = os.path.join(dirs["work"], "%s_barcode" % base_name)
    nomatch_file = "%s_unmatched_1_fastq.txt" % base_name
    metrics_file = "%s_bc.metrics" % base_name
    with utils.chdir(bc_dir):
        if not os.path.exists(nomatch_file) and not os.path.exists(metrics_file):
            tag_file = _make_tag_file(multiplex)
            cl = [config["program"]["barcode"], tag_file,
                  "%s_--b--_--r--_fastq.txt" % base_name,
                  fastq1]
            if fastq2:
                cl.append(fastq2)
            cl.append("--mismatch=%s" % config["algorithm"]["bc_mismatch"])
            cl.append("--metrics=%s" % metrics_file)
            if int(config["algorithm"]["bc_read"]) == 2:
                cl.append("--second")
            if int(config["algorithm"]["bc_position"]) == 5:
                cl.append("--five")
            if config["algorithm"].get("bc_allow_indels", True) is False:
                cl.append("--noindel")
            subprocess.check_call(cl)
    out_files = []
    for info in multiplex:
        fq_fname = lambda x: os.path.join(bc_dir, "%s_%s_%s_fastq.txt" %
                             (base_name, info["barcode_id"], x))
        bc_file1 = fq_fname("1")
        bc_file2 = fq_fname("2") if fastq2 else None
        out_files.append((info["barcode_id"], info["name"], bc_file1, bc_file2))
    return out_files

def _make_tag_file(barcodes):
    tag_file = "%s-barcodes.cfg" % barcodes[0].get("barcode_type", "barcode")
    barcodes = _adjust_illumina_tags(barcodes)
    with open(tag_file, "w") as out_handle:
        for bc in barcodes:
            out_handle.write("%s %s\n" % (bc["barcode_id"], bc["sequence"]))
    return tag_file

def _adjust_illumina_tags(barcodes):
    """Handle additional trailing A in Illumina barocdes.

    Illumina barcodes are listed as 6bp sequences but have an additional
    A base when coming off on the sequencer. This checks for this case and
    adjusts the sequences appropriately if needed.
    """
    illumina_size = 7
    all_illumina = True
    need_a = False
    for bc in barcodes:
        if bc.get("barcode_type", "illumina").lower().find("illumina") == -1:
            all_illumina = False
        if (not bc["sequence"].upper().endswith("A") or
            len(bc["sequence"]) < illumina_size):
            need_a = True
    if all_illumina and need_a:
        new = []
        for bc in barcodes:
            new_bc = copy.deepcopy(bc)
            new_bc["sequence"] = "%sA" % new_bc["sequence"]
            new.append(new_bc)
        barcodes = new
    return barcodes

def add_multiplex_across_lanes(run_items, fastq_dir, fc_name):
    """Add multiplex information to control and non-multiplexed lanes.

    Illumina runs include barcode reads for non-multiplex lanes, and the
    control, when run on a multiplexed flow cell. This checks for this
    situation and adds details to trim off the extra bases.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    # determine if we have multiplexes and collect expected size
    fastq_sizes = []
    tag_sizes = []
    has_barcodes = False
    for item in run_items:
        if item.get("multiplex", None):
            has_barcodes = True
            tag_sizes.extend([len(b["sequence"]) for b in item["multiplex"]])
            fastq_sizes.append(_get_fastq_size(item, fastq_dir, fc_name))
    if not has_barcodes: # nothing to worry about
        return run_items
    fastq_sizes = list(set(fastq_sizes))

    # discard 0 sizes to handle the case where lane(s) are empty or failed
    try:
        fastq_sizes.remove(0)
    except ValueError: pass

    tag_sizes = list(set(tag_sizes))
    final_items = []
    for item in run_items:
        if item.get("multiplex", None) is None:
            assert len(fastq_sizes) == 1, \
                   "Multi and non-multiplex reads with multiple sizes"
            expected_size = fastq_sizes[0]
            assert len(tag_sizes) == 1, \
                   "Expect identical tag size for a flowcell"
            tag_size = tag_sizes[0]
            this_size = _get_fastq_size(item, fastq_dir, fc_name)
            if this_size == expected_size:
                item["multiplex"] = [{"name" : item.get("name", item["description"]),
                                      "barcode_id": "trim",
                                      "sequence" : "N" * tag_size}]
            else:
                assert this_size == expected_size - tag_size, \
                       "Unexpected non-multiplex sequence"
        final_items.append(item)
    return final_items

def _get_fastq_size(item, fastq_dir, fc_name):
    """Retrieve the size of reads from the first flowcell sequence.
    """
    (fastq1, _) = get_fastq_files(fastq_dir, item, fc_name)
    with open(fastq1) as in_handle:
        try:
            rec = SeqIO.parse(in_handle, "fastq").next()
            size = len(rec.seq)
        except StopIteration:
            size = 0
    return size

