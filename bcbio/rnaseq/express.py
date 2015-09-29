import os
import shutil
import tempfile
from bcbio.utils import file_exists
from bcbio.bam import is_paired, _get_sort_order, sort
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.rnaseq import gtf
from bcbio.log import logger

def run(data):
    """Quantitaive isoforms expression by eXpress"""
    name = dd.get_sample_name(data)
    in_bam = dd.get_transcriptome_bam(data)
    config = data['config']
    if not in_bam:
        logger.info("Transcriptome-mapped BAM file not found, skipping eXpress.")
        return data
    out_dir = os.path.join(dd.get_work_dir(data), "express", name)
    out_file = os.path.join(out_dir, name + ".xprs")
    express = config_utils.get_program("express", data['config'])
    strand = _set_stranded_flag(in_bam, data)
    if not file_exists(out_file):
        gtf_fasta = gtf.gtf_to_fasta(dd.get_gtf_file(data), dd.get_ref_file(data))
        with tx_tmpdir(data) as tmp_dir:
            with file_transaction(out_dir) as tx_out_dir:
                bam_file = _prepare_bam_file(in_bam, tmp_dir, config)
                cmd = ("{express} --no-update-check -o {tx_out_dir} {strand} {gtf_fasta} {bam_file}")
                do.run(cmd.format(**locals()), "Run express on %s." % in_bam, {})
            shutil.move(os.path.join(out_dir, "results.xprs"), out_file)
    eff_count_file = _get_column(out_file, out_file.replace(".xprs", "_eff.counts"), 7)
    tpm_file = _get_column(out_file, out_file.replace("xprs", "tpm"), 14)
    fpkm_file = _get_column(out_file, out_file.replace("xprs", "fpkm"), 10)
    data = dd.set_express_counts(data, eff_count_file)
    data = dd.set_express_tpm(data, tpm_file)
    data = dd.set_express_fpkm(data, fpkm_file)
    return data

def _get_column(in_file, out_file, column):
    """Subset one column from a file
    """
    with file_transaction(out_file) as tx_out_file:
        with open(in_file) as in_handle:
            with open(tx_out_file, 'w') as out_handle:
                for line in in_handle:
                    cols = line.strip().split("\t")
                    if line.find("eff_count") > 0:
                        continue
                    number = cols[column]
                    if column == 7:
                        number = int(round(float(number), 0))
                    out_handle.write("%s\t%s\n" % (cols[1], number))
    return out_file

def _set_stranded_flag(bam_file, data):
    strand_flag = {"unstranded": "",
                   "firststrand": "--rf-stranded",
                   "secondstrand": "--fr-stranded",
                   "firststrand-s": "--r-stranded",
                   "secondstrand-s": "--f-stranded"}
    stranded = dd.get_strandedness(data)
    assert stranded in strand_flag, ("%s is not a valid strandedness value. "
            "Valid values are 'firststrand', "
            "'secondstrand' and 'unstranded" % (stranded))
    if stranded != "unstranded" and not is_paired(bam_file):
        stranded += "-s"
    flag = strand_flag[stranded]
    return flag

def _prepare_bam_file(bam_file, tmp_dir, config):
    """
    Pipe sort by name cmd in case sort by coordinates
    """
    sort_mode = _get_sort_order(bam_file, config)
    if sort_mode != "queryname":
        bam_file = sort(bam_file, config, "queryname")
    return bam_file

def isoform_to_gene_name(gtf_file, out_file=None):
    """
    produce a table of isoform -> gene mappings for loading into EBSeq
    """
    if not out_file:
         out_file = tempfile.NamedTemporaryFile(delete=False).name
    if file_exists(out_file):
        return out_file
    db = gtf.get_gtf_db(gtf_file)
    line_format = "{transcript}\t{gene}\n"
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for feature in db.features_of_type('transcript'):
                transcript = feature['transcript_id'][0]
                gene = feature['gene_id'][0]
                out_handle.write(line_format.format(**locals()))
    return out_file
