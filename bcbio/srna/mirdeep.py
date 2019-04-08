import os
import sys
import os.path as op

import pysam

from bcbio.log import logger
from bcbio.utils import file_exists, safe_makedir, chdir, get_perl_exports
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd


def run(data):
    config = data[0][0]['config']
    work_dir = dd.get_work_dir(data[0][0])
    genome = dd.get_ref_file(data[0][0])
    mirdeep2 = os.path.join(os.path.dirname(sys.executable), "miRDeep2.pl")
    perl_exports = get_perl_exports()
    hairpin, mature, species = "none", "none", "na"
    rfam_file = dd.get_mirdeep2_file(data[0][0])
    if file_exists(dd.get_mirbase_hairpin(data[0][0])):
        species = dd.get_species(data[0][0])
        hairpin = dd.get_mirbase_hairpin(data[0][0])
        mature = dd.get_mirbase_mature(data[0][0])

    logger.debug("Preparing for mirdeep2 analysis.")
    bam_file = op.join(work_dir, "align", "seqs.bam")
    seqs_dir = op.join(work_dir, "seqcluster", "prepare")
    collapsed = op.join(seqs_dir, "seqs.ma")
    out_dir = op.join(work_dir, "mirdeep2")
    out_file = op.join(out_dir, "result_res.csv")
    safe_makedir(out_dir)
    if not file_exists(rfam_file):
        logger.warning("mirdeep2 Rfam file not instaled. Skipping...")
        return None
    if not file_exists(mirdeep2):
        logger.warning("mirdeep2 executable file not found. Skipping...")
        return None
    with chdir(out_dir):
        collapsed, bam_file = _prepare_inputs(collapsed, bam_file, out_dir)
        cmd = ("{perl_exports} && perl {mirdeep2} {collapsed} {genome} {bam_file} {mature} none {hairpin} -f {rfam_file} -r simple -c -P -t {species} -z res -g -1").format(**locals())
        if not file_exists(out_file):
            try:
                do.run(cmd.format(**locals()), "Running mirdeep2.")
            except:
                logger.warning("mirdeep2 failed. Please report the error to https://github.com/lpantano/mirdeep2_core/issues.")
        if file_exists(out_file):
            novel_db = _parse_novel(out_file, dd.get_species(data[0][0]))
            return novel_db

def _prepare_inputs(ma_fn, bam_file, out_dir):
    """
    Convert to fastq with counts
    """
    fixed_fa = os.path.join(out_dir, "file_reads.fa")
    count_name =dict()
    with file_transaction(fixed_fa) as out_tx:
        with open(out_tx, 'w') as out_handle:
            with open(ma_fn) as in_handle:
                h = next(in_handle)
                for line in in_handle:
                    cols = line.split("\t")
                    name_with_counts = "%s_x%s" % (cols[0], sum(map(int, cols[2:])))
                    count_name[cols[0]] = name_with_counts
                    out_handle.write(">%s\n%s\n" % (name_with_counts, cols[1]))
    fixed_bam = os.path.join(out_dir, "align.bam")
    bam_handle = pysam.AlignmentFile(bam_file, "rb")
    with pysam.AlignmentFile(fixed_bam, "wb", template=bam_handle) as out_handle:
        for read in bam_handle.fetch():
            read.query_name = count_name[read.query_name]
            out_handle.write(read)

    return fixed_fa, fixed_bam

def _parse_novel(csv_file, sps="new"):
    """Create input of novel miRNAs from miRDeep2"""
    read = 0
    seen = set()
    safe_makedir("novel")
    with open("novel/hairpin.fa", "w") as fa_handle, open("novel/miRNA.str", "w") as str_handle:
        with open(csv_file) as in_handle:
            for line in in_handle:
                if line.startswith("mature miRBase miRNAs detected by miRDeep2"):
                    break
                if line.startswith("novel miRNAs predicted"):
                    read = 1
                    line = next(in_handle)
                    continue
                if read and line.strip():
                    cols = line.strip().split("\t")
                    name, start, score = cols[0], cols[16], cols[1]
                    if float(score) < 1:
                        continue
                    m5p, m3p, pre = cols[13], cols[14], cols[15].replace('u', 't').upper()
                    m5p_start = cols[15].find(m5p) + 1
                    m3p_start = cols[15].find(m3p) + 1
                    m5p_end = m5p_start + len(m5p) - 1
                    m3p_end = m3p_start + len(m3p) - 1
                    if m5p in seen:
                        continue
                    fa_handle.write(">{sps}-{name} {start}\n{pre}\n".format(**locals()))
                    str_handle.write(">{sps}-{name} ({score}) [{sps}-{name}-5p:{m5p_start}-{m5p_end}] [{sps}-{name}-3p:{m3p_start}-{m3p_end}]\n".format(**locals()))
                    seen.add(m5p)
    return op.abspath("novel")
