"""Next-gen alignments with BWA (http://bio-bwa.sourceforge.net/)
"""
import os
import subprocess

from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import novoalign

galaxy_location_file = "bwa_index.loc"

def align_bam(in_bam, ref_file, names, align_dir, config):
    """Perform direct alignment of an input BAM file with BWA using pipes.

    This avoids disk IO by piping between processes:
     - samtools sort of input BAM to queryname
     - bedtools conversion to interleaved FASTQ
     - bwa-mem alignment
     - samtools conversion to BAM
     - samtools sort to coordinate
    """
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(names["lane"]))
    samtools = config_utils.get_program("samtools", config)
    bedtools = config_utils.get_program("bedtools", config)
    bwa = config_utils.get_program("bwa", config)
    resources = config_utils.get_resources("samtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    max_mem = resources.get("memory", "768M")
    rg_info = novoalign.get_rg_info(names)
    if not utils.file_exists(out_file):
        with utils.curdir_tmpdir() as work_dir:
            with file_transaction(out_file) as tx_out_file:
                tx_out_prefix = os.path.splitext(tx_out_file)[0]
                prefix1 = "%s-in1" % tx_out_prefix
                cmd = ("{samtools} sort -n -o -l 0 -@ {num_cores} -m {max_mem} {in_bam} {prefix1} "
                       "| {bedtools} bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout "
                       "| {bwa} mem -p -M -t {num_cores} -R '{rg_info}' -v 1 {ref_file} - "
                       "| {samtools} view -b -S -u - "
                       "| {samtools} sort -@ {num_cores} -m {max_mem} - {tx_out_prefix}")
                print cmd.format(**locals())
                subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          rg_name=None):
    """Perform a BWA alignment, generating a SAM file.
    """
    sai1_file = os.path.join(align_dir, "%s_1.sai" % out_base)
    sai2_file = (os.path.join(align_dir, "%s_2.sai" % out_base)
                 if pair_file else None)
    sam_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not utils.file_exists(sam_file):
        if not utils.file_exists(sai1_file):
            with file_transaction(sai1_file) as tx_sai1_file:
                _run_bwa_align(fastq_file, ref_file, tx_sai1_file, config)
        if sai2_file and not utils.file_exists(sai2_file):
            with file_transaction(sai2_file) as tx_sai2_file:
                _run_bwa_align(pair_file, ref_file, tx_sai2_file, config)
        align_type = "sampe" if sai2_file else "samse"
        sam_cl = [config_utils.get_program("bwa", config), align_type, ref_file, sai1_file]
        if sai2_file:
            sam_cl.append(sai2_file)
        sam_cl.append(fastq_file)
        if sai2_file:
            sam_cl.append(pair_file)
        with file_transaction(sam_file) as tx_sam_file:
            with open(tx_sam_file, "w") as out_handle:
                logger.info(" ".join(sam_cl))
                subprocess.check_call(sam_cl, stdout=out_handle)
    return sam_file

def _bwa_args_from_config(config):
    num_cores = config["algorithm"].get("num_cores", 1)
    core_flags = ["-t", str(num_cores)] if num_cores > 1 else []
    qual_format = config["algorithm"].get("quality_format", "").lower()
    qual_flags = ["-I"] if qual_format == "illumina" else []
    return core_flags + qual_flags

def _run_bwa_align(fastq_file, ref_file, out_file, config):
    aln_cl = [config_utils.get_program("bwa", config), "aln",
              "-n %s" % config["algorithm"]["max_errors"],
              "-k %s" % config["algorithm"]["max_errors"]]
    aln_cl += _bwa_args_from_config(config)
    aln_cl += [ref_file, fastq_file]
    with open(out_file, "w") as out_handle:
        logger.info(" ".join(aln_cl))
        subprocess.check_call(aln_cl, stdout=out_handle)
