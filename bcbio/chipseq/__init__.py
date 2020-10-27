import os
import shutil
import subprocess
import sys
import toolz as tz
from bcbio import utils
from bcbio import bam
from bcbio.pipeline import config_utils
import bcbio.pipeline.datadict as dd
from bcbio.ngsalign import bowtie2, bwa
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.log import logger
from bcbio.heterogeneity.chromhacks import (get_mitochondrial_chroms,
                                            get_nonmitochondrial_chroms)
from bcbio.chipseq import atac

def clean_chipseq_alignment(data):
    # lcr_bed = utils.get_in(data, ("genome_resources", "variation", "lcr"))
    method = dd.get_chip_method(data)
    if method == "atac":
        data = shift_ATAC(data)
    work_bam = dd.get_work_bam(data)
    work_bam = bam.sort(work_bam, dd.get_config(data))
    bam.index(work_bam, dd.get_config(data))
    # an unfiltered BAM file is useful for calculating some metrics later
    data = tz.assoc_in(data, ['chipseq', 'align', "unfiltered"],  work_bam)
    clean_bam = remove_nonassembled_chrom(work_bam, data)
    clean_bam = remove_mitochondrial_reads(clean_bam, data)
    data = atac.calculate_complexity_metrics(clean_bam, data)
    if not dd.get_keep_multimapped(data):
        clean_bam = remove_multimappers(clean_bam, data)
    if not dd.get_keep_duplicates(data):
        clean_bam = bam.remove_duplicates(clean_bam, data)
    data["work_bam"] = clean_bam
    # for ATAC-seq, brewak alignments into NF, mono/di/tri nucleosome BAM files
    if method == "atac":
        data  = atac.split_ATAC(data)
    encode_bed = tz.get_in(["genome_resources", "variation", "encode_blacklist"], data)
    if encode_bed:
        data["work_bam"] = remove_blacklist_regions(dd.get_work_bam(data), encode_bed, data['config'])
        bam.index(data["work_bam"], data['config'])
    try:
        data["bigwig"] = _normalized_bam_coverage(dd.get_sample_name(data),
                                                  dd.get_work_bam(data), data)
    except subprocess.CalledProcessError:
        logger.warning(f"{dd.get_work_bam(data)} was too sparse to normalize, "
                       f" falling back to non-normalized coverage.")
        data["bigwig"] = _bam_coverage(dd.get_sample_name(data),
                                       dd.get_work_bam(data), data)
    return [[data]]

def remove_mitochondrial_reads(bam_file, data):
    mito = get_mitochondrial_chroms(data)
    if not mito:
        logger.info(f"Mitochondrial chromosome not identified, skipping removal of "
                    "mitochondrial reads from {bam_file}.")
        return bam_file
    nonmito = get_nonmitochondrial_chroms(data)
    mito_bam = os.path.splitext(bam_file)[0] + "-noMito.bam"
    if utils.file_exists(mito_bam):
        return mito_bam
    samtools = config_utils.get_program("samtools", dd.get_config(data))
    nonmito_flag = " ".join(nonmito)
    num_cores = dd.get_num_cores(data)
    with file_transaction(mito_bam) as tx_out_bam:
        cmd = (f"{samtools} view -bh -@ {num_cores} {bam_file} {nonmito_flag} "
               f"> {tx_out_bam}")
        message = f"Removing mitochondrial reads on {','.join(mito)} from {bam_file}."
        do.run(cmd, message)
    return mito_bam

def remove_multimappers(bam_file, data):
    aligner = dd.get_aligner(data)
    if aligner:
        if aligner == "bowtie2":
            filterer = bowtie2.filter_multimappers
        elif aligner == "bwa":
            filterer = bwa.filter_multimappers
        else:
            logger.error("ChIP-seq only supported for bowtie2 and bwa.")
            sys.exit(-1)
        unique_bam = filterer(bam_file, data)
    else:
        unique_bam = bam_file
        logger.warn("When a BAM file is given as input, bcbio skips removal of "
                    "multimappers.")
    return unique_bam

def remove_nonassembled_chrom(bam_file, data):
    """Remove non-assembled contigs from the BAM file"""
    ref_file =  dd.get_ref_file(data)
    config = dd.get_config(data)
    bam.index(bam_file, config)
    fai = "%s.fai" % ref_file
    chrom = []
    with open(fai) as inh:
        for line in inh:
            c = line.split("\t")[0]
            if c.find("_") < 0:
                chrom.append(c)
    chroms = " ".join(chrom)
    out_file = utils.append_stem(bam_file, '_chrom')
    samtools = config_utils.get_program("samtools", config)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out:
            cmd = "{samtools} view -b {bam_file} {chroms} > {tx_out}"
            do.run(cmd.format(**locals()), "Remove contigs from %s" % bam_file)
        bam.index(out_file, config)
    return out_file

def remove_blacklist_regions(bam_file, bed_file, config):
    """Remove blacklist regions from a BAM file"""
    if not bam_file or not bed_file:
        return bam_file
    out_file = utils.append_stem(bam_file, '_filter')
    bedtools = config_utils.get_program("bedtools", config)
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out:
            cmd = "{bedtools} subtract -nonamecheck -A -a {bam_file} -b {bed_file} > {tx_out}"
            do.run(cmd.format(**locals()), "Remove blacklist regions from %s" % bam_file)
    return out_file

def _bam_coverage(name, bam_input, data):
    """Run bamCoverage from deeptools"""
    cmd = ("{bam_coverage} --bam {bam_input} --outFileName {bw_output} "
          "--binSize 20 --effectiveGenomeSize {size} "
          "--smoothLength 60 --extendReads 150 --centerReads -p {cores} ")
    size = bam.fasta.total_sequence_length(dd.get_ref_file(data))
    cores = dd.get_num_cores(data)
    try:
        bam_coverage = config_utils.get_program("bamCoverage", data)
    except config_utils.CmdNotFound:
        logger.info("No bamCoverage found, skipping bamCoverage.")
        return None
    resources = config_utils.get_resources("bamCoverage", data["config"])
    if resources:
        options = resources.get("options")
        if options:
            cmd += " %s" % " ".join([str(x) for x in options])
    bw_output = os.path.join(os.path.dirname(bam_input), "%s.bw" % name)
    if utils.file_exists(bw_output):
        return bw_output
    with file_transaction(bw_output) as out_tx:
        do.run(cmd.format(**locals()), "Run bamCoverage in %s" % name)
    return bw_output

def _normalized_bam_coverage(name, bam_input, data):
    """Run bamCoverage from deeptools but produce normalized bigWig files"""
    cmd = ("{bam_coverage} --bam {bam_input} --outFileName {bw_output} "
          "--binSize 20 --effectiveGenomeSize {size} "
          "--smoothLength 60 --extendReads 150 --centerReads -p {cores} ")
    size = bam.fasta.total_sequence_length(dd.get_ref_file(data))
    cores = dd.get_num_cores(data)
    try:
        bam_coverage = config_utils.get_program("bamCoverage", data)
    except config_utils.CmdNotFound:
        logger.info("No bamCoverage found, skipping bamCoverage.")
        return None
    method = dd.get_chip_method(data)
    cmd += "--normalizeUsing CPM "
    toignore = get_mitochondrial_chroms(data)
    if toignore:
        ignorenormflag = f"--ignoreForNormalization {' '.join(toignore)} "
        cmd += ignorenormflag
    resources = config_utils.get_resources("bamCoverage", data["config"])
    if resources:
        options = resources.get("options")
        if options:
            cmd += " %s" % " ".join([str(x) for x in options])
    bw_output = os.path.join(os.path.dirname(bam_input), "%s.bw" % name)
    if utils.file_exists(bw_output):
        return bw_output
    with file_transaction(bw_output) as out_tx:
        do.run(cmd.format(**locals()), "Run bamCoverage in %s" % name)
    return bw_output

def shift_ATAC(data):
    """
    shift the ATAC-seq alignments
    """
    MAX_FRAG_LENGTH = 100
    sieve = config_utils.get_program("alignmentSieve", data)
    work_bam = dd.get_work_bam(data)
    num_cores = dd.get_num_cores(data)
    out_file = os.path.splitext(work_bam)[0] + "-shifted.bam"
    log_file = os.path.splitext(work_bam)[0] + "-shifted.log"
    if utils.file_exists(out_file):
        data["work_bam"] = out_file
        return data

    unsorted_bam = os.path.splitext(out_file)[0] + ".unsorted.bam"
    # shifting removes all reads if the BAM file is not paired
    shiftflag = "--ATACshift" if bam.is_paired(work_bam) else ""
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file, \
            file_transaction(log_file) as tx_log_file:
            tx_unsorted_file = os.path.splitext(tx_out_file)[0] + ".tmp.bam"
            cmd = (
                f"{sieve} --verbose --bam {work_bam} --outFile {tx_unsorted_file} "
                f"{shiftflag} "
                f"--numberOfProcessors {num_cores} --maxFragmentLength 0 "
                f"--minFragmentLength 0 "
                f"--minMappingQuality 10 "
                f"--filterMetrics {tx_log_file} ")
            do.run(cmd, f"Shifting ATAC-seq alignments in {work_bam} to {tx_unsorted_file}.")
            # shifting can cause the file to become unsorted
            sorted_file = bam.sort(tx_unsorted_file, dd.get_config(data), force=True)
            shutil.move(sorted_file, tx_out_file)
    bam.index(out_file, dd.get_config(data))
    data["work_bam"] = out_file
    return data
