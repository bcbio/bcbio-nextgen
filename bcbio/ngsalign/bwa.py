"""Next-gen alignments with BWA (http://bio-bwa.sourceforge.net/)
"""
import os
import subprocess

from bcbio.pipeline import config_utils
from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.ngsalign import alignprep, novoalign
from bcbio.provenance import do

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
    # adjust memory for samtools since used for input and output
    max_mem = config_utils.adjust_memory(resources.get("memory", "1G"),
                                         3, "decrease")
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
                cmd = cmd.format(**locals())
                do.run(cmd, "bwa mem alignment from BAM: %s" % names["sample"], None,
                       [do.file_nonempty(tx_out_file), do.file_reasonable_size(tx_out_file, in_bam)])
    return out_file

def can_pipe(fastq_file, data):
    """bwa-mem handle longer (> 70bp) reads with improved piping.
    Randomly samples 5000 reads from the first two million.
    Default to no piping if more than 75% of the sampled reads are small.
    """
    min_size = 70
    thresh = 0.75
    head_count = 8000000
    tocheck = 5000
    seqtk = config_utils.get_program("seqtk", data["config"])
    gzip_cmd = "zcat {fastq_file}" if fastq_file.endswith(".gz") else "cat {fastq_file}"
    cmd = (gzip_cmd + " | head -n {head_count} | "
           "{seqtk} sample -s42 - {tocheck} | "
           "awk '{{if(NR%4==2) print length($1)}}' | sort | uniq -c")
    count_out = subprocess.check_output(cmd.format(**locals()), shell=True,
                                        executable="/bin/bash")
    shorter = 0
    for count, size in (l.strip().split() for l in count_out.strip().split("\n")):
        if int(size) < min_size:
            shorter += int(count)
    return (float(shorter) / float(tocheck)) <= thresh

def align_pipe(fastq_file, pair_file, ref_file, names, align_dir, data):
    """Perform piped alignment of fastq input files, generating sorted output BAM.
    """
    pair_file = pair_file if pair_file else ""
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(names["lane"]))
    qual_format = data["config"]["algorithm"].get("quality_format", "").lower()
    if data.get("align_split"):
        final_file = out_file
        out_file, data = alignprep.setup_combine(final_file, data)
        fastq_file = alignprep.split_namedpipe_cl(fastq_file, data)
        if pair_file:
            pair_file = alignprep.split_namedpipe_cl(pair_file, data)
    else:
        final_file = None
        if qual_format == "illumina":
            fastq_file = alignprep.fastq_convert_pipe_cl(fastq_file, data)
            if pair_file:
                pair_file = alignprep.fastq_convert_pipe_cl(pair_file, data)
    samtools = config_utils.get_program("samtools", data["config"])
    bwa = config_utils.get_program("bwa", data["config"])
    resources = config_utils.get_resources("samtools", data["config"])
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    # adjust memory for samtools since used alongside alignment
    max_mem = config_utils.adjust_memory(resources.get("memory", "2G"),
                                         3, "decrease")
    rg_info = novoalign.get_rg_info(names)
    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        # If we cannot do piping, use older bwa aln approach
        if not can_pipe(fastq_file, data):
            return align(fastq_file, pair_file, ref_file, names, align_dir, data)
        else:
            with utils.curdir_tmpdir() as work_dir:
                with file_transaction(out_file) as tx_out_file:
                    tx_out_prefix = os.path.splitext(tx_out_file)[0]
                    cmd = ("{bwa} mem -M -t {num_cores} -R '{rg_info}' -v 1 {ref_file} "
                           "{fastq_file} {pair_file} "
                           "| {samtools} view -b -S -u - "
                           "| {samtools} sort -@ {num_cores} -m {max_mem} - {tx_out_prefix}")
                    cmd = cmd.format(**locals())
                    do.run(cmd, "bwa mem alignment from fastq: %s" % names["sample"], None,
                           [do.file_nonempty(tx_out_file), do.file_reasonable_size(tx_out_file, fastq_file)])
    data["work_bam"] = out_file
    return data

def align(fastq_file, pair_file, ref_file, names, align_dir, data):
    """Perform a BWA alignment, generating a SAM file.
    """
    assert not data.get("align_split"), "Do not handle split alignments with non-piped bwa"
    config = data["config"]
    sai1_file = os.path.join(align_dir, "%s_1.sai" % names["lane"])
    sai2_file = (os.path.join(align_dir, "%s_2.sai" % names["lane"])
                 if pair_file else None)
    sam_file = os.path.join(align_dir, "%s.sam" % names["lane"])
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
            cmd = "{cl} > {out_file}".format(cl=" ".join(sam_cl), out_file=tx_sam_file)
            do.run(cmd, "bwa {align_type}".format(**locals()), None)
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
    cmd = "{cl} > {out_file}".format(cl=" ".join(aln_cl), out_file=out_file)
    do.run(cmd, "bwa aln: {f}".format(f=os.path.basename(fastq_file)), None)
