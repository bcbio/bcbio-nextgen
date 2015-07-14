"""Next-gen alignments with BWA (http://bio-bwa.sourceforge.net/)
"""
import os
import subprocess

import toolz as tz

from bcbio.pipeline import config_utils
from bcbio import bam, utils
from bcbio.distributed import objectstore
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.ngsalign import alignprep, novoalign, postalign
from bcbio.provenance import do
from bcbio.rnaseq import gtf
import bcbio.pipeline.datadict as dd
from bcbio.bam import fastq
from bcbio.log import logger

galaxy_location_file = "bwa_index.loc"

def align_bam(in_bam, ref_file, names, align_dir, data):
    """Perform direct alignment of an input BAM file with BWA using pipes.

    This avoids disk IO by piping between processes:
     - samtools sort of input BAM to queryname
     - bedtools conversion to interleaved FASTQ
     - bwa-mem alignment
     - samtools conversion to BAM
     - samtools sort to coordinate
    """
    config = data["config"]
    out_file = os.path.join(align_dir, "{0}-sort.bam".format(names["lane"]))
    samtools = config_utils.get_program("samtools", config)
    bedtools = config_utils.get_program("bedtools", config)
    resources = config_utils.get_resources("samtools", config)
    num_cores = config["algorithm"].get("num_cores", 1)
    # adjust memory for samtools since used for input and output
    max_mem = config_utils.adjust_memory(resources.get("memory", "1G"),
                                         3, "decrease").upper()
    if not utils.file_exists(out_file):
        with tx_tmpdir(data) as work_dir:
            with postalign.tobam_cl(data, out_file, bam.is_paired(in_bam)) as (tobam_cl, tx_out_file):
                bwa_cmd = _get_bwa_mem_cmd(data, out_file, ref_file, "-")
                tx_out_prefix = os.path.splitext(tx_out_file)[0]
                prefix1 = "%s-in1" % tx_out_prefix
                cmd = ("{samtools} sort -n -o -l 1 -@ {num_cores} -m {max_mem} {in_bam} {prefix1} "
                       "| {bedtools} bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout "
                       "| {bwa_cmd} | ")
                cmd = cmd.format(**locals()) + tobam_cl
                do.run(cmd, "bwa mem alignment from BAM: %s" % names["sample"], None,
                       [do.file_nonempty(tx_out_file), do.file_reasonable_size(tx_out_file, in_bam)])
    return out_file

def _get_bwa_mem_cmd(data, out_file, ref_file, fastq1, fastq2=""):
    """Perform piped bwa mem mapping potentially with alternative alleles in GRCh38 + HLA typing.

    Commands for HLA post-processing:
       base=TEST
       run-HLA $base.hla > $base.hla.top
       cat $base.hla.HLA*.gt | grep ^GT | cut -f2- > $base.hla.all
       rm -f $base.hla.HLA*gt
       rm -f $base.hla.HLA*gz
    """
    alt_file = ref_file + ".alt"
    if utils.file_exists(alt_file):
        bwakit_dir = os.path.dirname(os.path.realpath(utils.which("run-bwamem")))
        hla_base = os.path.join(utils.safe_makedir(os.path.join(os.path.dirname(out_file), "hla")),
                                os.path.basename(out_file) + ".hla")
        alt_cmd = (" | {bwakit_dir}/k8 {bwakit_dir}/bwa-postalt.js -p {hla_base} {alt_file}")
    else:
        alt_cmd = ""
    bwa = config_utils.get_program("bwa", data["config"])
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    bwa_resources = config_utils.get_resources("bwa", data["config"])
    bwa_params = (" ".join([str(x) for x in bwa_resources.get("options", [])])
                  if "options" in bwa_resources else "")
    rg_info = novoalign.get_rg_info(data["rgnames"])
    pairing = "-p" if not fastq2 else ""
    # Restrict seed occurances to 1/2 of default, manage memory usage for centromere repeats in hg38
    # https://sourceforge.net/p/bio-bwa/mailman/message/31514937/
    # http://ehc.ac/p/bio-bwa/mailman/message/32268544/
    mem_usage = "-c 250"
    bwa_cmd = ("{bwa} mem {pairing} {mem_usage} -M -t {num_cores} {bwa_params} -R '{rg_info}' -v 1 "
               "{ref_file} {fastq1} {fastq2} ")
    return (bwa_cmd + alt_cmd).format(**locals())

def _can_use_mem(fastq_file, data):
    """bwa-mem handle longer (> 70bp) reads with improved piping.
    Randomly samples 5000 reads from the first two million.
    Default to no piping if more than 75% of the sampled reads are small.
    """
    min_size = 70
    thresh = 0.75
    head_count = 8000000
    tocheck = 5000
    seqtk = config_utils.get_program("seqtk", data["config"])
    fastq_file = objectstore.cl_input(fastq_file)
    gzip_cmd = "zcat {fastq_file}" if fastq_file.endswith(".gz") else "cat {fastq_file}"
    cmd = (gzip_cmd + " | head -n {head_count} | "
           "{seqtk} sample -s42 - {tocheck} | "
           "awk '{{if(NR%4==2) print length($1)}}' | sort | uniq -c")
    count_out = subprocess.check_output(cmd.format(**locals()), shell=True,
                                        executable="/bin/bash", stderr=open("/dev/null", "w"))
    if not count_out.strip():
        raise IOError("Failed to check fastq file sizes with: %s" % cmd.format(**locals()))
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
    rg_info = novoalign.get_rg_info(names)
    if not utils.file_exists(out_file) and (final_file is None or not utils.file_exists(final_file)):
        # If we cannot do piping, use older bwa aln approach
        if ("bwa-mem" in tz.get_in(["config", "algorithm", "tools_off"], data, [])
              or not _can_use_mem(fastq_file, data)):
            out_file = _align_backtrack(fastq_file, pair_file, ref_file, out_file,
                                        names, rg_info, data)
        else:
            out_file = _align_mem(fastq_file, pair_file, ref_file, out_file,
                                  names, rg_info, data)
    data["work_bam"] = out_file
    return data

def _align_mem(fastq_file, pair_file, ref_file, out_file, names, rg_info, data):
    """Perform bwa-mem alignment on supported read lengths.
    """
    with postalign.tobam_cl(data, out_file, pair_file != "") as (tobam_cl, tx_out_file):
        cmd = "%s | %s" % (_get_bwa_mem_cmd(data, out_file, ref_file, fastq_file, pair_file), tobam_cl)
        do.run(cmd, "bwa mem alignment from fastq: %s" % names["sample"], None,
                [do.file_nonempty(tx_out_file), do.file_reasonable_size(tx_out_file, fastq_file)])
    return out_file

def _align_backtrack(fastq_file, pair_file, ref_file, out_file, names, rg_info, data):
    """Perform a BWA alignment using 'aln' backtrack algorithm.
    """
    assert not data.get("align_split"), "Do not handle split alignments with non-piped bwa"
    bwa = config_utils.get_program("bwa", data["config"])
    config = data["config"]
    sai1_file = "%s_1.sai" % os.path.splitext(out_file)[0]
    sai2_file = "%s_2.sai" % os.path.splitext(out_file)[0] if pair_file else ""
    if not utils.file_exists(sai1_file):
        with file_transaction(data, sai1_file) as tx_sai1_file:
            _run_bwa_align(fastq_file, ref_file, tx_sai1_file, config)
    if sai2_file and not utils.file_exists(sai2_file):
        with file_transaction(data, sai2_file) as tx_sai2_file:
            _run_bwa_align(pair_file, ref_file, tx_sai2_file, config)
    with postalign.tobam_cl(data, out_file, pair_file != "") as (tobam_cl, tx_out_file):
        align_type = "sampe" if sai2_file else "samse"
        cmd = ("{bwa} {align_type} -r '{rg_info}' {ref_file} {sai1_file} {sai2_file} "
               "{fastq_file} {pair_file} | ")
        cmd = cmd.format(**locals()) + tobam_cl
        do.run(cmd, "bwa %s" % align_type, data)
    return out_file

def _bwa_args_from_config(config):
    num_cores = config["algorithm"].get("num_cores", 1)
    core_flags = ["-t", str(num_cores)] if num_cores > 1 else []
    return core_flags

def _run_bwa_align(fastq_file, ref_file, out_file, config):
    aln_cl = [config_utils.get_program("bwa", config), "aln",
              "-n 2", "-k 2"]
    aln_cl += _bwa_args_from_config(config)
    aln_cl += [ref_file, fastq_file]
    cmd = "{cl} > {out_file}".format(cl=" ".join(aln_cl), out_file=out_file)
    do.run(cmd, "bwa aln: {f}".format(f=os.path.basename(fastq_file)), None)

def index_transcriptome(gtf_file, ref_file, data):
    """
    use a GTF file and a reference FASTA file to index the transcriptome
    """
    gtf_fasta = gtf.gtf_to_fasta(gtf_file, ref_file)
    bwa = config_utils.get_program("bwa", data["config"])
    cmd = "{bwa} index {gtf_fasta}".format(**locals())
    message = "Creating transcriptome index of %s with bwa." % (gtf_fasta)
    do.run(cmd, message)
    return gtf_fasta

def align_transcriptome(fastq_file, pair_file, ref_file, data):
    """
    bwa mem with settings for aligning to the transcriptome for eXpress/RSEM/etc
    """
    work_bam = dd.get_work_bam(data)
    base, ext = os.path.splitext(work_bam)
    out_file = base + ".transcriptome" + ext
    if utils.file_exists(out_file):
        data = dd.set_transcriptome_bam(data, out_file)
        return data
    # bwa mem needs phred+33 quality, so convert if it is Illumina
    if dd.get_quality_format(data).lower() == "illumina":
        logger.info("bwa mem does not support the phred+64 quality format, "
                    "converting %s and %s to phred+33.")
        fastq_file = fastq.groom(fastq_file, in_qual="fastq-illumina", data=data)
        if pair_file:
            pair_file = fastq.groom(pair_file, in_qual="fastq-illumina", data=data)
    bwa = config_utils.get_program("bwa", data["config"])
    gtf_file = dd.get_gtf_file(data)
    gtf_fasta = index_transcriptome(gtf_file, ref_file, data)
    args = " ".join(_bwa_args_from_config(data["config"]))
    num_cores = data["config"]["algorithm"].get("num_cores", 1)
    cmd = ("{bwa} mem {args} -a -t {num_cores} {gtf_fasta} {fastq_file} "
           "{pair_file} | samtools view -bhS - > {tx_out_file}")

    with file_transaction(out_file) as tx_out_file:
        message = "Aligning %s and %s to the transcriptome." % (fastq_file, pair_file)
        do.run(cmd.format(**locals()), message)
    data = dd.set_transcriptome_bam(data, out_file)
    return data
