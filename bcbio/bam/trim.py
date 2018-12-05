"""atropos cutadapt-like multicore trimming of input reads from Fastq or BAM files.
"""
import os
import sys

from Bio.Seq import Seq

from bcbio import utils
from bcbio.bam.fastq import is_fastq
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.distributed import objectstore
from bcbio.pipeline import config_utils
from bcbio.provenance import do
import bcbio.pipeline.datadict as dd

SUPPORTED_ADAPTERS = {
    "illumina": ["AACACTCTTTCCCT", "AGATCGGAAGAGCG"],
    "truseq": ["AGATCGGAAGAG"],
    "polya": ["AAAAAAAAAAAAA"],
    "nextera": ["AATGATACGGCGA", "CAAGCAGAAGACG"],
    "truseq2": ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"], # 3'only: first read, second read
    "nextera2": ["CTGTCTCTTATACACATCT", "AGATGTGTATAAGAGACAG"] # Second read in pair 3', 5
}

def trim_adapters(data):
    to_trim = [x for x in data["files"] if x is not None and is_fastq(x)]
    if not to_trim:
        return data["files"]

    logger.info("Trimming low quality ends and read through adapter "
                "sequence from %s." % (", ".join(to_trim)))
    out_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(data), "trimmed"))
    return _trim_adapters(to_trim, out_dir, data)

def _trim_adapters(fastq_files, out_dir, data):
    """
    for small insert sizes, the read length can be longer than the insert
    resulting in the reverse complement of the 3' adapter being sequenced.
    this takes adapter sequences and trims the only the reverse complement
    of the adapter

    MYSEQUENCEAAAARETPADA -> MYSEQUENCEAAAA (no polyA trim)
    """
    to_trim = _get_sequences_to_trim(data["config"], SUPPORTED_ADAPTERS)
    if dd.get_trim_reads(data) == "fastp":
        out_files, report_file = _fastp_trim(fastq_files, to_trim, out_dir, data)
    else:
        out_files, report_file = _atropos_trim(fastq_files, to_trim, out_dir, data)
    # quality_format = _get_quality_format(data["config"])
    # out_files = replace_directory(append_stem(fastq_files, "_%s.trimmed" % name), out_dir)
    # log_file = "%s_log_cutadapt.txt" % splitext_plus(out_files[0])[0]
    # out_files = _cutadapt_trim(fastq_files, quality_format, to_trim, out_files, log_file, data)
    # if file_exists(log_file):
    #     content = open(log_file).read().replace(fastq_files[0], name)
    #     if len(fastq_files) > 1:
    #         content = content.replace(fastq_files[1], name)
    #     open(log_file, 'w').write(content)
    return out_files

# ## Atropos trimming

def _atropos_trim(fastq_files, adapters, out_dir, data):
    """Perform multicore trimming with atropos.
    """
    report_file = os.path.join(out_dir, "%s-report.json" % utils.splitext_plus(os.path.basename(fastq_files[0]))[0])
    out_files = [os.path.join(out_dir, "%s-trimmed.fq.gz" % utils.splitext_plus(os.path.basename(x))[0])
                 for x in fastq_files]
    if not utils.file_exists(out_files[0]):
        with file_transaction(data, *[report_file] + out_files) as tx_out:
            tx_report_file, tx_out1 = tx_out[:2]
            if len(tx_out) > 2:
                tx_out2 = tx_out[2]
            # polyX trimming, anchored to the 3' ends of reads
            if "polyx" in dd.get_adapters(data):
                adapters += ["A{200}", "C{200}", "G{200}", "T{200}"]
            adapters_args = " ".join(["-a '%s'" % a for a in adapters])
            adapters_args += " --overlap 8"  # Avoid very short internal matches (default is 3)
            adapters_args += " --no-default-adapters --no-cache-adapters"  # Prevent GitHub queries and saving pickles
            aligner_args = "--aligner adapter"
            if len(fastq_files) == 1:
                cores = dd.get_num_cores(data)
                input_args = "-se %s" % objectstore.cl_input(fastq_files[0])
                output_args = "-o >(bgzip --threads {cores} -c > {tx_out1})".format(**locals())
            else:
                assert len(fastq_files) == 2, fastq_files
                cores = max(1, dd.get_num_cores(data) // 2)
                adapters_args = adapters_args + " " + " ".join(["-A '%s'" % a for a in adapters])
                input_args = "-pe1 %s -pe2 %s" % tuple([objectstore.cl_input(x) for x in fastq_files])
                output_args = ("-o >(bgzip --threads {cores} -c > {tx_out1}) "
                               "-p >(bgzip --threads {cores} -c > {tx_out2})").format(**locals())
            quality_base = "64" if dd.get_quality_format(data).lower() == "illumina" else "33"
            sample_name = dd.get_sample_name(data)
            report_args = "--report-file %s --report-formats json --sample-id %s" % (tx_report_file,
                                                                                     dd.get_sample_name(data))
            ropts = " ".join(str(x) for x in
                             config_utils.get_resources("atropos", data["config"]).get("options", []))
            extra_opts = []
            for k, alt_ks, v, want in [("--quality-cutoff", ["-q "], "5", True),
                                       ("--minimum-length", ["-m "], str(dd.get_min_read_length(data)), True),
                                       ("--nextseq-trim", [], "25", ("polyx" in dd.get_adapters(data) or
                                                                     "polyg" in dd.get_adapters(data)))]:
                if k not in ropts and not any(alt_k in ropts for alt_k in alt_ks):
                    if want:
                        extra_opts.append("%s=%s" % (k, v))
            extra_opts = " ".join(extra_opts)
            thread_args = ("--threads %s" % cores if cores > 1 else "")
            cmd = ("atropos trim {ropts} {thread_args} --quality-base {quality_base} --format fastq "
                   "{adapters_args} {input_args} {output_args} {report_args} {extra_opts}")
            do.run(cmd.format(**locals()), "Trimming with atropos: %s" % dd.get_sample_name(data))
    return out_files, report_file

# ## fastp trimming

def _fastp_trim(fastq_files, adapters, out_dir, data):
    """Perform multicore trimming with fastp (https://github.com/OpenGene/fastp)
    """
    report_file = os.path.join(out_dir, "%s-report.json" % utils.splitext_plus(os.path.basename(fastq_files[0]))[0])
    out_files = [os.path.join(out_dir, "%s-trimmed.fq.gz" % utils.splitext_plus(os.path.basename(x))[0])
                 for x in fastq_files]
    if not utils.file_exists(out_files[0]):
        with file_transaction(data, *[report_file] + out_files) as tx_out:
            tx_report = tx_out[0]
            tx_out_files = tx_out[1:]
            cmd = ["fastp", "--thread", dd.get_num_cores(data)]
            if dd.get_quality_format(data).lower() == "illumina":
                cmd += ["--phred64"]
            for i, (inf, outf) in enumerate(zip(fastq_files, tx_out_files)):
                if i == 0:
                    cmd += ["-i", inf, "-o", outf]
                else:
                    cmd += ["-I", inf, "-O", outf]
            cmd += ["--cut_by_quality3", "--cut_mean_quality", "5",
                    "--length_required", str(dd.get_min_read_length(data)),
                    "--disable_quality_filtering"]
            if "polyx" in dd.get_adapters(data):
                cmd += ["--trim_poly_x", "--poly_x_min_len", "8"]
            if "polyx" in dd.get_adapters(data) or "polyg" in dd.get_adapters(data):
                cmd += ["--trim_poly_g", "--poly_g_min_len", "8"]
            for a in adapters:
                cmd += ["--adapter_sequence", a]
            if not adapters:
                cmd += ["--disable_adapter_trimming"]
            cmd += ["--json", report_file, "--report_title", dd.get_sample_name(data)]
            do.run(cmd, "Trimming with fastp: %s" % dd.get_sample_name(data))
    return out_files, report_file

def _get_sequences_to_trim(config, builtin):
    builtin_adapters = _get_builtin_adapters(config, builtin)
    polya = builtin_adapters.get("polya", [None])[0]
    # allow for trimming of custom sequences for advanced users
    custom_trim = config["algorithm"].get("custom_trim", [])
    builtin_adapters = {k: v for k, v in builtin_adapters.items() if
                        k != "polya"}
    trim_sequences = custom_trim
    # for unstranded RNA-seq, libraries, both polyA and polyT can appear
    # at the 3' end as well
    if polya:
        trim_sequences += [polya, str(Seq(polya).reverse_complement())]

    # also trim the reverse complement of the adapters
    for _, v in builtin_adapters.items():
        trim_sequences += [str(Seq(sequence)) for sequence in v]
        trim_sequences += [str(Seq(sequence).reverse_complement()) for
                           sequence in v]
    out = []
    for trim in trim_sequences:
        if utils.file_exists(trim):
            out.append("file:%s" % trim)
        else:
            out.append(trim)
    return out

def _get_quality_format(config):
    SUPPORTED_FORMATS = ["illumina", "standard"]
    quality_format = dd.get_quality_format(data).lower()
    if quality_format not in SUPPORTED_FORMATS:
        logger.error("quality_format is set to an unsupported format. "
                     "Supported formats are %s."
                     % (", ".join(SUPPORTED_FORMATS)))
        exit(1)
    return quality_format

def _get_builtin_adapters(config, builtin):
    chemistries = config["algorithm"].get("adapters", [])
    adapters = {chemistry: builtin[chemistry] for
                chemistry in chemistries if chemistry in builtin}
    return adapters

# Older cutadapt approach, to remove after review

def _cutadapt_trim(fastq_files, quality_format, adapters, out_files, log_file, data):
    """Trimming with cutadapt.
    """
    if all([utils.file_exists(x) for x in out_files]):
        return out_files
    cmd = _cutadapt_trim_cmd(fastq_files, quality_format, adapters, out_files, data)
    if len(fastq_files) == 1:
        of = [out_files[0], log_file]
        message = "Trimming %s in single end mode with cutadapt." % (fastq_files[0])
        with file_transaction(data, of) as of_tx:
            of1_tx, log_tx = of_tx
            do.run(cmd.format(**locals()), message)
    else:
        of = out_files + [log_file]
        with file_transaction(data, of) as tx_out_files:
            of1_tx, of2_tx, log_tx = tx_out_files
            tmp_fq1 = utils.append_stem(of1_tx, ".tmp")
            tmp_fq2 = utils.append_stem(of2_tx, ".tmp")
            singles_file = of1_tx + ".single"
            message = "Trimming %s and %s in paired end mode with cutadapt." % (fastq_files[0],
                                                                                fastq_files[1])
            do.run(cmd.format(**locals()), message)
    return out_files

def _cutadapt_trim_cmd(fastq_files, quality_format, adapters, out_files, data):
    """Trimming with cutadapt, using version installed with bcbio-nextgen.
    """
    if all([utils.file_exists(x) for x in out_files]):
        return out_files
    if quality_format == "illumina":
        quality_base = "64"
    else:
        quality_base = "33"

    # --times=2 tries twice remove adapters which will allow things like:
    # realsequenceAAAAAAadapter to remove both the poly-A and the adapter
    # this behavior might not be what we want; we could also do two or
    # more passes of cutadapt
    cutadapt = os.path.join(os.path.dirname(sys.executable), "cutadapt")
    adapter_cmd = " ".join(map(lambda x: "-a " + x, adapters))
    ropts = " ".join(str(x) for x in
                     config_utils.get_resources("cutadapt", data["config"]).get("options", []))
    base_cmd = ("{cutadapt} {ropts} --times=2 --quality-base={quality_base} "
                "--quality-cutoff=5 --format=fastq "
                "{adapter_cmd} ").format(**locals())
    if len(fastq_files) == 2:
        # support for the single-command paired trimming introduced in
        # cutadapt 1.8
        adapter_cmd = adapter_cmd.replace("-a ", "-A ")
        base_cmd += "{adapter_cmd} ".format(adapter_cmd=adapter_cmd)
        return _cutadapt_pe_cmd(fastq_files, out_files, quality_format, base_cmd, data)
    else:
        return _cutadapt_se_cmd(fastq_files, out_files, base_cmd, data)

def _cutadapt_se_cmd(fastq_files, out_files, base_cmd, data):
    """
    this has to use the -o option, not redirect to stdout in order for
    gzipping to be supported
    """
    min_length = dd.get_min_read_length(data)
    cmd = base_cmd + " --minimum-length={min_length} ".format(**locals())
    fq1 = objectstore.cl_input(fastq_files[0])
    of1 = out_files[0]
    cmd += " -o {of1_tx} " + str(fq1)
    cmd = "%s | tee > {log_tx}" % cmd
    return cmd

def _cutadapt_pe_cmd(fastq_files, out_files, quality_format, base_cmd, data):
    """
    run cutadapt in paired end mode
    """
    fq1, fq2 = [objectstore.cl_input(x) for x in fastq_files]
    of1, of2 = out_files
    base_cmd += " --minimum-length={min_length} ".format(min_length=dd.get_min_read_length(data))
    first_cmd = base_cmd + " -o {of1_tx} -p {of2_tx} " + fq1 + " " + fq2
    return first_cmd + "| tee > {log_tx};"
