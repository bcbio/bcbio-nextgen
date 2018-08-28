"""Structural variant detection with GRIDSS

The Genomic Rearrangement IDentification Software Suite
https://github.com/PapenfussLab/gridss
"""
import os

import toolz as tz

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.pipeline import config_utils
from bcbio.provenance import do
from bcbio.structural import shared as sshared
from bcbio.variation import effects, vcfutils

def run(items, background=None):
    """Perform detection of structural variations with Manta.
    """
    paired = vcfutils.get_paired(items)
    if paired:
        inputs = [paired.tumor_data]
        background = [paired.normal_data] if paired.normal_bam else []
    else:
        assert not background
        inputs, background = sshared.find_case_control(items)
    work_dir = _sv_workdir(inputs[0])
    variant_file = _run_gridss(inputs, background, work_dir)
    out = []
    for data in items:
        sample_file = variant_file
        if "sv" not in data:
            data["sv"] = []
        effects_vcf, _ = effects.add_to_vcf(sample_file, data, "snpeff")
        data["sv"].append({"variantcaller": "gridss",
                           "vrn_file": effects_vcf or sample_file})
        out.append(data)
    return out

def _run_gridss(inputs, background, work_dir):
    out_file = os.path.join(work_dir, "%s-gridss.sv.vcf" % (dd.get_batch(inputs[0]) or
                                                            dd.get_sample_name(inputs[0])))
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(inputs[0], out_file) as tx_out_file:
            htsjdk_opts = ["-Dsamjdk.create_index=true", "-Dsamjdk.use_async_io_read_samtools=true",
                           "-Dsamjdk.use_async_io_write_samtools=true", "-Dsamjdk.use_async_io_write_tribble=true"]
            cores = dd.get_cores(inputs[0])
            resources = config_utils.get_resources("gridss", inputs[0]["config"])
            jvm_opts = resources.get("jvm_opts", ["-Xms750m", "-Xmx4g"])
            jvm_opts = config_utils.adjust_opts(jvm_opts, {"algorithm": {"memory_adjust":
                                                                         {"direction": "increase",
                                                                          "magnitude": cores}}})
            jvm_opts = _finalize_memory(jvm_opts)
            tx_ref_file = _setup_reference_files(inputs[0], os.path.dirname(tx_out_file))
            blacklist_bed = sshared.prepare_exclude_file(inputs + background, out_file)
            cmd = ["gridss"] + jvm_opts + htsjdk_opts + ["gridss.CallVariants"] + \
                  ["THREADS=%s" % cores,
                   "TMP_DIR=%s" % os.path.dirname(tx_out_file), "WORKING_DIR=%s" % os.path.dirname(tx_out_file),
                   "OUTPUT=%s" % tx_out_file,
                   "ASSEMBLY=%s" % tx_out_file.replace(".sv.vcf", ".gridss.assembly.bam"),
                   "REFERENCE_SEQUENCE=%s" % tx_ref_file, "BLACKLIST=%s" % blacklist_bed]
            for data in inputs + background:
                cmd += ["INPUT=%s" % dd.get_align_bam(data), "INPUT_LABEL=%s" % dd.get_sample_name(data)]
            exports = utils.local_path_export()
            cmd = exports + " ".join(cmd)
            do.run(cmd, "GRIDSS SV analysis")
    return vcfutils.bgzip_and_index(out_file, inputs[0]["config"])

def _finalize_memory(jvm_opts):
    """GRIDSS does not recommend setting memory between 32 and 48Gb.

    https://github.com/PapenfussLab/gridss#memory-usage
    """
    avoid_min = 32
    avoid_max = 48
    out_opts = []
    for opt in jvm_opts:
        if opt.startswith("-Xmx"):
            spec = opt[4:]
            val = int(spec[:-1])
            mod = spec[-1]
            if mod.upper() == "M":
                adjust = 1024
                min_val = avoid_min * 1024
                max_val = avoid_max * 1024
            else:
                adjust = 1
                min_val, max_val = avoid_min, avoid_max
            if val >= min_val and val < max_val:
                val = min_val - adjust
            opt = "%s%s%s" % (opt[:4], val, mod)
        out_opts.append(opt)
    return out_opts

def _setup_reference_files(data, tx_out_dir):
    """Create a reference directory with fasta and bwa indices.

    GRIDSS requires all files in a single directory, so setup with symlinks.
    This needs bwa aligner indices available, which we ensure with `get_aligner_with_aliases`
    during YAML sample setup.
    """
    aligner = dd.get_aligner(data) or "bwa"
    out_dir = utils.safe_makedir(os.path.join(tx_out_dir, aligner))
    ref_fasta = dd.get_ref_file(data)
    ref_files = ["%s%s" % (utils.splitext_plus(ref_fasta)[0], ext) for ext in [".fa", ".fa.fai", ".dict"]]
    for orig_file in ref_files + tz.get_in(("reference", aligner, "indexes"), data):
        utils.symlink_plus(orig_file, os.path.join(out_dir, os.path.basename(orig_file)))
    return os.path.join(out_dir, os.path.basename(ref_fasta))

def _sv_workdir(data):
    return os.path.join(
        data["dirs"]["work"], "structural", dd.get_sample_name(data), "gridss")
