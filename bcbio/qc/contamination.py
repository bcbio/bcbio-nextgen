"""Estimate contamination using VerifyBamID2.
"""
import os
import shutil
import subprocess
import sys

from bcbio import utils
from bcbio.log import logger
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.qc import samtools as qc_samtools
from bcbio.variation import samtools

def run(bam_file, data, out_dir):
    out_base = os.path.join(utils.safe_makedir(out_dir),
                            "%s-verifybamid" % (dd.get_sample_name(data)))
    out_file = out_base + ".selfSM"
    failed_file = out_base + ".failed"
    exts = [".out"]
    out = {}
    if not utils.file_exists(out_file) and not utils.file_exists(failed_file):
        _generate_estimates(bam_file, out_base, failed_file, exts, data)
    if utils.file_exists(out_file):
        out["base"] = out_file
        out["secondary"] = [out_base + e for e in exts if os.path.exists(out_base + e)]
    return out

def _generate_estimates(bam_file, out_base, failed_file, exts, data):
    background = {"dataset": "1000g.phase3",
                  "nvars": "100k",
                  "build":"b38" if dd.get_genome_build(data) == "hg38" else "b37"}
    with file_transaction(data, out_base) as tx_out_base:
        cmd = ["verifybamid2", background["dataset"], background["nvars"], background["build"],
               "--Reference", dd.get_ref_file(data), "--Output", tx_out_base]
        cmd += _get_input_args(bam_file, data, out_base, background)
        try:
            do.run(cmd, "VerifyBamID contamination checks")
        except subprocess.CalledProcessError as msg:
            def allowed_errors(l):
                return (l.find("Insufficient Available markers") >= 0 or
                        l.find("No reads found in any of the regions") >= 0)
            if any([allowed_errors(l) for l in str(msg).split("\n")]):
                logger.info("Skipping VerifyBamID, not enough overlapping markers found: %s" %
                            dd.get_sample_name(data))
                with open(failed_file, "w") as out_handle:
                    out_handle.write(str(msg))
            else:
                logger.warning(str(msg))
                raise
        else:
            # Fix any sample name problems, for pileups
            shutil.move(tx_out_base + ".selfSM", tx_out_base + ".selfSM.orig")
            with open(tx_out_base + ".selfSM.orig") as in_handle:
                with open(tx_out_base + ".selfSM", "w") as out_handle:
                    sample_name = None
                    for line in in_handle:
                        if line.startswith("DefaultSampleName"):
                            line = line.replace("DefaultSampleName", dd.get_sample_name(data))
                        # work around bug in finding SM from BAM RG at end of line
                        if len(line.strip().split("\t")) == 1:
                            sample_name = line.strip()
                            line = None
                        elif sample_name:
                            parts = line.split("\t")
                            parts[0] = sample_name
                            line = "\t".join(parts)
                            sample_name = None
                        if line:
                            out_handle.write(line)
            for e in exts + [".selfSM"]:
                if os.path.exists(tx_out_base + e):
                    shutil.copy(tx_out_base + e, out_base + e)

def _get_input_args(bam_file, data, out_base, background):
    """Retrieve input args, depending on genome build.

    VerifyBamID2 only handles GRCh37 (1, 2, 3) not hg19, so need to generate
    a pileup for hg19 and fix chromosome naming.
    """
    if dd.get_genome_build(data) in ["hg19"]:
        return ["--PileupFile", _create_pileup(bam_file, data, out_base, background)]
    else:
        return ["--BamFile", bam_file]

def _create_pileup(bam_file, data, out_base, background):
    """Create pileup calls in the regions of interest for hg19 -> GRCh37 chromosome mapping.
    """
    out_file = "%s-mpileup.txt" % out_base
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            background_bed = os.path.normpath(os.path.join(
                os.path.dirname(os.path.realpath(utils.which("verifybamid2"))),
                "resource", "%s.%s.%s.vcf.gz.dat.bed" % (background["dataset"],
                                                         background["nvars"], background["build"])))
            local_bed = os.path.join(os.path.dirname(out_base),
                                     "%s.%s-hg19.bed" % (background["dataset"], background["nvars"]))
            if not utils.file_exists(local_bed):
                with file_transaction(data, local_bed) as tx_local_bed:
                    with open(background_bed) as in_handle:
                        with open(tx_local_bed, "w") as out_handle:
                            for line in in_handle:
                                out_handle.write("chr%s" % line)
            mpileup_cl = samtools.prep_mpileup([bam_file], dd.get_ref_file(data), data["config"], want_bcf=False,
                                                target_regions=local_bed)
            cl = ("{mpileup_cl} | sed 's/^chr//' > {tx_out_file}")
            do.run(cl.format(**locals()), "Create pileup from BAM input")
    return out_file
