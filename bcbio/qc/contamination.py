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
from bcbio.variation import samtools

def run(bam_file, data, out_dir):
    out_base = os.path.join(utils.safe_makedir(out_dir),
                            "%s-verifybamid" % (dd.get_sample_name(data)))
    out_file = out_base + ".selfSM"
    failed_file = out_base + ".failed"
    exts = [".out"]
    out = {}
    if not utils.file_exists(out_file) and not utils.file_exists(failed_file):
        with file_transaction(data, out_base) as tx_out_base:
            cmd = ["verifybamid2", "1000g.phase3", "100k", "b38" if dd.get_genome_build(data) == "hg38" else "b37",
                   "--Reference", dd.get_ref_file(data), "--Output", tx_out_base]
            cmd += _get_input_args(bam_file, data, out_base)
            try:
                do.run(cmd, "VerifyBamID contamination checks")
            except subprocess.CalledProcessError, msg:
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
    if utils.file_exists(out_file):
        out["base"] = out_file
        out["secondary"] = [out_base + e for e in exts if os.path.exists(out_base + e)]
    return out

def _get_input_args(bam_file, data, out_base):
    """Retrieve input args, depending on genome build.

    VerifyBamID2 only handles GRCh37 (1, 2, 3) not hg19, so need to generate
    a pileup for hg19 and fix chromosome naming.
    """
    if dd.get_genome_build(data) in ["hg19"]:
        out_file = "%s-mpileup.txt" % out_base
        if not utils.file_exists(out_file):
            with file_transaction(data, out_file) as tx_out_file:
                mpileup_cl = samtools.prep_mpileup([bam_file], dd.get_ref_file(data), data["config"], want_bcf=False,
                                                   target_regions=_get_autosomal_bed(data, tx_out_file))
                cl = ("{mpileup_cl} | sed 's/^chr//' > {tx_out_file}")
                do.run(cl.format(**locals()), "Create pileup from BAM input")
        return ["--PileupFile", out_file]
    else:
        return ["--BamFile", bam_file]

def _get_autosomal_bed(data, base_file):
    out_file = "%s-stdchroms.bed" % utils.splitext_plus(base_file)[0]
    with open(out_file, "w") as out_handle:
        for r in ref.file_contigs(dd.get_ref_file(data)):
            if chromhacks.is_autosomal(r.name):
                out_handle.write("%s\t0\t%s\n" % (r.name, r.size))
    return out_file

