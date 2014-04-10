"""Structural variant calling with Delly

https://github.com/tobiasrausch/delly
"""
import os
import subprocess

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import vcfutils

def _run_delly(bam_files, sv_type, ref_file, work_dir):
    """Run delly, calling structural variations for the specified type.
    """
    out_file = os.path.join(work_dir, "%s-svs%s.vcf"
                            % (os.path.splitext(os.path.basename(bam_files[0]))[0], sv_type))
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(out_file) as tx_out_file:
            cmd = ["delly", "-g", ref_file, "-o", tx_out_file] + bam_files
            try:
                do.run(cmd, "delly structural variant")
            except subprocess.CalledProcessError, msg:
                # delly returns an error exit code if there are no variants
                if "No structural variants found" in str(msg):
                    vcfutils.write_empty_vcf(out_file)
                else:
                    raise
    return out_file

def run(items):
    """Perform detection of structural variations with delly.
    """
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural",
                                               items[0]["name"][-1], "delly"))
    work_bams = [data["work_bam"] for data in items]
    ref_file = utils.get_in(items[0], ("reference", "fasta", "base"))
    bytype_vcfs = [_run_delly(work_bams, sv_type, ref_file, work_dir) for sv_type in ["DEL", "DUP", "INV", "TRA"]]
    out_file = "%s.vcf.gz" % os.path.commonprefix(bytype_vcfs)
    delly_vcf = vcfutils.combine_variant_files(bytype_vcfs, out_file, ref_file, items[0]["config"])
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = {}
        data["sv"]["delly"] = delly_vcf
        out.append(data)
    return out
