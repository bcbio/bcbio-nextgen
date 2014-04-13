"""Structural variant calling with Delly

https://github.com/tobiasrausch/delly
"""
import os
import subprocess

from bcbio import utils
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import vcfutils

def get_sv_exclude_file(items):
    """Retrieve SV file of regions to exclude.
    """
    sv_bed = utils.get_in(items[0], ("genome_resources", "variation", "sv_repeat"))
    if sv_bed and os.path.exists(sv_bed):
        return sv_bed

@utils.map_wrap
@zeromq_aware_logging
def _run_delly(bam_files, sv_type, ref_file, work_dir, items):
    """Run delly, calling structural variations for the specified type.
    """
    out_file = os.path.join(work_dir, "%s-svs%s.vcf"
                            % (os.path.splitext(os.path.basename(bam_files[0]))[0], sv_type))
    cores = min(utils.get_in(items[0], ("config", "algorithm", "num_cores"), 1),
                len(bam_files))
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(out_file) as tx_out_file:
            sv_exclude_bed = get_sv_exclude_file(items)
            exclude = ["-x", sv_exclude_bed] if sv_exclude_bed else []
            cmd = ["delly", "-t", sv_type, "-g", ref_file, "-o", tx_out_file] + exclude + bam_files
            multi_cmd = "export OMP_NUM_THREADS=%s && " % cores
            try:
                do.run(multi_cmd + " ".join(cmd), "delly structural variant")
            except subprocess.CalledProcessError, msg:
                # delly returns an error exit code if there are no variants
                if "No structural variants found" in str(msg):
                    vcfutils.write_empty_vcf(out_file)
                else:
                    raise
    return [out_file]

def run(items):
    """Perform detection of structural variations with delly.
    """
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural",
                                               items[0]["name"][-1], "delly"))
    work_bams = [data["align_bam"] for data in items]
    ref_file = utils.get_in(items[0], ("reference", "fasta", "base"))
    bytype_vcfs = run_multicore(_run_delly, [(work_bams, sv_type, ref_file, work_dir, items)
                                             for sv_type in ["DEL", "DUP", "INV", "TRA"]],
                                items[0]["config"])
    out_file = "%s.vcf.gz" % os.path.commonprefix(bytype_vcfs)
    delly_vcf = vcfutils.combine_variant_files(bytype_vcfs, out_file, ref_file, items[0]["config"])
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = {}
        data["sv"]["delly"] = delly_vcf
        out.append(data)
    return out
