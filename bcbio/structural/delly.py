"""Structural variant calling with Delly

https://github.com/tobiasrausch/delly
"""
import os
from contextlib import closing
import copy
import itertools
import subprocess

try:
    import pybedtools
except ImportError:
    pybedtools = None
import pysam
import toolz as tz

from bcbio import utils
from bcbio.bam import callable
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.provenance import do
from bcbio.variation import vcfutils

def _get_sv_exclude_file(items):
    """Retrieve SV file of regions to exclude.
    """
    sv_bed = utils.get_in(items[0], ("genome_resources", "variation", "sv_repeat"))
    if sv_bed and os.path.exists(sv_bed):
        return sv_bed

def prepare_exclude_file(items, base_file, chrom=None):
    """Prepare a BED file for exclusion, incorporating variant regions and chromosome.

    Excludes locally repetitive regions (if `remove_lcr` is set) and
    centromere regions, both of which contribute to long run times and
    false positive structural variant calls.
    """
    out_file = "%s-exclude.bed" % utils.splitext_plus(base_file)[0]
    all_vrs = filter(lambda x: x is not None,
                     [tz.get_in(("config", "algorithm", "variant_regions"), data)
                      for data in items
                      if tz.get_in(["config", "algorithm", "coverage_interval"], data) != "genome"])
    # Get a bedtool for the full region if no variant regions
    if len(all_vrs) == 0:
        want_bedtool = callable.get_ref_bedtool(tz.get_in(["reference", "fasta", "base"], items[0]),
                                                items[0]["config"], chrom)
        lcr_bed = shared.get_lcr_bed(items)
        if lcr_bed:
            want_bedtool = want_bedtool.subtract(pybedtools.BedTool(lcr_bed))
    else:
        ready_region = shared.subset_variant_regions(tz.first(all_vrs), chrom, base_file, items)
        want_bedtool = pybedtools.BedTool(ready_region).saveas()
    sv_exclude_bed = _get_sv_exclude_file(items)
    if sv_exclude_bed and len(want_bedtool) > 0:
        want_bedtool = want_bedtool.subtract(sv_exclude_bed).saveas()
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(out_file) as tx_out_file:
            full_bedtool = callable.get_ref_bedtool(tz.get_in(["reference", "fasta", "base"], items[0]),
                                                    items[0]["config"])
            if len(want_bedtool) > 0:
                full_bedtool.subtract(want_bedtool).saveas(tx_out_file)
            else:
                full_bedtool.saveas(tx_out_file)
    return out_file

@utils.map_wrap
@zeromq_aware_logging
def _run_delly(bam_files, chrom, sv_type, ref_file, work_dir, items):
    """Run delly, calling structural variations for the specified type.
    """
    out_file = os.path.join(work_dir, "%s-svs%s-%s.vcf"
                            % (os.path.splitext(os.path.basename(bam_files[0]))[0], sv_type, chrom))
    cores = min(utils.get_in(items[0], ("config", "algorithm", "num_cores"), 1),
                len(bam_files))
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(out_file) as tx_out_file:
            exclude = ["-x", prepare_exclude_file(items, out_file, chrom)]
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
    return [_clean_bgzip_delly(out_file)]

def _clean_bgzip_delly(in_file):
    """Provide bgzipped input files, removing problem GL specifications from output.

    GATK does not like missing GLs like '.,.,.'. This converts them to the recognized '.'
    """
    out_file = "%s.gz" % in_file
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cmd = "sed 's/\.,\.,\././g' {in_file} | bgzip -c > {tx_out_file}"
            do.run(cmd.format(**locals()), "Clean and bgzip delly output")
    return out_file

def run(items):
    """Perform detection of structural variations with delly.
    """
    work_dir = utils.safe_makedir(os.path.join(items[0]["dirs"]["work"], "structural",
                                               items[0]["name"][-1], "delly"))
    work_bams = [data["align_bam"] for data in items]
    ref_file = utils.get_in(items[0], ("reference", "fasta", "base"))
    # Add core request for delly
    config = copy.deepcopy(items[0]["config"])
    delly_config = utils.get_in(config, ("resources", "delly"), {})
    delly_config["cores"] = len(items)
    config["resources"]["delly"] = delly_config
    parallel = {"type": "local", "cores": config["algorithm"].get("num_cores", 1),
                "progs": ["delly"]}
    sv_types = ["DEL", "DUP", "INV"]  # "TRA" has invalid VCF END specifications that GATK doesn't like
    with closing(pysam.Samfile(work_bams[0], "rb")) as pysam_work_bam:
        bytype_vcfs = run_multicore(_run_delly, [(work_bams, chrom, sv_type, ref_file, work_dir, items)
                                                 for (chrom, sv_type)
                                                 in itertools.product(pysam_work_bam.references, sv_types)],
                                    config, parallel)
    out_file = "%s.vcf.gz" % os.path.commonprefix(bytype_vcfs)
    delly_vcf = vcfutils.combine_variant_files(bytype_vcfs, out_file, ref_file, items[0]["config"])
    out = []
    for data in items:
        if "sv" not in data:
            data["sv"] = {}
        data["sv"]["delly"] = delly_vcf
        out.append(data)
    return out
