"""Examine callable regions following genome mapping of short reads.

Identifies callable analysis regions surrounded by larger regions lacking
aligned bases. This allows parallelization of smaller chromosome chunks
through post-processing and variant calling, with each sub-section
mapping handled separately.
"""
import contextlib
import os

import pybedtools
import pysam

from bcbio import utils, broad
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction

def calc_callable_loci(in_bam, ref_file, config):
    """Determine callable bases for input BAM using Broad's CallableLoci walker.

    http://www.broadinstitute.org/gatk/gatkdocs/
    org_broadinstitute_sting_gatk_walkers_coverage_CallableLoci.html
    """
    broad_runner = broad.runner_from_config(config)
    out_bed = "%s-callable.bed" % os.path.splitext(in_bam)[0]
    out_summary = "%s-callable-summary.txt" % os.path.splitext(in_bam)[0]
    variant_regions = config["algorithm"].get("variant_regions", None)
    if not utils.file_exists(out_bed):
        with file_transaction(out_bed) as tx_out_file:
            broad_runner.run_fn("picard_index", in_bam)
            params = ["-T", "CallableLoci",
                      "-R", ref_file,
                      "-I", in_bam,
                      "--out", tx_out_file,
                      "--summary", out_summary]
            if variant_regions:
                params += ["-L", variant_regions]
            broad_runner.run_gatk(params)
    return out_bed

def get_ref_bedtool(ref_file, config):
    """Retrieve a pybedtool BedTool object with reference sizes from input reference.
    """
    broad_runner = broad.runner_from_config(config)
    ref_dict = broad_runner.run_fn("picard_index_ref", ref_file)
    ref_lines = []
    with contextlib.closing(pysam.Samfile(ref_dict, "r")) as ref_sam:
        for sq in ref_sam.header["SQ"]:
            ref_lines.append("%s\t%s\t%s" % (sq["SN"], 0, sq["LN"]))
    return pybedtools.BedTool("\n".join(ref_lines), from_string=True)

def _get_nblock_regions(in_file, min_n_size):
    """Retrieve coordinates of regions in reference genome with no mapping.
    These are potential breakpoints for parallelizing analysis.
    """
    out_lines = []
    with open(in_file) as in_handle:
        for line in in_handle:
            contig, start, end, ctype = line.rstrip().split()
            if (ctype in ["REF_N", "NO_COVERAGE"] and
                  int(end) - int(start) > min_n_size):
                out_lines.append("%s\t%s\t%s\n" % (contig, start, end))
    return pybedtools.BedTool("\n".join(out_lines), from_string=True)

def block_regions(in_bam, ref_file, config):
    """Find blocks of regions for analysis from mapped input BAM file.

    Identifies islands of callable regions, surrounding by regions
    with no read support, that can be analyzed independently.
    """
    min_n_size = int(config["algorithm"].get("nomap_split_size", 2000))
    callable_bed = calc_callable_loci(in_bam, ref_file, config)
    ref_regions = get_ref_bedtool(ref_file, config)
    nblock_regions = _get_nblock_regions(callable_bed, min_n_size)
    return [(r.chrom, int(r.start), int(r.stop)) for r in ref_regions.subtract(nblock_regions)]

def combine_sample_regions(samples):
    """Combine islands of callable regions from multiple samples.
    Creates a global set of callable samples usable across a
    project with multi-sample calling.
    """
    min_n_size = int(samples[0]["config"]["algorithm"].get("nomap_split_size", 2000))
    final_regions = None
    for regions in (x["regions"] for x in samples):
        bed_lines = ["%s\t%s\t%s" % (c, s, e) for (c, s, e) in regions]
        bed_regions = pybedtools.BedTool("\n".join(bed_lines), from_string=True)
        if final_regions is None:
            final_regions = bed_regions
        else:
            final_regions = final_regions.merge(bed_regions, d=min_n_size)
    return [(r.chrom, int(r.start), int(r.stop)) for r in final_regions]
