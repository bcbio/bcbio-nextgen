"""Tumor only somatic calling with Pisces.

https://github.com/Illumina/Pisces
"""
import os
import shutil

import pysam

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, ploidy, vcfutils

def run(align_bams, items, ref_file, assoc_files, region=None, out_file=None):
    """Run tumor only pisces calling

    Handles bgzipping output file and fixing VCF sample naming to match BAM sample.
    """
    paired = vcfutils.get_paired_bams(align_bams, items)
    assert paired and not paired.normal_bam, ("Pisces supports tumor-only variant calling: %s" %
                                              (",".join([dd.get_sample_name(d) for d in items])))
    vrs = bedutils.population_variant_regions(items)
    target = shared.subset_variant_regions(vrs, region,
                                            out_file, items=items, do_merge=True)
    min_af = float(dd.get_min_allele_fraction(paired.tumor_data)) / 100.0
    if not utils.file_exists(out_file):
        base_out_name = utils.splitext_plus(os.path.basename(paired.tumor_bam))[0]
        raw_file = "%s.vcf" % utils.splitext_plus(out_file)[0]
        with file_transaction(paired.tumor_data, raw_file) as tx_out_file:
            ref_dir = _prep_genome(os.path.dirname(tx_out_file), paired.tumor_data)
            out_dir = os.path.dirname(tx_out_file)
            cores = dd.get_num_cores(paired.tumor_data)
            emit_min_af = min_af / 10.0
            cmd = ("pisces --bampaths {paired.tumor_bam} --genomepaths {ref_dir} --intervalpaths {target} "
                   "--maxthreads {cores} --minvf {emit_min_af} --vffilter {min_af} "
                   "--ploidy somatic --gvcf false -o {out_dir}")
            # Recommended filtering for low frequency indels
            # https://github.com/bcbio/bcbio-nextgen/commit/49d0cbb1f6dcbea629c63749e2f9813bd06dcee3#commitcomment-29765373
            cmd += " -RMxNFilter 5,9,0.35"
            # For low frequency UMI tagged variants, set higher variant thresholds
            # https://github.com/Illumina/Pisces/issues/14#issuecomment-399756862
            if min_af < (1.0 / 100.0):
                cmd += " --minbasecallquality 30"
            do.run(cmd.format(**locals()), "Pisces tumor-only somatic calling")
            shutil.move(os.path.join(out_dir, "%s.vcf" % base_out_name),
                        tx_out_file)
        vcfutils.bgzip_and_index(raw_file, paired.tumor_data["config"],
                                 prep_cmd="sed 's#%s.bam#%s#' | %s" %
                                 (base_out_name, dd.get_sample_name(paired.tumor_data),
                                  vcfutils.add_contig_to_header_cl(dd.get_ref_file(paired.tumor_data), out_file)))
    return vcfutils.bgzip_and_index(out_file, paired.tumor_data["config"])

def _prep_genome(out_dir, data):
    """Create prepped reference directory for pisces.

    Requires a custom GenomeSize.xml file present.
    """
    genome_name = utils.splitext_plus(os.path.basename(dd.get_ref_file(data)))[0]
    out_dir = utils.safe_makedir(os.path.join(out_dir, genome_name))
    ref_file = dd.get_ref_file(data)
    utils.symlink_plus(ref_file, os.path.join(out_dir, os.path.basename(ref_file)))
    with open(os.path.join(out_dir, "GenomeSize.xml"), "w") as out_handle:
        out_handle.write('<sequenceSizes genomeName="%s">' % genome_name)
        for c in pysam.AlignmentFile("%s.dict" % utils.splitext_plus(ref_file)[0]).header["SQ"]:
            cur_ploidy = ploidy.get_ploidy([data], region=[c["SN"]])
            out_handle.write('<chromosome fileName="%s" contigName="%s" totalBases="%s" knownBases="%s" '
                             'isCircular="false" ploidy="%s" md5="%s"/>' %
                             (os.path.basename(ref_file), c["SN"], c["LN"], c["LN"], cur_ploidy, c["M5"]))
        out_handle.write('</sequenceSizes>')
    return out_dir
