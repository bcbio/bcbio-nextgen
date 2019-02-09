"""Low frequency somatic variant calling with smCounter2.

https://github.com/qiaseq/qiaseq-smcounter-v2
"""
import glob
import os
import shutil

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bedutils, vcfutils


def run(align_bams, items, ref_file, assoc_files, region=None, out_file=None):
    """Run tumor only smCounter2 calling.
    """
    paired = vcfutils.get_paired_bams(align_bams, items)
    assert paired and not paired.normal_bam, ("smCounter2 supports tumor-only variant calling: %s" %
                                              (",".join([dd.get_sample_name(d) for d in items])))
    vrs = bedutils.population_variant_regions(items)
    target = shared.subset_variant_regions(vrs, region,
                                            out_file, items=items, do_merge=True)
    out_file = out_file.replace(".vcf.gz", ".vcf")
    out_prefix = utils.splitext_plus(os.path.basename(out_file))[0]
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            cmd = ["smCounter2", "--runPath", os.path.dirname(tx_out_file),
                   "--outPrefix", out_prefix,
                   "--bedTarget", target, "--refGenome", ref_file,
                   "--bamFile", paired.tumor_bam, "--bamType", "consensus",
                   "--nCPU", dd.get_num_cores(paired.tumor_data)]
            do.run(cmd, "smcounter2 variant calling")
            for fname in glob.glob(os.path.join(os.path.dirname(tx_out_file), "*.smCounter*")):
                shutil.move(fname, os.path.join(os.path.dirname(out_file), os.path.basename(fname)))
            utils.symlink_plus(os.path.join(os.path.dirname(out_file),
                                            "%s.smCounter.cut.vcf" % out_prefix),
                               out_file)
    return vcfutils.bgzip_and_index(out_file, paired.tumor_data["config"], remove_orig=False,
                                    prep_cmd="sed 's#FORMAT\t%s#FORMAT\t%s#' | %s" %
                                    (out_prefix, dd.get_sample_name(paired.tumor_data),
                                     vcfutils.add_contig_to_header_cl(dd.get_ref_file(paired.tumor_data), out_file)))
