"""Utilities for manipulating variant files in standard VCF format.
"""
import os

from bcbio import broad
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.variation import bamprep

def split_snps_indels(broad_runner, orig_file, ref_file):
    """Split a variant call file into SNPs and INDELs for processing.
    """
    base, ext = os.path.splitext(orig_file)
    snp_file = "{base}-snp{ext}".format(base=base, ext=ext)
    indel_file = "{base}-indel{ext}".format(base=base, ext=ext)
    params = ["-T", "SelectVariants",
              "-R", ref_file,
              "--variant", orig_file]
    for out_file, select_type in [(snp_file, ["SNP"]),
                                  (indel_file, ["INDEL", "MIXED", "MNP",
                                                "SYMBOLIC", "NO_VARIATION"])]:
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                cur_params = params + ["--out", tx_out_file]
                for x in select_type:
                    cur_params += ["--selectTypeToInclude", x]
                broad_runner.run_gatk(cur_params)
    return snp_file, indel_file

def combine_variant_files(orig_files, out_file, ref_file, config,
                          quiet_out=True):
    """Combine multiple VCF files into a single output file.

    Handles complex merging of samples and other tricky issues using GATK.
    """
    broad_runner = broad.runner_from_config(config)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "CombineVariants",
                      "-R", ref_file,
                      "--out", tx_out_file]
            priority_order = []
            for i, orig_file in enumerate(orig_files):
                name = "v%s" % i
                params.extend(["--variant:{name}".format(name=name), orig_file])
                priority_order.append(name)
            params.extend(["--rod_priority_list", ",".join(priority_order)])
            if quiet_out:
                params.extend(["--suppressCommandLineHeader", "--setKey", "null"])
            variant_regions = config["algorithm"].get("variant_regions", None)
            if variant_regions:
                params += ["-L", bamprep.region_to_gatk(variant_regions),
                           "--interval_set_rule", "INTERSECTION"]
            broad_runner.run_gatk(params)
    return out_file

def concat_variant_files(orig_files, out_file, ref_file, config):
    """Concatenate multiple variant files from regions into a single output file.

    Lightweight approach to merging VCF files split by regions.
    """
    # XXX defaults to more general approach until implemented
    return combine_variant_files(orig_files, out_file, ref_file, config)
