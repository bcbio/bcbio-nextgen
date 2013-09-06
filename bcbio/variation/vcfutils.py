"""Utilities for manipulating variant files in standard VCF format.
"""
import contextlib
import copy
import itertools
import os
import shutil

import pysam

from bcbio import broad, utils
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import shared
from bcbio.variation import bamprep


def is_sample_pair(align_bams, items):

    """Determine if bams are from a sample pair or  not"""

    return (len(align_bams) == 2 and all(item["metadata"].get("phenotype")
                                    is not None for item in items))


def get_paired_bams(align_bams, items):

    """Split aligned bams imn tumor / normal pairs."""

    tumor_bam = None
    normal_bam = None

    for bamfile, item in itertools.izip(align_bams, items):

        metadata = item["metadata"]

        if metadata["phenotype"] == "normal":
            normal_bam = bamfile
            normal_sample_name = item["name"][1]
        elif metadata["phenotype"] == "tumor":
            tumor_bam = bamfile
            tumor_sample_name = item["name"][1]

    if tumor_bam is None or normal_bam is None:
        raise ValueError("Missing phenotype definition (tumor or normal) "
                         "in samples")

    return (tumor_bam, tumor_sample_name, normal_bam, normal_sample_name)


def write_empty_vcf(out_file):
    with open(out_file, "w") as out_handle:
        out_handle.write("##fileformat=VCFv4.1\n"
                         "## No variants; no reads aligned in region\n"
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def split_snps_indels(orig_file, ref_file, config):
    """Split a variant call file into SNPs and INDELs for processing.
    """
    broad_runner = broad.runner_from_config(config)
    base, ext = os.path.splitext(orig_file)
    snp_file = "{base}-snp{ext}".format(base=base, ext=ext)
    indel_file = "{base}-indel{ext}".format(base=base, ext=ext)
    params = ["-T", "SelectVariants",
              "-R", ref_file,
              "--variant", orig_file]
    for out_file, select_type in [(snp_file, ["SNP"]),
                                  (indel_file, ["INDEL", "MIXED", "MNP",
                                                "SYMBOLIC", "NO_VARIATION"])]:
        if not utils.file_exists(out_file):
            with file_transaction(out_file) as tx_out_file:
                cur_params = params + ["--out", tx_out_file]
                for x in select_type:
                    cur_params += ["--selectTypeToInclude", x]
                broad_runner.run_gatk(cur_params)
    return snp_file, indel_file

def _get_exclude_samples(in_file, to_exclude):
    """Identify samples in the exclusion list which are actually in the VCF.
    """
    out = []
    with open(in_file) as in_handle:
        for line in in_handle:
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                for s in parts[9:]:
                    if s in to_exclude:
                        out.append(s)
                break
    return out

def exclude_samples(in_file, out_file, to_exclude, ref_file, config):
    """Exclude specific samples from an input VCF file using GATK SelectVariants.
    """
    to_exclude = _get_exclude_samples(in_file, to_exclude)
    # can use the input sample, all exclusions already gone
    if len(to_exclude) == 0:
        out_file = in_file
    elif not utils.file_exists(out_file):
        broad_runner = broad.runner_from_config(config)
        with file_transaction(out_file) as tx_out_file:
            params = ["-T", "SelectVariants",
                      "-R", ref_file,
                      "--variant", in_file,
                      "--out", tx_out_file]
            for x in to_exclude:
                params += ["-xl_sn", x]
            broad_runner.run_gatk(params)
    return out_file

def combine_variant_files(orig_files, out_file, ref_file, config,
                          quiet_out=True, region=None):
    """Combine multiple VCF files into a single output file.

    Handles complex merging of samples and other tricky issues using GATK.
    """
    in_pipeline = False
    if isinstance(orig_files, dict):
        file_key = config["file_key"]
        in_pipeline = True
        orig_files = orig_files[file_key]
    broad_runner = broad.runner_from_config(config)
    if not utils.file_exists(out_file):
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
            cur_region = shared.subset_variant_regions(variant_regions, region, out_file)
            if cur_region:
                params += ["-L", bamprep.region_to_gatk(cur_region),
                           "--interval_set_rule", "INTERSECTION"]
            broad_runner.run_gatk(params)
    if in_pipeline:
        return [{file_key: out_file, "region": region, "sam_ref": ref_file, "config": config}]
    else:
        return out_file

def ref_file_contigs(ref_file, config):
    """Iterator of sequence contigs from a reference file.
    """
    broad_runner = broad.runner_from_config(config)
    ref_dict = broad_runner.run_fn("picard_index_ref", ref_file)
    with contextlib.closing(pysam.Samfile(ref_dict, "r")) as ref_sam:
        for sq in ref_sam.header["SQ"]:
            yield sq

def _sort_by_region(fnames, regions, ref_file, config):
    """Sort a set of regionally split files by region for ordered output.
    """
    contig_order = {}
    for i, sq in enumerate(ref_file_contigs(ref_file, config)):
        contig_order[sq["SN"]] = i
    sitems = []
    for region, fname in zip(regions, fnames):
        if isinstance(region, (list, tuple)):
            c, s, e = region
        else:
            c = region
            s, e = 0, 0
        sitems.append(((contig_order[c], s, e), fname))
    sitems.sort()
    return [x[1] for x in sitems]

def vcf_has_variants(in_file):
    if os.path.exists(in_file):
        with open(in_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    return True
    return False

def concat_variant_files(orig_files, out_file, regions, ref_file, config):
    """Concatenate multiple variant files from regions into a single output file.

    Lightweight approach to merging VCF files split by regions with same
    sample information so no complex merging needed.
    """
    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                sorted_files = _sort_by_region(orig_files, regions, ref_file, config)
                has_variants = False
                for i, orig_file in enumerate(f for f in sorted_files if vcf_has_variants(f)):
                    has_variants = True
                    with open(orig_file) as in_handle:
                        for line in in_handle:
                            if line.startswith("#"):
                                if i == 0:
                                    out_handle.write(line)
                            else:
                                out_handle.write(line)
                # if all empty, copy the (empty of calls) first file
                if not has_variants:
                    shutil.copyfile(sorted_files[0], tx_out_file)
    return out_file

# ## Parallel VCF file combining

def parallel_combine_variants(orig_files, out_file, ref_file, config, run_parallel):
    """Combine variants in parallel by chromosome, concatenating final outputs.
    """
    file_key = "vcf_files"
    items = [[{file_key: orig_files}]]
    def split_by_region(data):
        base, ext = os.path.splitext(os.path.basename(out_file))
        args = []
        for region in [x["SN"] for x in ref_file_contigs(ref_file, config)]:
            region_out = os.path.join(os.path.dirname(out_file), "%s-regions" % base,
                                      "%s-%s%s" % (base, region, ext))
            utils.safe_makedir(os.path.dirname(region_out))
            args.append((region_out, ref_file, config, True, region))
        return out_file, args
    config = copy.deepcopy(config)
    config["file_key"] = file_key
    parallel_split_combine(items, split_by_region, run_parallel,
                           "combine_variant_files", "concat_variant_files",
                           file_key, ["region", "sam_ref", "config"], split_outfile_i=0)
    return out_file
