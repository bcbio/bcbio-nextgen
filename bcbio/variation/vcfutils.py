"""Utilities for manipulating variant files in standard VCF format.
"""

from collections import namedtuple, defaultdict
import copy
import os
import pprint
import shutil
import subprocess

import toolz as tz

import six
from six.moves import zip

from bcbio import broad, utils
from bcbio.bam import ref
from bcbio.distributed.multi import run_multicore, zeromq_aware_logging
from bcbio.distributed.split import parallel_split_combine
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils, tools
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do

# ## Tumor/normal paired cancer analyses

PairedData = namedtuple("PairedData", ["tumor_bam", "tumor_name",
                                       "normal_bam", "normal_name", "normal_panel",
                                       "tumor_config", "tumor_data", "normal_data"])

def is_paired_analysis(align_bams, items):
    """Determine if BAMs are from a tumor/normal paired analysis.
    """
    return get_paired_bams(align_bams, items) is not None


def somatic_batches(items):
    """Group items into somatic calling batches (tumor-only or tumor/normal).

    Returns batches, where a data item may be in pairs, and somatic and non_somatic
    (which are the original list of items).
    """
    non_somatic = []
    somatic = []
    data_by_batches = defaultdict(list)
    for data in items:
        if not get_paired_phenotype(data):
            non_somatic.append(data)
        else:
            somatic.append(data)
            batches = dd.get_batches(data)
            if batches:
                for batch in batches:
                    data_by_batches[batch].append(data)
    return data_by_batches.values(), somatic, non_somatic

def get_paired(items):
    return get_paired_bams([dd.get_align_bam(d) for d in items], items)

def get_paired_bams(align_bams, items):
    """Split aligned bams into tumor / normal pairs if this is a paired analysis.
    Allows cases with only tumor BAMs to handle callers that can work without
    normal BAMs or with normal VCF panels.
    """
    tumor_bam, tumor_name, normal_bam, normal_name, normal_panel, tumor_config, normal_data = (None,) * 7
    for bamfile, item in zip(align_bams, items):
        phenotype = get_paired_phenotype(item)
        if phenotype == "normal":
            normal_bam = bamfile
            normal_name = dd.get_sample_name(item)
            normal_data = item
        elif phenotype == "tumor":
            tumor_bam = bamfile
            tumor_name = dd.get_sample_name(item)
            tumor_data = item
            tumor_config = item["config"]
            normal_panel = dd.get_background_variant(item)
    if tumor_bam or tumor_name:
        return PairedData(tumor_bam, tumor_name, normal_bam,
                          normal_name, normal_panel, tumor_config,
                          tumor_data, normal_data)

def get_somatic_variantcallers(items):
    """Retrieve all variant callers for somatic calling, handling somatic/germline.
    """
    out = []
    for data in items:
        vcs = dd.get_variantcaller(data)
        if isinstance(vcs, dict) and "somatic" in vcs:
            vcs = vcs["somatic"]
        if not isinstance(vcs, (list, tuple)):
            vcs = [vcs]
        out += vcs
    return set(vcs)

def check_paired_problems(items):
    """Check for incorrectly paired tumor/normal samples in a batch.
    """
    # ensure we're in a paired batch
    if not get_paired(items):
        return
    num_tumor = len([x for x in items if dd.get_phenotype(x).lower() == "tumor"])
    if num_tumor > 1:
        raise ValueError("Unsupported configuration: found multiple tumor samples in batch %s: %s" %
                         (tz.get_in(["metadata", "batch"], items[0]),
                          [dd.get_sample_name(data) for data in items]))
    elif num_tumor == 0 and any(dd.get_phenotype(data).lower() == "normal" for data in items):
        raise ValueError("Found normal sample without tumor in batch %s: %s" %
                         (tz.get_in(["metadata", "batch"], items[0]),
                          [dd.get_sample_name(data) for data in items]))
    else:
        vcs = get_somatic_variantcallers(items)
        if "mutect" in vcs or "mutect2" in vcs or "strelka2" in vcs:
            paired = get_paired(items)
            if not (paired.normal_data or paired.normal_panel):
                raise ValueError("MuTect, MuTect2 and Strelka2 somatic calling requires normal sample or panel: %s" %
                                 [dd.get_sample_name(data) for data in items])

def get_paired_phenotype(data):
    """Retrieve the phenotype for a paired tumor/normal analysis.
    """
    allowed_names = set(["tumor", "normal"])
    p = tz.get_in(["metadata", "phenotype"], data)
    return p if p in allowed_names else None

# ## General utilities

def fix_ambiguous_cl(column=4):
    """awk command to replace non-N ambiguous REF bases with N.

    Some callers include these if present in the reference genome but GATK does
    not like them.
    """
    return r"""awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDXkmryswbvhdx]/, "N", $%s) } {print}'""" % column

def remove_dup_cl():
    """awk command line to remove duplicate alleles where the ref and alt are the same.
    """
    return r""" awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {next} {print}'"""

def get_indelcaller(d_or_c):
    """Retrieve string for indelcaller to use, or empty string if not specified.
    """
    config = d_or_c if isinstance(d_or_c, dict) and "config" in d_or_c else d_or_c
    indelcaller = config["algorithm"].get("indelcaller", "")
    if not indelcaller:
        indelcaller = ""
    if isinstance(indelcaller, (list, tuple)):
        indelcaller = indelcaller[0] if (len(indelcaller) > 0) else ""
    return indelcaller

def write_empty_vcf(out_file, config=None, samples=None):
    needs_bgzip = False
    if out_file.endswith(".vcf.gz"):
        needs_bgzip = True
        out_file = out_file.replace(".vcf.gz", ".vcf")
    with open(out_file, "w") as out_handle:
        format_samples = ("\tFORMAT\t" + "\t".join(samples)) if samples else ""
        out_handle.write("##fileformat=VCFv4.1\n"
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO%s\n" % (format_samples))
    if needs_bgzip:
        return bgzip_and_index(out_file, config or {})
    else:
        return out_file

def to_standardonly(in_file, ref_file, data):
    """Subset a VCF input file to standard chromosomes (1-22,X,Y,MT).
    """
    from bcbio.heterogeneity import chromhacks
    out_file = "%s-stdchrs.vcf.gz" % utils.splitext_plus(in_file)[0]
    if not utils.file_exists(out_file):
        stds = []
        for c in ref.file_contigs(ref_file):
            if chromhacks.is_nonalt(c.name):
                stds.append(c.name)
        if stds:
            with file_transaction(data, out_file) as tx_out_file:
                stds = ",".join(stds)
                in_file = bgzip_and_index(in_file, data["config"])
                cmd = "bcftools view -o {tx_out_file} -O z {in_file} {stds}"
                do.run(cmd.format(**locals()), "Subset to standard chromosomes")
    return bgzip_and_index(out_file, data["config"]) if utils.file_exists(out_file) else in_file

def split_snps_indels(orig_file, ref_file, config):
    """Split a variant call file into SNPs and INDELs for processing.
    """
    base, ext = utils.splitext_plus(orig_file)
    snp_file = "{base}-snp{ext}".format(base=base, ext=ext)
    indel_file = "{base}-indel{ext}".format(base=base, ext=ext)
    for out_file, select_arg in [(snp_file, "--types snps"),
                                 (indel_file, "--exclude-types snps")]:
        if not utils.file_exists(out_file):
            with file_transaction(config, out_file) as tx_out_file:
                bcftools = config_utils.get_program("bcftools", config)
                output_type = "z" if out_file.endswith(".gz") else "v"
                cmd = "{bcftools} view -O {output_type} {orig_file} {select_arg} > {tx_out_file}"
                do.run(cmd.format(**locals()), "Subset to SNPs and indels")
        if out_file.endswith(".gz"):
            bgzip_and_index(out_file, config)
    return snp_file, indel_file

def get_normal_sample(in_file):
    """Retrieve normal sample if normal/turmor
    """
    with utils.open_gzipsafe(in_file) as in_handle:
        for line in in_handle:
            if line.startswith("##PEDIGREE"):
                parts = line.strip().split("Original=")[1][:-1]
                return parts

def get_samples(in_file):
    """Retrieve samples present in a VCF file
    """
    with utils.open_gzipsafe(in_file) as in_handle:
        for line in in_handle:
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                return parts[9:]
    raise ValueError("Did not find sample header in VCF file %s" % in_file)

def _get_exclude_samples(in_file, to_exclude):
    """Identify samples in the exclusion list which are actually in the VCF.
    """
    include, exclude = [], []
    to_exclude = set(to_exclude)
    for s in get_samples(in_file):
        if s in to_exclude:
            exclude.append(s)
        else:
            include.append(s)
    return include, exclude

def exclude_samples(in_file, out_file, to_exclude, ref_file, config, filters=None):
    """Exclude specific samples from an input VCF file.
    """
    include, exclude = _get_exclude_samples(in_file, to_exclude)
    # can use the input sample, all exclusions already gone
    if len(exclude) == 0:
        out_file = in_file
    elif not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            bcftools = config_utils.get_program("bcftools", config)
            output_type = "z" if out_file.endswith(".gz") else "v"
            include_str = ",".join(include)
            filter_str = "-f %s" % filters if filters is not None else ""  # filters could be e.g. 'PASS,.'
            cmd = "{bcftools} view -O {output_type} -s {include_str} {filter_str} {in_file} > {tx_out_file}"
            do.run(cmd.format(**locals()), "Exclude samples: {}".format(to_exclude))
    return out_file

def select_sample(in_file, sample, out_file, config, filters=None):
    """Select a single sample from the supplied multisample VCF file.
    """
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            if len(get_samples(in_file)) == 1:
                shutil.copy(in_file, tx_out_file)
            else:
                if in_file.endswith(".gz"):
                    bgzip_and_index(in_file, config)
                bcftools = config_utils.get_program("bcftools", config)
                output_type = "z" if out_file.endswith(".gz") else "v"
                filter_str = "-f %s" % filters if filters is not None else ""  # filters could be e.g. 'PASS,.'
                cmd = "{bcftools} view -O {output_type} {filter_str} {in_file} -s {sample} > {tx_out_file}"
                do.run(cmd.format(**locals()), "Select sample: %s" % sample)
    if out_file.endswith(".gz"):
        bgzip_and_index(out_file, config)
    return out_file

def vcf_has_variants(in_file):
    if os.path.exists(in_file):
        with utils.open_gzipsafe(in_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    return True
    return False

def vcf_has_nonfiltered_variants(in_file):
    if os.path.exists(in_file):
        with utils.open_gzipsafe(in_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    parts = line.split("\t")
                    if parts[6] in set(["PASS", "."]):
                        return True
    return False

# ## Merging of variant files

def merge_variant_files(orig_files, out_file, ref_file, config, region=None):
    """Combine multiple VCF files with different samples into a single output file.

    Uses bcftools merge on bgzipped input files, handling both tricky merge and
    concatenation of files. Does not correctly handle files with the same
    sample (use combine_variant_files instead).
    """
    in_pipeline = False
    if isinstance(orig_files, dict):
        file_key = config["file_key"]
        in_pipeline = True
        orig_files = orig_files[file_key]
    out_file = _do_merge(orig_files, out_file, config, region)
    if in_pipeline:
        return [{file_key: out_file, "region": region, "sam_ref": ref_file, "config": config}]
    else:
        return out_file

def _do_merge(orig_files, out_file, config, region):
    """Do the actual work of merging with bcftools merge.
    """
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            _check_samples_nodups(orig_files)
            prep_files = run_multicore(p_bgzip_and_index, [[x, config] for x in orig_files], config)
            input_vcf_file = "%s-files.txt" % utils.splitext_plus(out_file)[0]
            with open(input_vcf_file, "w") as out_handle:
                for fname in prep_files:
                    out_handle.write(fname + "\n")
            bcftools = config_utils.get_program("bcftools", config)
            output_type = "z" if out_file.endswith(".gz") else "v"
            region_str = "-r {}".format(region) if region else ""
            cmd = "{bcftools} merge -O {output_type} {region_str} `cat {input_vcf_file}` > {tx_out_file}"
            do.run(cmd.format(**locals()), "Merge variants")
    if out_file.endswith(".gz"):
        bgzip_and_index(out_file, config)
    return out_file

def _check_samples_nodups(fnames):
    """Ensure a set of input VCFs do not have duplicate samples.
    """
    counts = defaultdict(int)
    for f in fnames:
        for s in get_samples(f):
            counts[s] += 1
    duplicates = [s for s, c in counts.items() if c > 1]
    if duplicates:
        raise ValueError("Duplicate samples found in inputs %s: %s" % (duplicates, fnames))

def _sort_by_region(fnames, regions, ref_file, config):
    """Sort a set of regionally split files by region for ordered output.
    """
    contig_order = {}
    for i, sq in enumerate(ref.file_contigs(ref_file, config)):
        contig_order[sq.name] = i
    sitems = []
    assert len(regions) == len(fnames), (regions, fnames)
    added_fnames = set([])
    for region, fname in zip(regions, fnames):
        if fname not in added_fnames:
            if isinstance(region, (list, tuple)):
                c, s, e = region
            elif isinstance(region, six.string_types) and region.find(":") >= 0:
                c, coords = region.split(":")
                s, e = [int(x) for x in coords.split("-")]
            else:
                c = region
                s, e = 0, 0
            sitems.append(((contig_order[c], s, e), c, fname))
            added_fnames.add(fname)
    sitems.sort()
    return [(x[1], x[2]) for x in sitems]

def concat_variant_files(orig_files, out_file, regions, ref_file, config):
    """Concatenate multiple variant files from regions into a single output file.

    Uses GATK4's GatherVcfs, falling back to bcftools concat --naive if it fails.
    These both only combine samples and avoid parsing, allowing scaling to large
    file sizes.
    """
    if not utils.file_exists(out_file):
        input_file_list = _get_file_list(orig_files, out_file, regions, ref_file, config)
        try:
            out_file = _run_concat_variant_files_gatk4(input_file_list, out_file, config)
        except subprocess.CalledProcessError as msg:
            if ("We require all VCFs to have complete VCF headers" in str(msg)
                  or "Features added out of order" in str(msg)
                  or "The reference allele cannot be missing" in str(msg)):
                out_file = _run_concat_variant_files_bcftools(input_file_list, out_file, config, naive=True)
            else:
                print("## Original contigs")
                pprint.pprint(zip(regions, orig_files))
                print("## Ordered file list")
                with open(input_file_list) as in_handle:
                    print(in_handle.read())
                raise
    if out_file.endswith(".gz"):
        bgzip_and_index(out_file, config)
    return out_file

def _run_concat_variant_files_gatk4(input_file_list, out_file, config):
    """Use GATK4 GatherVcfs for concatenation of scattered VCFs.
    """
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            params = ["-T", "GatherVcfs", "-I", input_file_list, "-O", tx_out_file]
            # Use GATK4 for merging, tools_off: [gatk4] applies to variant calling
            config = utils.deepish_copy(config)
            if "gatk4" in dd.get_tools_off({"config": config}):
                config["algorithm"]["tools_off"].remove("gatk4")
            # Allow specification of verbosity in the unique style this tool uses
            resources = config_utils.get_resources("gatk", config)
            opts = [str(x) for x in resources.get("options", [])]
            if "--verbosity" in opts:
                params += ["--VERBOSITY:%s" % opts[opts.index("--verbosity") + 1]]
            broad_runner = broad.runner_from_config(config)
            broad_runner.run_gatk(params)
    return out_file

def _get_file_list(orig_files, out_file, regions, ref_file, config):
    """Create file with region sorted list of non-empty VCFs for concatenating.
    """
    sorted_files = _sort_by_region(orig_files, regions, ref_file, config)
    exist_files = [(c, x) for c, x in sorted_files if os.path.exists(x) and vcf_has_variants(x)]
    if len(exist_files) == 0:  # no non-empty inputs, merge the empty ones
        exist_files = [x for c, x in sorted_files if os.path.exists(x)]
    elif len(exist_files) > 1:
        exist_files = _fix_gatk_header(exist_files, out_file, config)
    else:
        exist_files = [x for c, x in exist_files]
    ready_files = run_multicore(p_bgzip_and_index, [[x, config] for x in exist_files], config)
    input_file_list = "%s-files.list" % utils.splitext_plus(out_file)[0]
    with open(input_file_list, "w") as out_handle:
        for fname in ready_files:
            out_handle.write(fname + "\n")
    return input_file_list

def _fix_gatk_header(exist_files, out_file, config):
    """Ensure consistent headers for VCF concatenation.

    Fixes problems for genomes that start with chrM by reheadering the first file.
    These files do haploid variant calling which lack the PID phasing key/value
    pair in FORMAT, so initial chrM samples cause errors during concatenation
    due to the lack of header merging. This fixes this by updating the first header.
    """
    from bcbio.variation import ploidy
    c, base_file = exist_files[0]
    replace_file = base_file
    items = [{"config": config}]
    if ploidy.get_ploidy(items, region=(c, 1, 2)) == 1:
        for c, x in exist_files[1:]:
            if ploidy.get_ploidy(items, (c, 1, 2)) > 1:
                replace_file = x
                break
    base_fix_file = os.path.join(os.path.dirname(out_file),
                                 "%s-fixheader%s" % utils.splitext_plus(os.path.basename(base_file)))
    with file_transaction(config, base_fix_file) as tx_out_file:
        header_file = "%s-header.vcf" % utils.splitext_plus(tx_out_file)[0]
        do.run("zgrep ^# %s > %s"
                % (replace_file, header_file), "Prepare header file for merging")
        resources = config_utils.get_resources("picard", config)
        ropts = []
        if "options" in resources:
            ropts += [str(x) for x in resources.get("options", [])]
        do.run("%s && picard FixVcfHeader HEADER=%s INPUT=%s OUTPUT=%s %s" %
               (utils.get_java_clprep(), header_file, base_file, base_fix_file, " ".join(ropts)),
               "Reheader initial VCF file in merge")
    bgzip_and_index(base_fix_file, config)
    return [base_fix_file] + [x for (c, x) in exist_files[1:]]

def concat_variant_files_bcftools(orig_files, out_file, config):
    if not utils.file_exists(out_file):
        exist_files = [x for x in orig_files if os.path.exists(x)]
        ready_files = run_multicore(p_bgzip_and_index, [[x, config] for x in exist_files], config)
        input_file_list = "%s-files.list" % utils.splitext_plus(out_file)[0]
        with open(input_file_list, "w") as out_handle:
            for fname in ready_files:
                out_handle.write(fname + "\n")
        return _run_concat_variant_files_bcftools(input_file_list, out_file, config)
    else:
        return bgzip_and_index(out_file, config)

def _run_concat_variant_files_bcftools(in_list, out_file, config, naive=False):
    """Concatenate variant files using bcftools concat, potentially using the fast naive option.
    """
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            bcftools = config_utils.get_program("bcftools", config)
            output_type = "z" if out_file.endswith(".gz") else "v"
            if naive:
                args = "--naive"
            else:
                args = "--allow-overlaps"
            cmd = "{bcftools} concat {args} -O {output_type} --file-list {in_list} -o {tx_out_file}"
            do.run(cmd.format(**locals()), "bcftools concat variants")
    if out_file.endswith(".gz"):
        bgzip_and_index(out_file, config)
    return out_file

def combine_variant_files(orig_files, out_file, ref_file, config,
                          quiet_out=True, region=None):
    """Combine VCF files from the same sample into a single output file.

    Handles cases where we split files into SNPs/Indels for processing then
    need to merge back into a final file.
    """
    in_pipeline = False
    if isinstance(orig_files, dict):
        file_key = config["file_key"]
        in_pipeline = True
        orig_files = orig_files[file_key]
    if not utils.file_exists(out_file):
        with file_transaction(config, out_file) as tx_out_file:
            exist_files = [x for x in orig_files if os.path.exists(x)]
            ready_files = run_multicore(p_bgzip_and_index, [[x, config] for x in exist_files], config)
            dict_file = "%s.dict" % utils.splitext_plus(ref_file)[0]
            cores = dd.get_num_cores({"config": config})
            memscale = {"magnitude": 0.9 * cores, "direction": "increase"} if cores > 1 else None
            cmd = ["picard"] + broad.get_picard_opts(config, memscale) + \
                  ["MergeVcfs", "D=%s" % dict_file, "O=%s" % tx_out_file] + \
                  ["I=%s" % f for f in ready_files]
            cmd = "%s && %s" % (utils.get_java_clprep(), " ".join(cmd))
            do.run(cmd, "Combine variant files")
    if out_file.endswith(".gz"):
        bgzip_and_index(out_file, config)
    if in_pipeline:
        return [{file_key: out_file, "region": region, "sam_ref": ref_file, "config": config}]
    else:
        return out_file

def sort_by_ref(vcf_file, data):
    """Sort a VCF file by genome reference and position, adding contig information.
    """
    out_file = "%s-prep.vcf.gz" % utils.splitext_plus(vcf_file)[0]
    if not utils.file_uptodate(out_file, vcf_file):
        with file_transaction(data, out_file) as tx_out_file:
            header_file = "%s-header.txt" % utils.splitext_plus(tx_out_file)[0]
            with open(header_file, "w") as out_handle:
                for region in ref.file_contigs(dd.get_ref_file(data), data["config"]):
                    out_handle.write("##contig=<ID=%s,length=%s>\n" % (region.name, region.size))
            cat_cmd = "zcat" if vcf_file.endswith("vcf.gz") else "cat"
            cmd = ("{cat_cmd} {vcf_file} | grep -v ^##contig | bcftools annotate -h {header_file} | "
                   "vt sort -m full -o {tx_out_file} -")
            with utils.chdir(os.path.dirname(tx_out_file)):
                do.run(cmd.format(**locals()), "Sort VCF by reference")
    return bgzip_and_index(out_file, data["config"])

def add_contig_to_header_cl(ref_file, out_file):
    """Add update ##contig lines to VCF header, required for bcftools/GATK compatibility.
    """
    header_file = "%s-contig_header.txt" % utils.splitext_plus(out_file)[0]
    with open(header_file, "w") as out_handle:
        for region in ref.file_contigs(ref_file, {}):
            out_handle.write("##contig=<ID=%s,length=%s>\n" % (region.name, region.size))
    return ("grep -v ^##contig | bcftools annotate -h %s" % header_file)

def add_contig_to_header(line, ref_file):
    """Streaming target to add contigs to a VCF file header.
    """
    if line.startswith("##fileformat=VCF"):
        out = [line]
        for region in ref.file_contigs(ref_file):
            out.append("##contig=<ID=%s,length=%s>" % (region.name, region.size))
        return "\n".join(out)
    else:
        return line

# ## Parallel VCF file combining

def parallel_combine_variants(orig_files, out_file, ref_file, config, run_parallel):
    """Combine variants in parallel by chromosome, concatenating final outputs.
    """
    file_key = "vcf_files"
    def split_by_region(data):
        base, ext = utils.splitext_plus(os.path.basename(out_file))
        args = []
        for region in [x.name for x in ref.file_contigs(ref_file, config)]:
            region_out = os.path.join(os.path.dirname(out_file), "%s-regions" % base,
                                      "%s-%s%s" % (base, region, ext))
            utils.safe_makedir(os.path.dirname(region_out))
            args.append((region_out, ref_file, config, region))
        return out_file, args
    config = copy.deepcopy(config)
    config["file_key"] = file_key
    prep_files = run_multicore(p_bgzip_and_index, [[x, config] for x in orig_files], config)
    items = [[{file_key: prep_files}]]
    parallel_split_combine(items, split_by_region, run_parallel,
                           "merge_variant_files", "concat_variant_files",
                           file_key, ["region", "sam_ref", "config"], split_outfile_i=0)
    return out_file

# ## VCF preparation

def move_vcf(orig_file, new_file):
    """Move a VCF file with associated index.
    """
    for ext in ["", ".idx", ".tbi"]:
        to_move = orig_file + ext
        if os.path.exists(to_move):
            shutil.move(to_move, new_file + ext)

def bgzip_and_index(in_file, config=None, remove_orig=True, prep_cmd="", tabix_args=None, out_dir=None):
    """bgzip and tabix index an input file, handling VCF and BED.
    """
    if config is None:
        config = {}
    out_file = in_file if in_file.endswith(".gz") else in_file + ".gz"
    if out_dir:
        remove_orig = False
        out_file = os.path.join(out_dir, os.path.basename(out_file))
    if (not utils.file_exists(out_file) or not os.path.lexists(out_file)
          or (utils.file_exists(in_file) and not utils.file_uptodate(out_file, in_file))):
        assert not in_file == out_file, "Input file is bgzipped but not found: %s" % in_file
        assert os.path.exists(in_file), "Input file %s not found" % in_file
        if not utils.file_uptodate(out_file, in_file):
            with file_transaction(config, out_file) as tx_out_file:
                bgzip = tools.get_bgzip_cmd(config)
                cat_cmd = "zcat" if in_file.endswith(".gz") else "cat"
                if prep_cmd:
                    prep_cmd = "| %s " % prep_cmd
                cmd = "{cat_cmd} {in_file} {prep_cmd} | {bgzip} -c > {tx_out_file}"
                try:
                    do.run(cmd.format(**locals()), "bgzip %s" % os.path.basename(in_file))
                except subprocess.CalledProcessError:
                    # Race conditions: ignore errors where file has been deleted by another
                    if os.path.exists(in_file) and not os.path.exists(out_file):
                        raise
            if remove_orig:
                try:
                    os.remove(in_file)
                except OSError:  # Handle cases where run in parallel and file has been deleted
                    pass
    tabix_index(out_file, config, tabix_args=tabix_args)
    return out_file

@utils.map_wrap
@zeromq_aware_logging
def p_bgzip_and_index(in_file, config):
    """Parallel-aware bgzip and indexing
    """
    return [bgzip_and_index(in_file, config)]

def _guess_preset(f):
    if f.lower().endswith(".vcf.gz"):
        return "vcf"
    elif f.lower().endswith(".bed.gz"):
        return "bed"
    elif f.lower().endswith(".gff.gz"):
        return "gff"
    else:
        raise ValueError("Unexpected tabix input: %s" % f)

def tabix_index(in_file, config, preset=None, tabix_args=None):
    """Index a file using tabix.
    """
    in_file = os.path.abspath(in_file)
    out_file = in_file + ".tbi"
    if not utils.file_exists(out_file) or not utils.file_uptodate(out_file, in_file):
        # Remove old index files to prevent linking into tx directory
        utils.remove_safe(out_file)
        with file_transaction(config, out_file) as tx_out_file:
            tabix = tools.get_tabix_cmd(config)
            tx_in_file = os.path.splitext(tx_out_file)[0]
            utils.symlink_plus(in_file, tx_in_file)
            if tabix_args:
                cmd = "{tabix} -f {tabix_args} {tx_in_file}"
            else:
                preset = _guess_preset(in_file) if preset is None else preset
                cmd = "{tabix} -f -p {preset} {tx_in_file}"
            do.run(cmd.format(**locals()), "tabix index %s" % os.path.basename(in_file))
    return out_file

def is_gvcf_file(in_file):
    """Check if an input file is raw gVCF
    """
    to_check = 100
    n = 0
    with utils.open_gzipsafe(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith("##"):
                if n > to_check:
                    break
                n += 1
                parts = line.split("\t")
                # GATK
                if parts[4] == "<NON_REF>":
                    return True
                # strelka2
                if parts[4] == "." and parts[7].startswith("BLOCKAVG"):
                    return True
                # freebayes
                if parts[4] == "<*>":
                    return True
                # platypue
                if parts[4] == "N" and parts[6] == "REFCALL":
                    return True

def cyvcf_add_filter(rec, name):
    """Add a FILTER value to a cyvcf2 record
    """
    if rec.FILTER:
        filters = rec.FILTER.split(";")
    else:
        filters = []
    if name not in filters:
        filters.append(name)
        rec.FILTER = filters
    return rec

def cyvcf_remove_filter(rec, name):
    """Remove filter with the given name from a cyvcf2 record
    """
    if rec.FILTER:
        filters = rec.FILTER.split(";")
    else:
        filters = []
    new_filters = [x for x in filters if not str(x) == name]
    if len(new_filters) == 0:
        new_filters = ["PASS"]
    rec.FILTER = new_filters
    return rec
