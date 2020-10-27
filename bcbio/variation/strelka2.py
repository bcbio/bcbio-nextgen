"""Germline and somatic calling with Strelka2: https://github.com/illumina/strelka
"""
import collections
import os
import six
import sys
import numpy as np

from bcbio import utils
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import shared
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import bamprep, bedutils, joint, ploidy, vcfutils

cyvcf2 = utils.LazyImport("cyvcf2")

def run(align_bams, items, ref_file, assoc_files, region, out_file):
    """Run strelka2 variant calling, either paired tumor/normal or germline calling.

    region can be a single region or list of multiple regions for multicore calling.
    """
    call_file = "%s-raw.vcf.gz" % utils.splitext_plus(out_file)[0]
    strelka_work_dir = "%s-work" % utils.splitext_plus(out_file)[0]
    paired = vcfutils.get_paired_bams(align_bams, items)
    if paired:
        assert paired.normal_bam, "Strelka2 requires a normal sample"
        call_file = _run_somatic(paired, ref_file, assoc_files, region, call_file, strelka_work_dir)
    else:
        call_file = _run_germline(align_bams, items, ref_file,
                                  assoc_files, region, call_file, strelka_work_dir)
    return _af_annotate_and_filter(paired, items, call_file, out_file)

def get_region_bed(region, items, out_file, want_gzip=True):
    """Retrieve BED file of regions to analyze, either single or multi-region.
    """
    variant_regions = bedutils.population_variant_regions(items, merged=True)
    target = shared.subset_variant_regions(variant_regions, region, out_file, items)
    if not target:
        raise ValueError("Need BED input for strelka2 regions: %s %s" % (region, target))
    if not isinstance(target, six.string_types) or not os.path.isfile(target):
        chrom, start, end = target
        target = "%s-regions.bed" % utils.splitext_plus(out_file)[0]
        with file_transaction(items[0], target) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write("%s\t%s\t%s\n" % (chrom, start, end))
    out_file = target
    if want_gzip:
        out_file = vcfutils.bgzip_and_index(out_file, items[0]["config"])
    if bedutils.has_regions(out_file):
        return out_file

def coverage_interval_from_bed(bed_file, per_chrom=True):
    """Calculate a coverage interval for the current region BED.

    This helps correctly work with cases of uneven coverage across an analysis
    genome. strelka2 and other model based callers have flags for targeted and non
    which depend on the local context.

    Checks coverage per chromosome, avoiding non-standard chromosomes, if per_chrom is set.
    Otherwise does a global check over all regions. The global check performs better for
    strelka2 but not for DeepVariant:

    https://github.com/bcbio/bcbio_validations/tree/master/deepvariant#deepvariant-v06-release-strelka2-stratification-and-initial-gatk-cnn
    """
    total_starts = {}
    total_ends = {}
    bed_bases = collections.defaultdict(int)
    with utils.open_gzipsafe(bed_file) as in_handle:
        for line in in_handle:
            parts = line.split()
            if len(parts) >= 3:
                chrom, start, end = parts[:3]
                if chromhacks.is_autosomal(chrom):
                    start = int(start)
                    end = int(end)
                    bed_bases[chrom] += (end - start)
                    total_starts[chrom] = min([start, total_starts.get(chrom, sys.maxsize)])
                    total_ends[chrom] = max([end, total_ends.get(chrom, 0)])
    # can check per chromosome -- any one chromosome with larger, or over all regions
    if per_chrom:
        freqs = [float(bed_bases[c]) / float(total_ends[c] - total_starts[c]) for c in sorted(bed_bases.keys())]
    elif len(bed_bases) > 0:
        freqs = [sum([bed_bases[c] for c in sorted(bed_bases.keys())]) /
                 sum([float(total_ends[c] - total_starts[c]) for c in sorted(bed_bases.keys())])]
    else:
        freqs = []
    # Should be importing GENOME_COV_THRESH but get circular imports
    if any([f >= 0.40 for f in freqs]):
        return "genome"
    else:
        return "targeted"

def _is_targeted_region(cur_bed, data):
    """Calculate if we should process region as a targeted or WGS.

    Currently always based on total coverage interval, as that validates best and
    is consistent between CWL (larger blocks) and non-CWL runs (smaller blocks).
    We can check core usage and provide a consistent report when moving to CWL
    exclusively.
    """
    cores = dd.get_num_cores(data)
    if cores > 0:  # Apply to all core setups now for consistency
        return dd.get_coverage_interval(data) not in ["genome"]
    else:
        return coverage_interval_from_bed(cur_bed, per_chrom=False) == "targeted"

def _get_ploidy(regions, items, base_file):
    samples = [dd.get_sample_name(d) for d in items]
    out_file = "%s-ploidy.vcf" % utils.splitext_plus(base_file)[0]
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(items[0], out_file) as tx_outfile:
            with open(tx_outfile, "w") as h:
                h.write("##fileformat=VCFv4.1\n")
                h.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
                h.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
                h.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
                for region in regions:
                    ploidies = [ploidy.get_ploidy([d], region) for d in items]
                    h.write("\t".join([region[0], str(region[1]), ".", "N", "<CNV>", ".", ".",
                                       "END=%s" % region[2], "CN"] + [str(x) for x in ploidies]) + "\n")
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _configure_germline(align_bams, items, ref_file, region, out_file, tx_work_dir):
    utils.safe_makedir(tx_work_dir)
    cmd = [utils.get_program_python("configureStrelkaGermlineWorkflow.py"),
           os.path.realpath(utils.which("configureStrelkaGermlineWorkflow.py"))]
    cur_bed = get_region_bed(region, items, out_file)
    if cur_bed:
        cmd += ["--referenceFasta=%s" % ref_file,
                "--callRegions=%s" % cur_bed,
                "--ploidy=%s" % _get_ploidy(shared.to_multiregion(region), items, out_file),
                "--runDir=%s" % tx_work_dir]
        cmd += ["--bam=%s" % b for b in align_bams]
        if _is_targeted_region(cur_bed, items[0]):
            cmd += ["--targeted"]
        do.run(cmd, "Configure Strelka2 germline calling: %s" % (", ".join([dd.get_sample_name(d) for d in items])))
        return os.path.join(tx_work_dir, "runWorkflow.py")

def _run_germline(align_bams, items, ref_file, assoc_files, region, out_file, work_dir):
    if not utils.file_exists(out_file):
        with file_transaction(items[0], work_dir) as tx_work_dir:
            workflow_file = _configure_germline(align_bams, items, ref_file, region, out_file, tx_work_dir)
            if workflow_file:
                has_variants = True
                _run_workflow(items[0], workflow_file, tx_work_dir)
            else:
                has_variants = False
                vcfutils.write_empty_vcf(out_file, items[0]["config"], [dd.get_sample_name(d) for d in items])
        if has_variants:
            raw_file = os.path.join(work_dir, "results", "variants",
                                    "genome.vcf.gz" if joint.want_gvcf(items) else "variants.vcf.gz")
            utils.copy_plus(raw_file, out_file)
            # Remove files with relative symlinks
            utils.remove_plus(os.path.join(work_dir, "results", "variants", "genome.vcf.gz"))
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _configure_somatic(paired, ref_file, region, out_file, tx_work_dir):
    utils.safe_makedir(tx_work_dir)
    cmd = [utils.get_program_python("configureStrelkaSomaticWorkflow.py"),
           os.path.realpath(utils.which("configureStrelkaSomaticWorkflow.py"))]
    cur_bed = get_region_bed(region, [paired.tumor_data, paired.normal_data], out_file)
    if cur_bed:
        cmd += ["--referenceFasta=%s" % ref_file,
                "--callRegions=%s" % cur_bed,
                "--runDir=%s" % tx_work_dir,
                "--normalBam=%s" % paired.normal_bam, "--tumorBam=%s" % paired.tumor_bam]
        if _is_targeted_region(cur_bed, paired.tumor_data):
            cmd += ["--targeted"]
        do.run(cmd, "Configure Strelka2 germline calling: %s" % paired.tumor_name)
        return os.path.join(tx_work_dir, "runWorkflow.py")

def _tumor_normal_genotypes(ref, alt, info, fname, coords):
    """Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.

    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention

    ref: The REF allele from a VCF line
    alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    info: The VCF INFO field
    fname, coords: not currently used, for debugging purposes
    """
    known_names = set(["het", "hom", "ref", "conflict"])
    def name_to_gt(val):
        if val.lower() == "het":
            return "0/1"
        elif val.lower() == "hom":
            return "1/1"
        elif val.lower() in set(["ref", "conflict"]):
            return "0/0"
        else:
            # Non-standard representations, het is our best imperfect representation
            # print(fname, coords, ref, alt, info, val)
            return "0/1"
    def alleles_to_gt(val):
        gt_indices = {gt.upper(): i for i, gt in enumerate([ref] + alt)}
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts) == 0:
                tumor_gt = "0/0"
            elif 0 in tumor_gts:
                tumor_gt = "0/%s" % min([x for x in tumor_gts if x > 0])
            else:
                tumor_gt = "%s/%s" % (min(tumor_gts), max(tumor_gts))
        else:
            tumor_gt = name_to_gt(val)
        return tumor_gt
    nt_val = [x.split("=")[-1] for x in info if x.startswith("NT=")][0]
    normal_gt = name_to_gt(nt_val)
    sgt_val = [x.split("=")[-1] for x in info if x.startswith("SGT=")]
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val = sgt_val[0].split("->")[-1]
        tumor_gt = alleles_to_gt(sgt_val)
    return tumor_gt, normal_gt

def _af_annotate_and_filter(paired, items, in_file, out_file):
    """Populating FORMAT/AF, and dropping variants with AF<min_allele_fraction

    Strelka2 doesn't report exact AF for a variant, however it can be calculated as alt_counts/dp from existing fields:
    somatic
      snps:    GT:DP:FDP:SDP:SUBDP:AU:CU:GU:TU                 dp=DP                {ALT}U[0] = alt_counts(tier1,tier2)
      indels:  GT:DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50  dp=DP                TIR = alt_counts(tier1,tier2)
    germline
      snps:    GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL(:PS)       dp=sum(alt_counts)   AD = ref_count,alt_counts
      indels:  GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL(:PS)             dp=sum(alt_counts)   AD = ref_count,alt_counts
    """
    data = paired.tumor_data if paired else items[0]
    min_freq = float(utils.get_in(data["config"], ("algorithm", "min_allele_fraction"), 10)) / 100.0
    logger.debug("Filtering Strelka2 calls with allele fraction threshold of %s" % min_freq)
    ungz_out_file = "%s.vcf" % utils.splitext_plus(out_file)[0]
    if not utils.file_exists(ungz_out_file) and not utils.file_exists(ungz_out_file + ".gz"):
        with file_transaction(data, ungz_out_file) as tx_out_file:
            vcf = cyvcf2.VCF(in_file)
            vcf.add_format_to_header({
                'ID': 'AF',
                'Description': 'Allele frequency, as calculated in bcbio: AD/DP (germline), <ALT>U/DP (somatic snps), '
                               'TIR/DPI (somatic indels)',
                'Type': 'Float',
                'Number': '.'})
            vcf.add_filter_to_header({
                'ID': 'MinAF',
                'Description': 'Allele frequency is lower than %s%% ' % (min_freq*100) + (
                    '(configured in bcbio as min_allele_fraction)'
                    if utils.get_in(data["config"], ("algorithm", "min_allele_fraction"))
                    else '(default threshold in bcbio; override with min_allele_fraction in the algorithm section)')})
            w = cyvcf2.Writer(tx_out_file, vcf)
            tumor_index = vcf.samples.index(data['description'])
            for rec in vcf:
                if paired:  # somatic?
                    if rec.is_snp:  # snps?
                        alt_counts = rec.format(rec.ALT[0] + 'U')[:,0]  # {ALT}U=tier1_depth,tier2_depth
                    else:  # indels
                        alt_counts = rec.format('TIR')[:,0]  # TIR=tier1_depth,tier2_depth
                    dp = rec.format('DP')[:,0]
                elif rec.format("AD") is not None:  # germline?
                    alt_counts = rec.format('AD')[:,1:]  # AD=REF,ALT1,ALT2,...
                    dp = np.sum(rec.format('AD')[:,0:], axis=1)[:, None]
                else: # germline gVCF record
                    alt_counts, dp = (None, None)
                if dp is not None:
                    with np.errstate(divide='ignore', invalid='ignore'):  # ignore division by zero and put AF=.0
                        af = np.true_divide(alt_counts, dp)
                        af[~np.isfinite(af)] = .0  # -inf inf NaN -> .0
                    rec.set_format('AF', af)
                    if paired and np.all(af[tumor_index] < min_freq):
                        vcfutils.cyvcf_add_filter(rec, 'MinAF')
                w.write_record(rec)
            w.close()
    return vcfutils.bgzip_and_index(ungz_out_file, data["config"])

def _postprocess_somatic(in_file, paired):
    """Post-process somatic calls to provide standard output.

    - Converts SGT and NT into standard VCF GT fields
    - Replace generic TUMOR NORMAL names in VCF with sample names.
    """
    out_file = in_file.replace(".vcf.gz", "-fixed.vcf")
    if not utils.file_exists(out_file) and not utils.file_exists(out_file + ".gz"):
        with file_transaction(paired.tumor_data, out_file) as tx_out_file:
            with utils.open_gzipsafe(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    added_gt = False
                    normal_index, tumor_index = (None, None)
                    for line in in_handle:
                        if line.startswith("##FORMAT") and not added_gt:
                            added_gt = True
                            out_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                            out_handle.write(line)
                        elif line.startswith("#CHROM"):
                            assert added_gt
                            parts = line.strip().split("\t")
                            normal_index = parts.index("NORMAL")
                            tumor_index = parts.index("TUMOR")
                            line = line.replace("NORMAL", paired.normal_name).replace("TUMOR", paired.tumor_name)
                            out_handle.write(line)
                        elif line.startswith("#"):
                            out_handle.write(line)
                        else:
                            parts = line.rstrip().split("\t")
                            tumor_gt, normal_gt = _tumor_normal_genotypes(parts[3], parts[4].split(","),
                                                                          parts[7].split(";"), in_file, parts[:2])
                            parts[8] = "GT:%s" % parts[8]
                            parts[normal_index] = "%s:%s" % (normal_gt, parts[normal_index])
                            parts[tumor_index] = "%s:%s" % (tumor_gt, parts[tumor_index])
                            out_handle.write("\t".join(parts) + "\n")
    return vcfutils.bgzip_and_index(out_file, paired.tumor_data["config"])

def _run_somatic(paired, ref_file, assoc_files, region, out_file, work_dir):
    if not utils.file_exists(out_file):
        with file_transaction(paired.tumor_data, work_dir) as tx_work_dir:
            workflow_file = _configure_somatic(paired, ref_file, region, out_file, tx_work_dir)
            if workflow_file:
                has_variants = True
                _run_workflow(paired.tumor_data, workflow_file, tx_work_dir)
            else:
                has_variants = False
                vcfutils.write_empty_vcf(out_file, paired.tumor_data["config"],
                                         [dd.get_sample_name(d) for d in [paired.tumor_data, paired.normal_data]])
        if has_variants:
            var_dir = os.path.join(work_dir, "results", "variants")
            vcfutils.combine_variant_files([_postprocess_somatic(os.path.join(var_dir, f), paired)
                                            for f in ["somatic.snvs.vcf.gz", "somatic.indels.vcf.gz"]],
                                           out_file, ref_file, paired.tumor_data["config"], region=region)
    return out_file

def _run_workflow(data, workflow_file, work_dir):
    """Run Strelka2 analysis inside prepared workflow directory.
    """
    utils.remove_safe(os.path.join(work_dir, "workspace"))
    cmd = [utils.get_program_python("configureStrelkaGermlineWorkflow.py"),
           workflow_file, "-m", "local", "-j", dd.get_num_cores(data), "--quiet"]
    do.run(cmd, "Run Strelka2: %s" % dd.get_sample_name(data))
    utils.remove_safe(os.path.join(work_dir, "workspace"))

# ## Joint calling

def run_gvcfgenotyper(data, orig_region, vrn_files, out_file):
    """Merge strelka2 and Illumina compatible gVCFs with gvcfgenotyper.

    https://github.com/Illumina/gvcfgenotyper

    Also need to explore GLnexus (https://github.com/dnanexus-rnd/GLnexus)
    """
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            regions = _find_gvcf_blocks(vrn_files[0], bamprep.region_to_gatk(orig_region),
                                        os.path.dirname(tx_out_file))
            if len(regions) == 1:
                _run_gvcfgenotyper(data, regions[0], vrn_files, tx_out_file)
            else:
                split_outs = [_run_gvcfgenotyper(data, r, vrn_files,
                                                 "%s-%s.vcf.gz" % (utils.splitext_plus(out_file)[0],
                                                                   r.replace(":", "_").replace("-", "_")))
                              for r in regions]
                vcfutils.concat_variant_files(split_outs, tx_out_file, regions,
                                              dd.get_ref_file(data), data["config"])
    return vcfutils.bgzip_and_index(out_file, data["config"])

def _run_gvcfgenotyper(data, region, vrn_files, out_file):
    """Run gvcfgenotyper on a single gVCF region in input file.
    """
    if not utils.file_exists(out_file):
        with file_transaction(data, out_file) as tx_out_file:
            input_file = "%s-inputs.txt" % utils.splitext_plus(tx_out_file)[0]
            with open(input_file, "w") as out_handle:
                out_handle.write("%s\n" % "\n".join(vrn_files))
            cmd = ["gvcfgenotyper", "-f", dd.get_ref_file(data), "-l", input_file,
                   "-r", region, "-O", "z", "-o", tx_out_file]
            do.run(cmd, "gvcfgenotyper: %s %s" % (dd.get_sample_name(data), region))
    return out_file

def _find_gvcf_blocks(vcf_file, region, tmp_dir):
    """Retrieve gVCF blocks within our current evaluation region.

    gvcfgenotyper does not support calling larger regions with individual
    coverage blocks, so we split our big region into potentially multiple.
    """
    region_file = os.path.join(tmp_dir, "cur_region.bed")
    with open(region_file, "w") as out_handle:
        chrom, coords = region.split(":")
        start, end = coords.split("-")
        out_handle.write("\t".join([chrom, start, end]) + "\n")
    final_file = os.path.join(tmp_dir, "split_regions.bed")
    cmd = "gvcf_regions.py {vcf_file} | bedtools intersect -a - -b {region_file} > {final_file}"
    do.run(cmd.format(**locals()))
    regions = []
    with open(final_file) as in_handle:
        for line in in_handle:
            chrom, start, end = line.strip().split("\t")
            regions.append("%s:%s-%s" % (chrom, start, end))
    return regions
