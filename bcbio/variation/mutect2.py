"""GATK variant calling -- MuTect2.
"""
from distutils.version import LooseVersion
import os

import numpy as np

from bcbio import bam, broad, utils
from bcbio.bam import is_paired
from bcbio.log import logger
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bamprep, bedutils, gatk, vcfutils, ploidy

cyvcf2 = utils.LazyImport("cyvcf2")

def _add_tumor_params(paired, items, gatk_type):
    """Add tumor/normal BAM input parameters to command line.
    """
    params = []
    if not paired:
        raise ValueError("Specified MuTect2 calling but 'tumor' phenotype not present in batch\n"
                         "https://bcbio-nextgen.readthedocs.org/en/latest/contents/"
                         "pipelines.html#cancer-variant-calling\n"
                         "for samples: %s" % ", " .join([dd.get_sample_name(x) for x in items]))
    if gatk_type == "gatk4":
        params += ["-I", paired.tumor_bam]
        params += ["--tumor-sample", paired.tumor_name]
    else:
        params += ["-I:tumor", paired.tumor_bam]
    if paired.normal_bam is not None:
        if gatk_type == "gatk4":
            params += ["-I", paired.normal_bam]
            params += ["--normal-sample", paired.normal_name]
        else:
            params += ["-I:normal", paired.normal_bam]
    if paired.normal_panel is not None:
        panel_dir = utils.safe_makedir(os.path.join(dd.get_work_dir(items[0]), "mutect2", "panels"))
        normal_panel = vcfutils.bgzip_and_index(paired.normal_panel, items[0]["config"], out_dir=panel_dir)
        if gatk_type == "gatk4":
            params += ["--panel-of-normals", normal_panel]
        else:
            params += ["--normal_panel", normal_panel]
    return params

def _add_region_params(region, out_file, items, gatk_type):
    """Add parameters for selecting by region to command line.
    """
    params = []
    variant_regions = bedutils.population_variant_regions(items)
    region = subset_variant_regions(variant_regions, region, out_file, items)
    if region:
        if gatk_type == "gatk4":
            params += ["-L", bamprep.region_to_gatk(region), "--interval-set-rule", "INTERSECTION"]
        else:
            params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule", "INTERSECTION"]
    params += gatk.standard_cl_params(items)
    return params

def _add_assoc_params(assoc_files):
    params = []
    if assoc_files.get("dbsnp"):
        params += ["--dbsnp", assoc_files["dbsnp"]]
    if assoc_files.get("cosmic"):
        params += ["--cosmic", assoc_files["cosmic"]]
    return params

def _prep_inputs(align_bams, ref_file, items):
    """Ensure inputs to calling are indexed as expected.
    """
    broad_runner = broad.runner_from_path("picard", items[0]["config"])
    broad_runner.run_fn("picard_index_ref", ref_file)
    for x in align_bams:
        bam.index(x, items[0]["config"])

def mutect2_caller(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Call variation with GATK's MuTect2.

    This requires the full non open-source version of GATK 3.5+.
    """
    if out_file is None:
        out_file = "%s-variants.vcf.gz" % utils.splitext_plus(align_bams[0])[0]
    if not utils.file_exists(out_file):
        paired = vcfutils.get_paired_bams(align_bams, items)
        broad_runner = broad.runner_from_config(items[0]["config"])
        gatk_type = broad_runner.gatk_type()
        f1r2_file = None
        _prep_inputs(align_bams, ref_file, items)
        with file_transaction(items[0], out_file) as tx_out_file:
            params = ["-T", "Mutect2" if gatk_type == "gatk4" else "MuTect2",
                      "--annotation", "ClippingRankSumTest",
                      "--annotation", "DepthPerSampleHC"]
            if gatk_type == "gatk4":
                params += ["--reference", ref_file]
            else:
                params += ["-R", ref_file]
            for a in annotation.get_gatk_annotations(items[0]["config"], include_baseqranksum=False):
                params += ["--annotation", a]
            # Avoid issues with BAM CIGAR reads that GATK doesn't like
            if gatk_type == "gatk4":
                params += ["--read-validation-stringency", "LENIENT"]
            params += _add_tumor_params(paired, items, gatk_type)
            params += _add_region_params(region, out_file, items, gatk_type)


            if all(is_paired(bam) for bam in align_bams) and (
                    "mutect2_readmodel" in utils.get_in(items[0], "config",
                                                        "tools_on")):
                orientation_filter = True
            else:
                orientation_filter = False

            if gatk_type == "gatk4" and orientation_filter:
                f1r2_file = "{}-f1r2.tar.gz".format(
                    utils.splitext_plus(out_file)[0])
                params += ["--f1r2-tar-gz", f1r2_file]

            # Avoid adding dbSNP/Cosmic so they do not get fed to variant filtering algorithm
            # Not yet clear how this helps or hurts in a general case.
            #params += _add_assoc_params(assoc_files)
            resources = config_utils.get_resources("mutect2", items[0]["config"])
            if "options" in resources:
                params += [str(x) for x in resources.get("options", [])]
            assert LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.5"), \
                "Require full version of GATK 3.5+ for mutect2 calling"
            broad_runner.new_resources("mutect2")
            gatk_cmd = broad_runner.cl_gatk(params, os.path.dirname(tx_out_file))
            if gatk_type == "gatk4":

                tx_raw_prefilt_file = "%s-raw%s" % utils.splitext_plus(out_file)
                tx_raw_file = "%s-raw-filt%s" % utils.splitext_plus(tx_out_file)

                if orientation_filter:
                    tx_f1r2_file = "{}-read-orientation-model.tar.gz"
                    tx_f1r2_file = tx_f1r2_file.format(
                        utils.splitext_plus(f1r2_file)[0])
                    tx_read_orient_cmd = _mutect2_read_filter(broad_runner,
                                                              f1r2_file,
                                                              tx_f1r2_file)

                    filter_cmd = _mutect2_filter(broad_runner, items,
                                                 tx_raw_prefilt_file,
                                                 tx_raw_file, ref_file,
                                                 tx_f1r2_file)
                else:
                    filter_cmd = _mutect2_filter(broad_runner, items,
                                                 tx_raw_prefilt_file,
                                                 tx_raw_file, ref_file)
                if orientation_filter:
                    cmd = "{gatk_cmd} -O {tx_raw_prefilt_file} && {tx_read_orient_cmd} && {filter_cmd}"
                else:
                    cmd = "{gatk_cmd} -O {tx_raw_prefilt_file} && {filter_cmd}"
            else:
                tx_raw_file = "%s-raw%s" % utils.splitext_plus(tx_out_file)
                cmd = "{gatk_cmd} > {tx_raw_file}"
            do.run(cmd.format(**locals()), "MuTect2")
            out_file = _af_filter(paired.tumor_data, tx_raw_file, out_file)
    return vcfutils.bgzip_and_index(out_file, items[0]["config"])

def _mutect2_read_filter(broad_runner, in_file, out_file):

    """Calculate and apply the Mutect2 read model to compensate for
    stand bias and artefacts such as those cause by FFPE."""

    params = ["-T", "LearnReadOrientationModel", "-I", in_file, "-O",
              out_file]
    return "{}".format(broad_runner.cl_gatk(params, os.path.dirname(out_file)))

def _mutect2_filter(broad_runner, items, in_file, out_file, ref_file, orient_file=None):
    """Filter of MuTect2 calls, a separate step in GATK4.

    Includes a pre-step to avoid stats information with zero callable reads, which
    cause MuTect2 filter errors when there are calls. The sed file increases this
    to 10 reads, which matches with actually having a call in the output file.
    """
    params = ["-T", "FilterMutectCalls", "--reference", ref_file, "--variant", in_file, "--output", out_file]
    resources = config_utils.get_resources("mutect2_filter", items[0]["config"])
    if "options" in resources:
        params += [str(x) for x in resources.get("options", [])]
    if orient_file is not None:
        params += ["--ob-priors", orient_file]
    avoid_zero_callable = r"sed -i 's/callable\t0.0/callable\t10.0/' %s.stats" % in_file
    return "%s && %s" % (avoid_zero_callable, broad_runner.cl_gatk(params, os.path.dirname(out_file)))

def _af_filter(data, in_file, out_file):
    """Soft-filter variants with AF below min_allele_fraction (appends "MinAF" to FILTER)
    """
    min_freq = float(utils.get_in(data["config"], ("algorithm", "min_allele_fraction"), 10)) / 100.0
    logger.debug("Filtering MuTect2 calls with allele fraction threshold of %s" % min_freq)
    ungz_out_file = "%s.vcf" % utils.splitext_plus(out_file)[0]
    if not utils.file_exists(ungz_out_file) and not utils.file_exists(ungz_out_file + ".gz"):
        with file_transaction(data, ungz_out_file) as tx_out_file:
            vcf = cyvcf2.VCF(in_file)
            vcf.add_filter_to_header({
                'ID': 'MinAF',
                'Description': 'Allele frequency is lower than %s%% ' % (min_freq*100) + (
                    '(configured in bcbio as min_allele_fraction)'
                    if utils.get_in(data["config"], ("algorithm", "min_allele_fraction"))
                    else '(default threshold in bcbio; override with min_allele_fraction in the algorithm section)')})
            w = cyvcf2.Writer(tx_out_file, vcf)
            # GATK 3.x can produce VCFs without sample names for empty VCFs
            try:
                tumor_index = vcf.samples.index(dd.get_sample_name(data))
            except ValueError:
                tumor_index = None
            for rec in vcf:
                if tumor_index is not None and np.all(rec.format('AF')[tumor_index] < min_freq):
                    vcfutils.cyvcf_add_filter(rec, 'MinAF')
                w.write_record(rec)
            w.close()
    return vcfutils.bgzip_and_index(ungz_out_file, data["config"])
