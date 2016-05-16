"""GATK variant calling -- MuTect2.
"""
from distutils.version import LooseVersion
import os
import sys

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import annotation, bamprep, bedutils, gatk, vcfutils, ploidy

def _add_tumor_params(paired):
    """Add tumor/normal BAM input parameters to command line.
    """
    params = []
    if not paired:
        raise ValueError("Specified MuTect2 calling but 'tumor' phenotype not present in batch\n"
                         "https://bcbio-nextgen.readthedocs.org/en/latest/contents/"
                         "pipelines.html#cancer-variant-calling\n"
                         "for samples: %s" % ", " .join([dd.get_sample_name(x) for x in items]))
    params += ["-I:tumor", paired.tumor_bam]
    if paired.normal_bam is not None:
        params += ["-I:normal", paired.normal_bam]
    if paired.normal_panel is not None:
        params += ["--normal_panel", paired.normal_panel]
    return params

def _add_region_params(region, out_file, items):
    """Add parameters for selecting by region to command line.
    """
    params = []
    variant_regions = bedutils.population_variant_regions(items)
    region = subset_variant_regions(variant_regions, region, out_file, items)
    if region:
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
        _prep_inputs(align_bams, ref_file, items)
        with file_transaction(items[0], out_file) as tx_out_file:
            params = ["-T", "MuTect2",
                      "-R", ref_file,
                      "--annotation", "ClippingRankSumTest",
                      "--annotation", "DepthPerSampleHC"]
            for a in annotation.get_gatk_annotations(items[0]["config"]):
                params += ["--annotation", a]
            paired = vcfutils.get_paired_bams(align_bams, items)
            params += _add_tumor_params(paired)
            params += _add_region_params(region, out_file, items)
            params += _add_assoc_params(assoc_files)
            params += ["-ploidy", str(ploidy.get_ploidy(items, region))]
            resources = config_utils.get_resources("mutect2", items[0]["config"])
            if "options" in resources:
                params += [str(x) for x in resources.get("options", [])]
            broad_runner = broad.runner_from_config(items[0]["config"])
            assert LooseVersion(broad_runner.gatk_major_version()) >= LooseVersion("3.5"), \
                "Require full version of GATK 3.5+ for mutect2 calling"
            broad_runner.new_resources("mutect2")
            gatk_cmd = " ".join(broad_runner.cl_gatk(params, os.path.dirname(tx_out_file)))
            pp_cmd = _post_process_cl(paired)
            cmd = "{gatk_cmd} | {pp_cmd} | bgzip -c > {tx_out_file}"
            do.run(cmd.format(**locals()), "MuTect2")
    out_file = vcfutils.bgzip_and_index(out_file, items[0]["config"])
    return out_file

def _post_process_cl(paired):
    py_cl = os.path.join(os.path.dirname(sys.executable), "py")
    normal_name = paired.normal_name or ""
    tumor_name = paired.tumor_name or ""
    cmd = ("{py_cl} -x 'bcbio.variation.mutect2.fix_mutect2_output(x,"
           """ "{tumor_name}", "{normal_name}")' """)
    return cmd.format(**locals())

def fix_mutect2_output(line, tumor_name, normal_name):
    """Provide corrective fixes for mutect2 output.
    """
    if line.startswith("#CHROM") and (tumor_name or normal_name):
        parts = line.split("\t")
        base_header = parts[:9]
        old_samples = parts[9:]
        if len(old_samples) > 0:
            mapping = {"NORMAL": normal_name, "TUMOR": tumor_name,
                       tumor_name: tumor_name, normal_name: normal_name}
            samples = [mapping[sample_name] for sample_name in old_samples]
            samples = [x for x in samples if x]
            assert len(old_samples) == len(samples)
            line = "\t".join(base_header + samples)
    return line
