"""Provide support for MuTect and other paired analysis tools."""

from distutils.version import LooseVersion
import os

import toolz as tz

from bcbio import bam, broad, utils
from bcbio.utils import file_exists, get_in, open_gzipsafe
from bcbio.distributed.transaction import file_transaction
from bcbio.heterogeneity import chromhacks
from bcbio.pipeline import config_utils
from bcbio.pipeline import datadict as dd
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.variation import bamprep, bedutils, gatk, vcfutils, scalpel
from bcbio.variation.realign import has_aligned_reads
from bcbio.variation.vcfutils import bgzip_and_index
from bcbio.log import logger

_PASS_EXCEPTIONS = set(["java.lang.RuntimeException: "
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!",
                        "java.lang.IllegalArgumentException: "
                        "Comparison method violates its general contract!"])

def _check_mutect_version(broad_runner):
    mutect_version = broad_runner.get_mutect_version()
    try:
        assert mutect_version is not None
    except AssertionError:
        logger.warn("WARNING")
        logger.warn("MuTect version could not be determined from jar file. "
                    "Please ensure you are using at least version 1.1.5, "
                    "as versions 1.1.4 and lower have known issues.")
        logger.warn("Proceeding but assuming correct version 1.1.5.")
    else:
        try:
            assert LooseVersion(mutect_version) >= LooseVersion("1.1.5")
        except AssertionError:
            message = ("MuTect 1.1.4 and lower is known to have incompatibilities "
                       "with Java < 7, and this may lead to problems in analyses. "
                       "Please use MuTect 1.1.5 or higher (note that it requires "
                       "Java 7).")
            raise ValueError(message)

def _config_params(base_config, assoc_files, region, out_file, items):
    """Add parameters based on configuration variables, associated files and genomic regions.
    """
    params = []
    dbsnp = assoc_files.get("dbsnp")
    if dbsnp:
        params += ["--dbsnp", dbsnp]
    cosmic = assoc_files.get("cosmic")
    if cosmic:
        params += ["--cosmic", cosmic]
    variant_regions = bedutils.population_variant_regions(items)
    region = subset_variant_regions(variant_regions, region, out_file)
    if region:
        params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule",
                   "INTERSECTION"]
    # set low frequency calling parameter if adjusted
    # to set other MuTect parameters on contamination, pass options to resources for mutect
    # --fraction_contamination --minimum_normal_allele_fraction
    min_af = tz.get_in(["algorithm", "min_allele_fraction"], base_config)
    if min_af:
        params += ["--minimum_mutation_cell_fraction", "%.2f" % (min_af / 100.0)]
    resources = config_utils.get_resources("mutect", base_config)
    if resources.get("options") is not None:
        params += [str(x) for x in resources.get("options", [])]
    # Output quality scores
    if "--enable_qscore_output" not in params:
        params.append("--enable_qscore_output")
    # drf not currently supported in MuTect to turn off duplicateread filter
    # params += gatk.standard_cl_params(items)
    return params

def _mutect_call_prep(align_bams, items, ref_file, assoc_files,
                       region=None, out_file=None):
    """Preparation work for MuTect.
    """
    base_config = items[0]["config"]
    broad_runner = broad.runner_from_path("picard", base_config)
    broad_runner.run_fn("picard_index_ref", ref_file)

    broad_runner = broad.runner_from_config(base_config, "mutect")
    _check_mutect_version(broad_runner)
    for x in align_bams:
        bam.index(x, base_config)

    paired = vcfutils.get_paired_bams(align_bams, items)
    if not paired:
        raise ValueError("Specified MuTect calling but 'tumor' phenotype not present in batch\n"
                         "https://bcbio-nextgen.readthedocs.org/en/latest/contents/"
                         "pipelines.html#cancer-variant-calling\n"
                         "for samples: %s" % ", " .join([dd.get_sample_name(x) for x in items]))
    params = ["-R", ref_file, "-T", "MuTect", "-U", "ALLOW_N_CIGAR_READS"]
    params += ["--read_filter", "NotPrimaryAlignment"]
    params += ["-I:tumor", paired.tumor_bam]
    params += ["--tumor_sample_name", paired.tumor_name]
    if paired.normal_bam is not None:
        params += ["-I:normal", paired.normal_bam]
        params += ["--normal_sample_name", paired.normal_name]
    if paired.normal_panel is not None:
        params += ["--normal_panel", paired.normal_panel]
    params += _config_params(base_config, assoc_files, region, out_file, items)
    return broad_runner, params

def mutect_caller(align_bams, items, ref_file, assoc_files, region=None,
                  out_file=None):
    """Run the MuTect paired analysis algorithm.
    """
    config = items[0]["config"]
    if out_file is None:
        out_file = "%s-paired-variants.vcf.gz" % os.path.splitext(align_bams[0])[0]
    if not file_exists(out_file):
        base_config = items[0]["config"]
        broad_runner = broad.runner_from_config(base_config, "mutect")
        out_file_mutect = (out_file.replace(".vcf", "-mutect.vcf")
                           if "vcf" in out_file else out_file + "-mutect.vcf")
        broad_runner, params = \
            _mutect_call_prep(align_bams, items, ref_file, assoc_files,
                                   region, out_file_mutect)
        if (not isinstance(region, (list, tuple)) and
              not all(has_aligned_reads(x, region) for x in align_bams)):
                vcfutils.write_empty_vcf(out_file)
                return
        out_file_orig = "%s-orig%s" % utils.splitext_plus(out_file_mutect)
        if not file_exists(out_file_orig):
            with file_transaction(config, out_file_orig) as tx_out_file:
                # Rationale: MuTect writes another table to stdout, which we don't need
                params += ["--vcf", tx_out_file, "-o", os.devnull]
                broad_runner.run_mutect(params)
        is_paired = "-I:normal" in params
        if not utils.file_uptodate(out_file_mutect, out_file_orig):
            out_file_mutect = _fix_mutect_output(out_file_orig, config, out_file_mutect, is_paired)
        indelcaller = vcfutils.get_indelcaller(base_config)
        if ("scalpel" in indelcaller.lower() and region and isinstance(region, (tuple, list))
              and chromhacks.is_autosomal_or_sex(region[0])):
            # Scalpel InDels
            out_file_indels = (out_file.replace(".vcf", "-somaticIndels.vcf")
                               if "vcf" in out_file else out_file + "-somaticIndels.vcf")
            if scalpel.is_installed(items[0]["config"]):
                if not is_paired:
                    vcfutils.check_paired_problems(items)
                    scalpel._run_scalpel_caller(align_bams, items, ref_file, assoc_files,
                                                region=region, out_file=out_file_indels)
                else:
                    scalpel._run_scalpel_paired(align_bams, items, ref_file, assoc_files,
                                                region=region, out_file=out_file_indels)
                out_file = vcfutils.combine_variant_files(orig_files=[out_file_mutect, out_file_indels],
                                                          out_file=out_file,
                                                          ref_file=items[0]["sam_ref"],
                                                          config=items[0]["config"],
                                                          region=region)
            else:
                utils.symlink_plus(out_file_mutect, out_file)
        elif "pindel" in indelcaller.lower():
            from bcbio.structural import pindel
            out_file_indels = (out_file.replace(".vcf", "-somaticIndels.vcf")
                               if "vcf" in out_file else out_file + "-somaticIndels.vcf")
            if pindel.is_installed(items[0]["config"]):
                pindel._run_tumor_pindel_caller(align_bams, items, ref_file, assoc_files, region=region,
                                          out_file=out_file_indels)
                out_file = vcfutils.combine_variant_files(orig_files=[out_file_mutect, out_file_indels],
                                                          out_file=out_file,
                                                          ref_file=ref_file,
                                                          config=items[0]["config"],
                                                          region=region)
            else:
                utils.symlink_plus(out_file_mutect, out_file)
        elif (("somaticindeldetector" in indelcaller.lower() or "sid" in indelcaller.lower())
              and "appistry" in broad_runner.get_mutect_version()):
            # SomaticIndelDetector InDels
            out_file_indels = (out_file.replace(".vcf", "-somaticIndels.vcf")
                               if "vcf" in out_file else out_file + "-somaticIndels.vcf")
            params_indels = _SID_call_prep(align_bams, items, ref_file, assoc_files,
                                           region, out_file_indels)
            with file_transaction(config, out_file_indels) as tx_out_file:
                params_indels += ["-o", tx_out_file]
                broad_runner.run_mutect(params_indels)
            out_file = vcfutils.combine_variant_files(orig_files=[out_file_mutect, out_file_indels],
                                                      out_file=out_file,
                                                      ref_file=items[0]["sam_ref"],
                                                      config=items[0]["config"],
                                                      region=region)
        else:
            utils.symlink_plus(out_file_mutect, out_file)
    return out_file

def _SID_call_prep(align_bams, items, ref_file, assoc_files, region=None, out_file=None):
    """Preparation work for SomaticIndelDetector.
    """
    base_config = items[0]["config"]
    for x in align_bams:
        bam.index(x, base_config)

    params = ["-R", ref_file, "-T", "SomaticIndelDetector", "-U", "ALLOW_N_CIGAR_READS"]
    # Limit per base read start count to between 200-10000, i.e. from any base
    # can no more 10000 new reads begin.
    # Further, limit maxNumberOfReads accordingly, otherwise SID discards
    # windows for high coverage panels.
    paired = vcfutils.get_paired_bams(align_bams, items)
    params += ["--read_filter", "NotPrimaryAlignment"]
    params += ["-I:tumor", paired.tumor_bam]
    min_af = float(get_in(paired.tumor_config, ("algorithm", "min_allele_fraction"), 10)) / 100.0
    if paired.normal_bam is not None:
        params += ["-I:normal", paired.normal_bam]
        # notice there must be at least 4 reads of coverage in normal
        params += ["--filter_expressions", "T_COV<6||N_COV<4||T_INDEL_F<%s||T_INDEL_CF<0.7" % min_af]
    else:
        params += ["--unpaired"]
        params += ["--filter_expressions", "COV<6||INDEL_F<%s||INDEL_CF<0.7" % min_af]
    if region:
        params += ["-L", bamprep.region_to_gatk(region), "--interval_set_rule",
                   "INTERSECTION"]
    return params

def _fix_mutect_output(orig_file, config, out_file, is_paired):
    """Adjust MuTect output to match other callers.

    - Rename allelic fraction field in mutect output from FA to FREQ to standarize with other tools
    - Remove extra 'none' samples introduced when calling tumor-only samples
    """
    out_file_noc = out_file.replace(".vcf.gz", ".vcf")
    none_index = -1
    with file_transaction(config, out_file_noc) as tx_out_file:
        with open_gzipsafe(orig_file) as in_handle:
            with open(tx_out_file, 'w') as out_handle:
                for line in in_handle:
                    if not is_paired and line.startswith("#CHROM"):
                        parts = line.rstrip().split("\t")
                        none_index = parts.index("none")
                        del parts[none_index]
                        line = "\t".join(parts) + "\n"
                    elif line.startswith("##FORMAT=<ID=FA"):
                        line = line.replace("=FA", "=FREQ")
                    elif not line.startswith("#"):
                        if none_index > 0:
                            parts = line.rstrip().split("\t")
                            del parts[none_index]
                            line = "\t".join(parts) + "\n"
                        line = line.replace("FA", "FREQ")
                    out_handle.write(line)
    return bgzip_and_index(out_file_noc, config)
