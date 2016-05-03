"""calling using Pindel

http://gmt.genome.wustl.edu/packages/pindel/
"""

from __future__ import print_function
import os
import time
import itertools
import shutil
from bcbio import bam, utils, broad
from bcbio.distributed.transaction import file_transaction, tx_tmpdir
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions, remove_lcr_regions
from bcbio.variation.vcfutils import bgzip_and_index, get_paired_bams
from bcbio.variation import annotation
from bcbio.provenance import do


def _pindel_options(items, config, out_file, region, tmp_path):
    """parse pindel options. Add region to cmd.
    :param items: (dict) information from yaml
    :param config: (dict) information from yaml (items[0]['config'])
    :param region: (str or tupple) region to analyze
    :param tmp_path: (str) temporal folder
    :returns: (list) options for pindel
    """
    variant_regions = utils.get_in(config, ("algorithm", "variant_regions"))
    target = subset_variant_regions(variant_regions, region, out_file, items)
    opts = ""
    if target:
        if isinstance(target, basestring) and os.path.isfile(target):
            target_bed = target
        else:
            target_bed = os.path.join(tmp_path, "tmp.bed")
            with file_transaction(config, target_bed) as tx_tmp_bed:
                if not isinstance(region, (list, tuple)):
                    message = ("Region must be a tuple - something odd just happened")
                    raise ValueError(message)
                chrom, start, end = region
                with open(tx_tmp_bed, "w") as out_handle:
                    print("%s\t%s\t%s" % (chrom, start, end), file=out_handle)
        opts = "-j " + remove_lcr_regions(target_bed, items)
    return opts


def is_installed(config):
    """Check for pindel installation on machine.
    :param config: (dict) information from yaml(items[0]['config'])
    :returns: (boolean) if pindel is installed
    """
    try:
        config_utils.get_program("pindel2vcf", config)
        config_utils.get_program("pindel", config)
        return True
    except config_utils.CmdNotFound:
        return False


def _run_tumor_pindel_caller(align_bams, items, ref_file, assoc_files,
                             region=None, out_file=None):
    """Detect indels with pindel in tumor/[normal] analysis.
    Only attempts to detect small insertion/deletions and not larger structural events.
    :param align_bam: (list) bam files
    :param items: (dict) information from yaml
    :param ref_file: (str) genome in fasta format
    :param assoc_file: (dict) files for annotation
    :param region: (str or tupple) region to analyze
    :param out_file: (str) final vcf file
    :returns: (str) final vcf file
    """
    config = items[0]["config"]
    paired = get_paired_bams(align_bams, items)
    if out_file is None:
        out_file = "%s-indels.vcf" % os.path.splitext(align_bams[0])[0]
    paired_bam = [paired.tumor_bam]
    paired_name = [paired.tumor_name]
    if paired.normal_bam:
        paired_bam.append(paired.normal_bam)
        paired_name.append(paired.normal_name)
    if not utils.file_exists(out_file):
        with tx_tmpdir(config) as tmp_path:
            for align_bam in align_bams:
                bam.index(align_bam, config)
            root_pindel = os.path.join(tmp_path, "pindelroot")
            pindel = config_utils.get_program("pindel", config)
            opts = _pindel_options(items, config, out_file, region, tmp_path)
            tmp_input = _create_tmp_input(paired_bam, paired_name, tmp_path, config)
            cmd = ("{pindel} -f {ref_file} -i {tmp_input} -o {root_pindel} " +
                   "{opts} --max_range_index 2 --IndelCorrection "
                   "--report_breakpoints false --report_interchromosomal_events false")
            do.run(cmd.format(**locals()), "Genotyping with pindel", {})
            out_file = _create_vcf(root_pindel, out_file, ref_file,
                                   items, paired)
    ann_file = annotation.annotate_nongatk_vcf(out_file, align_bams,
                                               assoc_files.get("dbsnp"),
                                               ref_file, config)
    return ann_file


def _create_tmp_input(input_bams, names, tmp_path, config):
    """Create input file for pindel. tab file: bam file, insert size, name
    :param input_bams: (list) bam files
    :param names: (list) names of samples
    :param tmp_path: (str) temporal dir
    :param config: (dict) information from yaml file(itmes[0]['config'])
    :returns: (str) input file for pindel
    """
    tmp_input = os.path.join(tmp_path, "pindel.txt")
    with open(tmp_input, 'w') as out_handle:
        for bam_file, name in itertools.izip(input_bams, names):
            print("%s\t%s\t%s\n" % (bam_file, 250, name), file=out_handle)
    return tmp_input


def _create_vcf(root_file, out_file, reference, items, paired=None):
    """use pindel2vcf to create vcf file from pindel format
    :param root_file: (str) prefix for pindel
    :param out_file: (str) final vcf file
    :param reference (str) genome in fasta format
    :param items: (dics) information fro yaml file
    :param paired: (tupple) bam files and name of tumor/normal samples
    :returns: (str) final vcf file
    """
    config = items[0]["config"]
    name_ref = items[0]["genome_build"]
    date = time.strftime("%Y%m%d")
    if not utils.file_exists(out_file):
        pindel2vcf = config_utils.get_program("pindel2vcf", config)
        vcf_file = out_file.replace(".gz", "")
        with file_transaction(items[0], vcf_file) as tx_out_file:
            cmd = ("{pindel2vcf} --gatk_compatible -P {root_file} -r {reference} -R {name_ref} "
                   "-d {date} -v {tx_out_file} --compact_output_limit 15")
            do.run(cmd.format(**locals()), "Converting to vcf", {})
            if paired.normal_name:
                _filter_paired(paired.tumor_name, paired.normal_name,
                               tx_out_file, reference, items[0])
        out_file = bgzip_and_index(vcf_file, config)
    return out_file


def _filter_paired(tumor, normal, out_file, reference, data):
    """filter paired vcf file with GATK
    :param    tumor: (str) sample name for tumor
    :param    normal: (str) sample name for normal
    :param    out_file: (str) final vcf file
    :param    reference: (str) genome in fasta format
    :param    data: (dict) information from yaml file(items[0])
    :returns: (str) name of final vcf file
    """
    in_file = utils.splitext_plus(out_file)[0] + "-tmp.vcf"
    shutil.move(out_file, in_file)
    config = data["config"]
    with file_transaction(data, out_file) as tx_out_file:
        params = ["-T", "SomaticPindelFilter", "-V", in_file, "-o",
                  tx_out_file, "-TID", tumor, "-NID", normal, "-R", reference]
        jvm_opts = broad.get_gatk_framework_opts(config)
        cmd = [config_utils.get_program("gatk-framework", config)] + jvm_opts + params
        do.run(cmd, "Filter pindel variants")
    return out_file
