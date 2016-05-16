"""
functions to access the data dictionary in a clearer way
"""

import os
import toolz as tz
from bcbio.utils import file_exists
from bcbio.log import logger
import sys

LOOKUPS = {
    "config": {"keys": ['config']},
    "tmp_dir": {"keys": ['config', 'resources', 'tmp', 'dir']},
    "num_cores": {"keys": ['config', 'algorithm', 'num_cores'],
                  "default": 1},
    "svprioritize": {"keys": ['config', 'algorithm', 'svprioritize']},
    "problem_region_dir": {"keys": ["config", "algorithm", "problem_region_dir"]},
    "genome_build": {"keys": ["genome_build"]},
    "gtf_file": {"keys": ['genome_resources', 'rnaseq', 'transcripts'],
                 "checker": file_exists},
    "transcriptome_fasta": {"keys": ["config", "algorithm", "transcriptome_fasta"],
                            "default": None},
    "tx2gene": {"keys": ["tx2gene"]},
    "srna_gtf_file": {"keys": ['genome_resources', 'srnaseq', 'srna-transcripts'],
                      "checker": file_exists},
    "srna_trna_file": {"keys": ['genome_resources', 'srnaseq', 'trna-fasta'],
                      "checker": file_exists},
    "mirdeep2_file": {"keys": ['genome_resources', 'srnaseq', 'mirdeep2-fasta'],
                      "checker": file_exists},
    "mirbase_hairpin": {"keys": ['genome_resources', 'srnaseq', 'mirbase-hairpin'],
                      "checker": file_exists},
    "mirbase_mature": {"keys": ['genome_resources', 'srnaseq', 'mirbase-mature'],
                      "checker": file_exists},
    "gene_bed": {"keys": ['genome_resources', 'rnaseq', 'gene_bed'],
                 "checker": file_exists},
    "work_dir": {"keys": ['dirs', 'work']},
    "sam_ref": {"keys": ["sam_ref"]},
    "disambiguate": {"keys": ["config", "algorithm", "disambiguate"],
                     "default": False},
    "lane": {"keys": ["rgnames", "lane"]},
    "cores": {"keys": ["config", "algorithm", "num_cores"], "default": 1},
    "sample_name": {"keys": ['rgnames', 'sample']},
    "strandedness": {"keys": ['config', 'algorithm', 'strandedness'],
                     "default": "unstranded"},
    "analysis": {"keys": ["analysis"]},
    "square_vcf": {"keys": ['square_vcf']},
    "ploidy": {"keys": ['config', 'algorithm', 'ploidy'], "default": 2},
    "gender": {"keys": ["metadata", "sex"], "default": ""},
    "batch": {"keys": ["metadata", "batch"]},
    "mark_duplicates": {"keys": ["config", "algorithm", "mark_duplicates"], "default": True},
    "phenotype": {"keys": ["metadata", "phenotype"], "default": ""},
    "hetcaller": {"keys": ["config", "algorithm", "hetcaller"]},
    "variantcaller": {"keys": ['config', 'algorithm', 'variantcaller']},
    "peakcaller": {"keys": ['config', 'algorithm', 'peakcaller'], "default": []},
    "work_bam": {"keys": ["work_bam"]},
    "count_file": {"keys": ["count_file"]},
    "combined_counts": {"keys": ["combined_counts"]},
    "mirna_counts": {"keys": ["mirna_counts"]},
    "isomir_counts": {"keys": ["isomir_counts"]},
    "novel_mirna_counts": {"keys": ["novel_mirna_counts"]},
    "novel_isomir_counts": {"keys": ["novel_isomir_counts"]},
    "annotated_combined_counts": {"keys": ["annotated_combined_counts"]},
    "ref_file": {"keys": ["reference", "fasta", "base"]},
    "dexseq_gff": {"keys": ['genome_resources', 'rnaseq', 'dexseq']},
    "combined_fpkm": {"keys": ['combined_fpkm']},
    "combined_fpkm_isoform": {"keys": ['combined_fpkm_isoform']},
    "express_fpkm": {"keys": ['express_fpkm']},
    "express_tpm": {"keys": ['express_tpm']},
    "express_counts": {"keys": ['express_counts']},
    "isoform_to_gene": {"keys": ['isoform_to_gene']},
    "fusion_mode": {"keys": ['config', 'algorithm', 'fusion_mode']},
    "dexseq_counts": {"keys": ['dexseq_counts']},
    "description": {"keys": ['description']},
    "aligner": {"keys": ['config', 'algorithm', 'aligner']},
    "platform": {"keys": ['config', 'algorithm', 'platform'],
                 "default": "illumina"},
    "quality_format": {"keys": ['config', 'algorithm', 'quality_format'],
                       "default": "standard"},
    "adapters": {"keys": ['config', 'algorithm', 'adapters'],
                 "default": []},
    "custom_trim": {"keys": ['config', 'algorithm', 'custom_trim'],
                 "default": []},
    "species": {"keys": ['config', 'algorithm', 'species'],
                 "default": None},
    "trim_reads": {"keys": ['config', 'algorithm', 'trim_reads'],
                 "default": None},
    "variation_resources": {"keys": ["genome_resources", "variation"], "default": {}},
    "qsig_file": {"keys": ['genome_resources', 'variation', 'qsignature'],
                  "checker": file_exists},
    "mixup_check": {"keys": ["config", "algorithm", "mixup_check"],
                    "default": False},
    "cufflinks_dir": {"keys": ['cufflinks_dir']},
    "rsem": {"keys": ["config", "algorithm", "rsem"], "default": False},
    "transcriptome_align": {"keys": ["config", "algorithm", "transcriptome_align"],
                            "default": False},
    "expression_caller": {"keys": ["config", "algorithm", "expression_caller"],
                          "default": []},
    "transcriptome_bam": {"keys": ["transcriptome_bam"]},
    "fpkm_isoform": {"keys": ["fpkm_isoform"]},
    "fpkm": {"keys": ["fpkm"]},
    "galaxy_dir": {"keys": ["dirs", "galaxy"]},
    "assembled_gtf": {"keys": ["assembled_gtf"], "default": []},
    "merged_gtf": {"keys": ["merged_gtf"], "default": None},
    "transcript_assembler": {"keys": ["config", "algorithm", "transcript_assembler"],
                             "default": []},
    "oncofuse_file": {"keys": ["oncofuse_file"]},
    "split_bam": {"keys": ["split_bam"]},
    "vrn_file": {"keys": ["vrn_file"]},
    "variant_regions": {"keys": ["config", "algorithm", "variant_regions"]},
    "variant_regions_merged": {"keys": ["config", "algorithm", "variant_regions_merged"]},
    "callable_regions": {"keys": ["regions", "callable"]},
    "callable_min_size": {"keys": ["config", "algorithm", "callable_min_size"],
                          "default": 1000000},
    "min_allele_fraction": {"keys": ["config", "algorithm", "min_allele_fraction"]},
    "offtarget_stats": {"keys": ["regions", "offtarget_stats"]},
    "salmon": {"keys": ["salmon"]},
    "umi_type": {"keys": ["config", "algorithm", "umi_type"]},
    "cellular_barcodes": {"keys": ["config", "algorithm", "cellular_barcodes"],
                          "default": []},
    "minimum_barcode_depth": {"keys": ["config", "algorithm", "minimum_barcode_depth"],
                              "default": 100000},
    "salmon_dir": {"keys": ["salmon_dir"]},
    "sailfish": {"keys": ["sailfish"]},
    "sailfish_dir": {"keys": ["sailfish_dir"]},
    "sailfish_tidy": {"keys": ["sailfish_tidy"]},
    "sailfish_transcript_tpm": {"keys": ["sailfish_transcript_tpm"]},
    "sailfish_gene_tpm": {"keys": ["sailfish_gene_tpm"]},
    "sample_callable": {"keys": ["regions", "sample_callable"]},
    "coverage_interval": {"keys": ["config", "algorithm", "coverage_interval"]},
    "coverage": {"keys": ["config", "algorithm", "coverage"]},
    "report": {"keys": ["config", "algorithm", "report"]},
    "coverage_depth_min": {"keys": ["config", "algorithm", "coverage_depth_min"],
                           "default": 4},
    "joint_group_size": {"keys": ["config", "algorithm", "joint_group_size"],
                         "default": 200},
    "coverage": {"keys": ["config", "algorithm", "coverage"]},
    "deduped_bam": {"keys": ["deduped_bam"]},
    "align_bam": {"keys": ["align_bam"]},
    "align_prep_method": {"keys": ["config", "algorithm", "align_prep_method"], "default": "grabix"},
    "tools_off": {"keys": ["config", "algorithm", "tools_off"], "default": []},
    "tools_on": {"keys": ["config", "algorithm", "tools_on"], "default": []},
    "cwl_reporting": {"keys": ["config", "algorithm", "cwl_reporting"]},
}

def get_batches(data):
    batches = get_batch(data)
    if batches:
        if not isinstance(batches, (list, tuple)):
            batches = [batches]
        return batches

def get_input_sequence_files(data, default=None):
    """
    returns the input sequencing files, these can be single or paired FASTQ
    files or BAM files
    """
    if "files" not in data:
        file1, file2 = None, None
    elif len(data["files"]) == 2:
        file1, file2 = data["files"]
    else:
        assert len(data["files"]) == 1, data["files"]
        file1, file2 = data["files"][0], None
    return file1, file2

def get_dexseq_gff(config, default=None):
    """
    some older versions of the genomes have the DEXseq gff file as
    gff instead of gff3, so this handles that by looking for either one
    """
    dexseq_gff = tz.get_in(tz.get_in(['dexseq_gff', 'keys'], LOOKUPS, {}),
                           config, None)
    if not dexseq_gff:
        return None
    gtf_file = get_gtf_file(config)
    if gtf_file:
        base_dir = os.path.dirname(gtf_file)
    else:
        base_dir = os.path.dirname(dexseq_gff)
    base, _ = os.path.splitext(dexseq_gff)
    gff_file = os.path.join(base_dir, base + ".gff")
    if file_exists(gff_file):
        return gff_file
    gtf_file = os.path.join(base_dir, base + ".gff3")
    if file_exists(gtf_file):
        return gtf_file
    else:
        return None

def getter(keys, global_default=None):
    def lookup(config, default=None):
        default = global_default if not default else default
        return tz.get_in(keys, config, default)
    return lookup

def setter(keys, checker):
    def update(config, value):
        if checker and not checker(value):
            logger.error("%s fails check %s." % (value, checker))
            sys.exit(1)
        return tz.update_in(config, keys, lambda x: value, default=value)
    return update

def is_setter(keys):
    def present(config):
        try:
            value = tz.get_in(keys, config, no_default=True)
        except:
            value = False
        return True if value else False
    return present

"""
generate the getter and setter functions but don't override any explicitly
defined
"""
_g = globals()
for k, v in LOOKUPS.items():
    keys = v['keys']
    getter_fn = 'get_' + k
    if getter_fn not in _g:
        _g["get_" + k] = getter(keys, v.get('default', None))
    setter_fn = 'set_' + k
    if setter_fn not in _g:
        _g["set_" + k] = setter(keys, v.get('checker', None))
    is_setter_fn = "is_set" + k
    if is_setter_fn not in _g:
        _g["is_set_" + k] = is_setter(keys)

def sample_data_iterator(samples):
    """
    for a list of samples, return the data dictionary of each sample
    """
    for sample in samples:
        yield sample[0]

def get_in_samples(samples, fn):
    """
    for a list of samples, return the value of a global option
    """
    for sample in samples:
        if fn(sample[0], None):
            return fn(sample[0])
    return None
