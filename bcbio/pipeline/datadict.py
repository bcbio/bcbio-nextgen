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
    "num_cores": {"keys": ['config', 'algorithm', 'num_cores'],
                  "default": 1},
    "priority_regions": {"keys": ['config', 'algorithm', 'priority_regions']},
    "gtf_file": {"keys": ['genome_resources', 'rnaseq', 'transcripts'],
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
    "square_vcf": {"keys": ['square_vcf']},
    "ploidy": {"keys": ['config', 'algorithm', 'ploidy'], "default": 2},
    "gender": {"keys": ["metadata", "sex"], "default": ""},
    "batch": {"keys": ["metadata", "batch"]},
    "hetcaller": {"keys": ["config", "algorithm", "hetcaller"]},
    "variantcaller": {"keys": ['config', 'algorithm', 'variantcaller']},
    "work_bam": {"keys": ["work_bam"]},
    "count_file": {"keys": ["count_file"]},
    "combined_counts": {"keys": ["combined_counts"]},
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
    "variation_resources": {"keys": ["genome_resources", "variation"], "default": {}},
    "qsig_file": {"keys": ['genome_resources', 'variation', 'qsignature'],
                  "checker": file_exists},
    "mixup_check": {"keys": ["config", "algorithm", "mixup_check"],
                    "default": False},
    "cufflinks_dir": {"keys": ['cufflinks_dir']},
    "rsem": {"keys": ["config", "algorithm", "rsem"], "default": False},
    "transcriptome_bam": {"keys": ["transcriptome_bam"]},
    "fpkm_isoform": {"keys": ["fpkm_isoform"]},
    "fpkm": {"keys": ["fpkm"]},
    "galaxy_dir": {"keys": ["dirs", "galaxy"]},
    "assembled_gtf": {"keys": ["assembled_gtf"]},
    "assemble_transcripts": {"keys": ["config", "algorithm", "assemble_transcripts"],
                             "default": False},
    "oncofuse_file": {"keys": ["oncofuse_file"]},
    "split_bam": {"keys": ["split_bam"]},
    "vrn_file": {"keys": ["vrn_file"]},
    "variant_regions": {"keys": ["config", "algorithm", "variant_regions"]},
    "callable_regions": {"keys": ["regions", "callable"]},
    "offtarget_stats": {"keys": ["regions", "offtarget_stats"]},
    "sample_callable": {"keys": ["regions", "sample_callable"]},
    "coverage_interval": {"keys": ["config", "algorithm", "coverage_interval"]},
    "coverage_regions": {"keys": ["config", "algorithm", "coverage"]},
    "deduped_bam": {"keys": ["deduped_bam"]},
    "align_bam": {"keys": ["align_bam"]},
    "tools_off": {"keys": ["config", "algorithm", "tools_off"], "default": []},
    "tools_on": {"keys": ["config", "algorithm", "tools_on"], "default": []},
}

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

def sample_data_iterator(samples):
    """
    for a list of samples, return the data dictionary of each sample
    """
    for sample in samples:
        yield sample[0]
