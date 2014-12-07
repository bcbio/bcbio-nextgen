"""
functions to access the data dictionary in a clearer way
"""

import os
import toolz as tz
from bcbio.utils import file_exists
from bcbio.log import logger
import sys

LOOKUPS = {
    "num_cores": {"keys": ['config', 'algorithm', 'num_cores'],
                  "default": 1},
    "gtf_file": {"keys": ['genome_resources', 'rnaseq', 'transcripts'],
                 "checker": file_exists},
    "work_dir": {"keys": ['dirs', 'work']},
    "sam_ref": {"keys": ["sam_ref"]},
    "lane": {"keys": ["rgnames", "lane"]},
    "cores": {"keys": ["config", "algorithm", "num_cores"], "default": 1},
    "sample_name": {"keys": ['rgnames', 'sample']},
    "strandedness": {"keys": ['config', 'algorithm', 'strandedness'],
                     "default": "unstranded"},
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
    "fusion_mode": {"keys": ['config', 'algorithm', 'fusion_mode']},
    "dexseq_counts": {"keys": ['dexseq_counts']},
    "description": {"keys": ['description']},
    "aligner": {"keys": ['config', 'algorithm', 'aligner']},
    "quality_format": {"keys": ['config', 'algorithm', 'quality_format'],
                       "default": "standard"},
    "adapters": {"keys": ['config', 'algorithm', 'adapters'],
                 "default": []},
    "qsig_file": {"keys": ['genome_resources', 'variation', 'qsignature'],
                  "checker": file_exists},
    "mixup_check": {"keys": ["config", "algorithm", "mixup_check"],
                    "default": False},
    "cufflinks_dir": {"keys": ['cufflinks_dir']},
    "rsem": {"keys": ["config", "algorithm", "rsem"], "default": False},
    "transcriptome_bam": {"keys": ["transcriptome_bam"]},
    "fpkm_isoform": {"keys": ["fpkm_isoform"]},
    "fpkm": {"keys": ["fpkm"]},
    "assembled_gtf": {"keys": ["assembled_gtf"]},
    "assemble_transcripts": {"keys": ["config", "algorithm", "assemble_transcripts"],
                             "default": False},
    "oncofuse_file": {"keys": ["oncofuse_file"]}
}

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
