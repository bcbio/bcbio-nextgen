"""
functions to access the data dictionary in a clearer way
"""

import toolz as tz
from bcbio.utils import file_exists
from bcbio.log import logger
import sys

LOOKUPS = {
    "gtf_file": {"keys": ['genome_resources', 'rnaseq', 'transcripts'],
                 "checker": file_exists},
    "work_dir": {"keys": ['dirs', 'work']},
    "sample_name": {"keys": ['rgnames', 'sample']},
    "strandedness": {"keys": ['config', 'algorithm', 'strandedness'],
                     "default": "unstranded"},
    "work_bam": {"keys": ["work_bam"]},
    "ref_file": {"keys": ["reference", "fasta", "base"]},
    "dexseq_gff": {"keys": ['genome_resources', 'rnaseq', 'dexseq']},
    "fusion_mode": {"keys": ['config', 'algorithm', 'fusion_mode']},
    "dexseq_counts": {"keys": ['dexseq_counts']},
    "description": {"keys": ['description']},
    "quality_format": {"keys": ['config', 'algorithm', 'quality_format'],
                       "default": "standard"},
    "adapters": {"keys": ['config', 'algorithm', 'adapters'],
                 "default": []},
    "qsig_file": {"keys": ['genome_resources', 'variation', 'qsignature'],
                  "checker": file_exists}
}

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

_g = globals()
for k, v in LOOKUPS.items():
    keys = v['keys']
    _g["get_" + k] = getter(keys, v.get('default', None))
    _g["set_" + k] = setter(keys, v.get('checker', None))
