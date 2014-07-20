"""
functions to access the data dictionary in a clearer way
"""

import toolz as tz
from bcbio.utils import file_exists
from bcbio.log import logger
from bcbio.pipeline import config_utils
import sys

LOOKUPS = {
    "gtf_file": {"keys": ['genome_resources', 'rnaseq', 'transcripts'],
                 "checker": file_exists},
    "work_dir": {"keys": ['dirs', 'work']},
    "sample_name": {"keys": ['rgnames', 'sample']},
    "strandedness": {"keys": ['config', 'algorithm', 'strandedness']},
    "work_bam": {"keys": ["work_bam"]},
    "dexseq_gff": {"keys": ['genome_resources', 'rnaseq', 'dexseq']}
}

def getter(keys):
    def lookup(config, default=None):
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
    checker = v.get('checker', None)
    _g["get_" + k] = getter(keys)
    _g["set_" + k] = setter(keys, checker)

