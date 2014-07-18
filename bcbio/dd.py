"""
functions to access the data dictionary in a clearer way
"""

import toolz as tz
from bcbio.utils import file_exists
from bcbio.log import logger
import sys

def make_access_methods(cls):
    for k, v in cls.LOOKUPS.items():
        keys = v['keys']
        checker = v.get('checker', None)
        setattr(cls, "get_" + k, getter(cls, keys))
        setattr(cls, "set_" + k, setter(cls, keys, checker))
    return cls

def getter(cls, keys):
    def lookup(cls, config, default=None):
        return tz.get_in(keys, config, default)
    return lookup

def setter(cls, keys, checker):
    print checker
    def update(cls, config, value):
        if checker and not checker(value):
            logger.error("%s fails check %s." % (value, checker))
            sys.exit(1)
        return tz.update_in(config, keys, lambda x: value, default=value)
    return update

@make_access_methods
class ConfigParser(object):
    LOOKUPS = {"gtf_file": {"keys": ['genome_resources', 'rnaseq', 'transcripts'],
                            "checker": file_exists}}

    def __init__(self):
        pass

dd = ConfigParser()
