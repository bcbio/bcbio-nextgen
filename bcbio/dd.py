"""
functions to access the data dictionary in a clearer way
"""

import toolz as tz

def make_access_methods(cls):
    for k, v in cls.LOOKUPS.items():
        setattr(cls, "get_" + k, getter(cls, v))
        setattr(cls, "set_" + k, setter(cls, v))
    return cls

def getter(cls, keys):
    def lookup(cls, config, default=None):
        return tz.get_in(keys, config, default)
    return lookup

def setter(cls, keys):
    def update(cls, config, value):
        return tz.update_in(config, keys, lambda x: value, default=value)
    return update

@make_access_methods
class ConfigParser(object):
    LOOKUPS = {"gtf_file": ['genome_resources', 'rnaseq', 'transcripts']}

    def __init__(self):
        pass

dd = ConfigParser()
