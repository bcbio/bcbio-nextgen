"""Loads configurations from .yaml files and expands environment variables.
"""
import copy
import collections
import glob
import math
import os
import sys
import yaml


class CmdNotFound(Exception):
    pass

# ## Generalized configuration

def update_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    name_remaps = {"variant": ["SNP calling", "variant", "variant2"],
                   "SNP calling": ["SNP calling", "variant", "variant2"],
                   "variant2": ["SNP calling", "variant", "variant2"]}
    config = copy.deepcopy(config)
    base_name = lane_info.get("analysis")
    for analysis_type in name_remaps.get(base_name, [base_name]):
        custom = config.get("custom_algorithms", {}).get(analysis_type)
        if custom:
            for key, val in custom.iteritems():
                config["algorithm"][key] = val
    # apply any algorithm details specified with the lane
    for key, val in lane_info.get("algorithm", {}).iteritems():
        config["algorithm"][key] = val
    # apply any resource details specified with the lane
    for prog, pkvs in lane_info.get("resources", {}).iteritems():
        if prog not in config["resources"]:
            config["resources"][prog] = {}
        for key, val in pkvs.iteritems():
            config["resources"][prog][key] = val
    return config

def add_cached_versions(config):
    """Add version information to configuration, avoiding multiple access during parallel runs.
    """
    from bcbio import broad
    # cache GATK version in sample information to avoid multiple retrievals later
    if "gatk" in config["resources"]:
        config["resources"]["gatk"]["version"] = broad.runner_from_config(config).get_gatk_version()
    return config

# ## Retrieval functions

def load_system_config(config_file):
    """Load bcbio_system.yaml configuration file, handling standard defaults.

    Looks for configuration file in default location within
    final base directory from a standard installation.
    """
    if not os.path.exists(config_file):
        base_dir = os.path.normpath(os.path.join(os.path.realpath(sys.executable), os.pardir, os.pardir, os.pardir))
        test_config = os.path.join(base_dir, "galaxy", config_file)
        if os.path.exists(test_config):
            config_file = test_config
        else:
            raise ValueError("Could not find input system configuration file %s, "
                             "including inside standard directory %s" %
                             (config_file, os.path.join(base_dir, "galaxy")))
    return load_config(config_file), config_file

def load_config(config_file):
    """Load YAML config file, replacing environmental variables.
    """
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    for field, setting in config.items():
        if isinstance(config[field], dict):
            for sub_field, sub_setting in config[field].items():
                config[field][sub_field] = expand_path(sub_setting)
        else:
            config[field] = expand_path(setting)
    return config

def expand_path(path):
    """ Combines os.path.expandvars with replacing ~ with $HOME.
    """
    try:
        return os.path.expandvars(path.replace("~", "$HOME"))
    except AttributeError:
        return path

def get_resources(name, config):
    """Retrieve resources for a program, pulling from multiple config sources.
    """
    return config.get("resources", {}).get(name, {})

def get_program(name, config, ptype="cmd", default=None):
    """Retrieve program information from the configuration.

    This handles back compatible location specification in input
    YAML. The preferred location for program information is in
    `resources` but the older `program` tag is also supported.
    """
    try:
        pconfig = config.get("resources", {})[name]
        # If have leftover old
    except KeyError:
        pconfig = {}
    old_config = config.get("program", {}).get(name, None)
    if old_config:
        for key in ["dir", "cmd"]:
            if not pconfig.has_key(key):
                pconfig[key] = old_config
    if ptype == "cmd":
        return _get_program_cmd(name, pconfig, default)
    elif ptype == "dir":
        return _get_program_dir(name, pconfig)
    else:
        raise ValueError("Don't understand program type: %s" % ptype)

def _get_check_program_cmd(fn):

    def wrap(name, config, default):
        program = expand_path(fn(name, config, default))
        is_ok = lambda f: os.path.isfile(f) and os.access(f, os.X_OK)
        if is_ok(program): return program

        for adir in os.environ['PATH'].split(":"):
            if is_ok(os.path.join(adir, program)):
                return os.path.join(adir, program)
        else:
            raise CmdNotFound(" ".join(map(repr, (fn.func_name, name, config, default))))
    return wrap

@_get_check_program_cmd
def _get_program_cmd(name, config, default):
    """Retrieve commandline of a program.
    """
    if config is None:
        return name
    elif isinstance(config, basestring):
        return config
    elif config.has_key("cmd"):
        return config["cmd"]
    elif default is not None:
        return default
    else:
        return name

def _get_program_dir(name, config):
    """Retrieve directory for a program (local installs/java jars).
    """
    if config is None:
        raise ValueError("Could not find directory in config for %s" % name)
    elif isinstance(config, basestring):
        return config
    elif config.has_key("dir"):
        return config["dir"]
    else:
        raise ValueError("Could not find directory in config for %s" % name)

def get_jar(base_name, dname):
    """Retrieve a jar in the provided directory
    """
    jars = glob.glob(os.path.join(expand_path(dname), "%s*.jar" % base_name))
    if len(jars) == 1:
        return jars[0]
    elif len(jars) > 1:
        raise ValueError("Found multiple jars for %s in %s. Need single jar: %s" %
                         (base_name, dname, jars))
    else:
        raise ValueError("Could not find java jar %s in %s" %
                         (base_name, dname))

def adjust_memory(val, magnitude, direction="increase"):
    """Adjust memory based on number of cores utilized.
    """
    modifier = val[-1:]
    amount = int(val[:-1])
    if direction == "decrease":
        new_amount = amount / magnitude
        # dealing with a specifier like 1G, need to scale to Mb
        if new_amount < 1:
            if modifier.upper().startswith("G"):
                new_amount = (amount * 1024) / magnitude
                modifier = "M" + modifier[1:]
            else:
                raise ValueError("Unexpected decrease in memory: %s by %s" % (val, magnitude))
        amount = new_amount
    elif direction == "increase":
        # for increases with multiple cores, leave small percentage of
        # memory for system to maintain process running resource and
        # avoid OOM killers
        adjuster = 0.91
        amount = int(math.ceil(amount * (adjuster * magnitude)))
    return "{amount}{modifier}".format(amount=amount, modifier=modifier)

def adjust_opts(in_opts, config):
    """Establish JVM opts, adjusting memory for the context if needed.

    This allows using less or more memory for highly parallel or multicore
    supporting processes, respectively.
    """
    memory_adjust = config["algorithm"].get("memory_adjust", {})
    out_opts = []
    for opt in in_opts:
        if opt.startswith(("-Xmx", "-Xms")):
            arg = opt[:4]
            opt = "{arg}{val}".format(arg=arg,
                                      val=adjust_memory(opt[4:],
                                                        memory_adjust.get("magnitude", 1),
                                                        memory_adjust.get("direction")))
        out_opts.append(opt)
    return out_opts

# specific program usage

def _get_coverage_params(alg):
    Cov = collections.namedtuple("Cov", ["interval", "depth"])
    return Cov(alg.get("coverage_interval", "exome").lower(),
               alg.get("coverage_depth", "high").lower())

def use_vqsr(algs):
    """Processing uses GATK's Variant Quality Score Recalibration.
    """
    for alg in algs:
        cov = _get_coverage_params(alg)
        callers = alg.get("variantcaller", "gatk")
        if isinstance(callers, basestring):
            callers = [callers]
        vqsr_supported_caller = False
        for c in callers:
            if c in ["gatk", "gatk-haplotype"]:
                vqsr_supported_caller = True
                break
        if (cov.interval not in ["regional", "exome"] and cov.depth != "low"
            and vqsr_supported_caller):
            return True
    return False


## functions for navigating through the standard galaxy directory of files

def get_transcript_gtf(genome_dir):
    out_file = os.path.join(genome_dir, "rnaseq", "ref-transcripts.gtf")
    return out_file

def get_rRNA_interval(genome_dir):
    return os.path.join(genome_dir, "rnaseq", "rRNA.interval")

def get_transcript_refflat(genome_dir):
    return os.path.join(genome_dir, "rnaseq", "ref-transcripts.refFlat")

def get_rRNA_sequence(genome_dir):
    return os.path.join(genome_dir, "rnaseq", "rRNA.fa")
