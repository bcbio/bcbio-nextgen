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
    if "algorithm" not in config:
        config["algorithm"] = {}
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

# ## Retrieval functions

def load_system_config(config_file, work_dir=None):
    """Load bcbio_system.yaml configuration file, handling standard defaults.

    Looks for configuration file in default location within
    final base directory from a standard installation. Handles both standard
    installs (galaxy/bcbio_system.yaml) and docker installs (config/bcbio_system.yaml).
    """
    docker_config = _get_docker_config()
    if not os.path.exists(config_file):
        base_dir = get_base_installdir()
        test_config = os.path.join(base_dir, "galaxy", config_file)
        if os.path.exists(test_config):
            config_file = test_config
        else:
            raise ValueError("Could not find input system configuration file %s, "
                             "including inside standard directory %s" %
                             (config_file, os.path.join(base_dir, "galaxy")))
    config = load_config(config_file)
    if docker_config:
        assert work_dir is not None, "Need working directory to merge docker config"
        config_file = os.path.join(work_dir, "%s-merged%s" % os.path.splitext(os.path.basename(config_file)))
        config = _merge_system_configs(config, docker_config, config_file)
    if "algorithm" not in config:
        config["algorithm"] = {}
    config["bcbio_system"] = config_file
    return config, config_file

def get_base_installdir():
    return os.path.normpath(os.path.join(os.path.realpath(sys.executable), os.pardir, os.pardir, os.pardir))

def _merge_system_configs(host_config, container_config, out_file=None):
    """Create a merged system configuration from external and internal specification.
    """
    out = copy.deepcopy(container_config)
    for k, v in host_config.iteritems():
        if k in set(["galaxy_config"]):
            out[k] = v
        elif k == "resources":
            for pname, resources in v.iteritems():
                if not isinstance(resources, dict) and pname not in out[k]:
                    out[k][pname] = resources
                else:
                    for rname, rval in resources.iteritems():
                        if rname in set(["cores", "jvm_opts", "memory"]):
                            if pname not in out[k]:
                                out[k][pname] = {}
                            out[k][pname][rname] = rval
    # Ensure final file is relocatable by mapping back to reference directory
    if "bcbio_system" in out and ("galaxy_config" not in out or not os.path.isabs(out["galaxy_config"])):
        out["galaxy_config"] = os.path.normpath(os.path.join(os.path.dirname(out["bcbio_system"]),
                                                             os.pardir, "galaxy",
                                                             "universe_wsgi.ini"))
    if out_file:
        with open(out_file, "w") as out_handle:
            yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out

def _get_docker_config():
    base_dir = get_base_installdir()
    docker_configfile = os.path.join(base_dir, "config", "bcbio_system.yaml")
    if os.path.exists(docker_configfile):
        return load_config(docker_configfile)

def merge_resources(args):
    """Merge docker local resources and global resource specification in a set of arguments.

    Finds the `data` object within passed arguments and updates the resources
    from a local docker configuration if present.
    """
    docker_config = _get_docker_config()
    if not docker_config:
        return args
    else:
        def _update_resources(config):
            config["resources"] = _merge_system_configs(config, docker_config)["resources"]
            return config
        return _update_config(args, _update_resources)

def load_config(config_file):
    """Load YAML config file, replacing environmental variables.
    """
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    config = _expand_paths(config)
    # lowercase resource names, the preferred way to specify, for back-compatibility
    newr = {}
    for k, v in config["resources"].iteritems():
        if k.lower() != k:
            newr[k.lower()] = v
    config["resources"].update(newr)
    return config

def _expand_paths(config):
    for field, setting in config.items():
        if isinstance(config[field], dict):
            config[field] = _expand_paths(config[field])
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
            if not key in pconfig:
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
    elif "cmd" in config:
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
    elif "dir" in config:
        return expand_path(config["dir"])
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

# ## Retrieval and update to configuration from arguments

def _dictdissoc(orig, k):
    """Imitates immutability: create a new dictionary with the key dropped.
    """
    v = orig.pop(k, None)
    new = copy.deepcopy(orig)
    orig[k] = v
    return new

def is_std_config_arg(x):
    return isinstance(x, dict) and "algorithm" in x and "resources" in x and not "files" in x

def is_nested_config_arg(x):
    return isinstance(x, dict) and "config" in x and is_std_config_arg(x["config"])

def get_algorithm_config(xs):
    """Flexibly extract algorithm configuration for a sample from any function arguments.
    """
    for x in xs:
        if is_std_config_arg(x):
            return x["algorithm"]
        elif is_nested_config_arg(x):
            return x["config"]["algorithm"]
    raise ValueError("Did not find algorithm configuration in items: {0}"
                     .format(xs))

def add_cores_to_config(args, cores_per_job, parallel=None):
    """Add information about available cores for a job to configuration.
    Ugly hack to update core information in a configuration dictionary.
    """
    def _update_cores(config):
        config["algorithm"]["num_cores"] = int(cores_per_job)
        if parallel:
            config["parallel"] = _dictdissoc(parallel, "view")
        return config
    return _update_config(args, _update_cores)

def _update_config(args, update_fn):
    """Update configuration, nested in argument list, with the provided update function.
    """
    new_i = None
    for i, arg in enumerate(args):
        if is_std_config_arg(arg) or is_nested_config_arg(arg):
            new_i = i
            break
    if new_i is None:
        raise ValueError("Could not find configuration in args: %s" % args)

    new_arg = copy.deepcopy(args[new_i])
    if is_nested_config_arg(new_arg):
        new_arg["config"] = update_fn(new_arg["config"])
    elif is_std_config_arg(new_arg):
        new_arg = update_fn(new_arg)
    else:
        raise ValueError("Unexpected configuration dictionary: %s" % new_arg)
    args = list(args)[:]
    args[new_i] = new_arg
    return args

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
        elif not callers:  # no variant calling, no VQSR
            continue
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
