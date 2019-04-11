"""Loads configurations from .yaml files and expands environment variables.
"""
import copy
import collections
import glob
import math
import os
import pprint
import sys
import yaml

import toolz as tz

import six


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
            for key, val in custom.items():
                config["algorithm"][key] = val
    # apply any algorithm details specified with the lane
    for key, val in lane_info.get("algorithm", {}).items():
        config["algorithm"][key] = val
    # apply any resource details specified with the lane
    for prog, pkvs in lane_info.get("resources", {}).items():
        if prog not in config["resources"]:
            config["resources"][prog] = {}
        for key, val in pkvs.items():
            config["resources"][prog][key] = val
    return config

# ## Retrieval functions

def load_system_config(config_file=None, work_dir=None, allow_missing=False):
    """Load bcbio_system.yaml configuration file, handling standard defaults.

    Looks for configuration file in default location within
    final base directory from a standard installation. Handles both standard
    installs (galaxy/bcbio_system.yaml) and docker installs (config/bcbio_system.yaml).
    """
    docker_config = _get_docker_config()
    if config_file is None:
        config_file = "bcbio_system.yaml"
    if not os.path.exists(config_file):
        base_dir = get_base_installdir()
        test_config = os.path.join(base_dir, "galaxy", config_file)
        if os.path.exists(test_config):
            config_file = test_config
        elif allow_missing:
            config_file = None
        else:
            raise ValueError("Could not find input system configuration file %s, "
                             "including inside standard directory %s" %
                             (config_file, os.path.join(base_dir, "galaxy")))
    config = load_config(config_file) if config_file else {}
    if docker_config:
        assert work_dir is not None, "Need working directory to merge docker config"
        config_file = os.path.join(work_dir, "%s-merged%s" % os.path.splitext(os.path.basename(config_file)))
        config = _merge_system_configs(config, docker_config, config_file)
    if "algorithm" not in config:
        config["algorithm"] = {}
    config["bcbio_system"] = config_file
    return config, config_file

def get_base_installdir(cmd=sys.executable):
    return os.path.normpath(os.path.join(os.path.realpath(cmd), os.pardir, os.pardir, os.pardir))

def _merge_system_configs(host_config, container_config, out_file=None):
    """Create a merged system configuration from external and internal specification.
    """
    out = copy.deepcopy(container_config)
    for k, v in host_config.items():
        if k in set(["galaxy_config"]):
            out[k] = v
        elif k == "resources":
            for pname, resources in v.items():
                if not isinstance(resources, dict) and pname not in out[k]:
                    out[k][pname] = resources
                else:
                    for rname, rval in resources.items():
                        if (rname in set(["cores", "jvm_opts", "memory"])
                              or pname in set(["gatk", "mutect"])):
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
            yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
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
        return _update_config(args, _update_resources, allow_missing=True)

def load_config(config_file):
    """Load YAML config file, replacing environmental variables.
    """
    with open(config_file) as in_handle:
        config = yaml.safe_load(in_handle)
    config = _expand_paths(config)
    if 'resources' not in config:
        config['resources'] = {}
    # lowercase resource names, the preferred way to specify, for back-compatibility
    newr = {}
    for k, v in config["resources"].items():
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
    return tz.get_in(["resources", name], config,
                     tz.get_in(["resources", "default"], config, {}))

def get_program(name, config, ptype="cmd", default=None):
    """Retrieve program information from the configuration.

    This handles back compatible location specification in input
    YAML. The preferred location for program information is in
    `resources` but the older `program` tag is also supported.
    """
    # support taking in the data dictionary
    config = config.get("config", config)
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
        return _get_program_cmd(name, pconfig, config, default)
    elif ptype == "dir":
        return _get_program_dir(name, pconfig)
    else:
        raise ValueError("Don't understand program type: %s" % ptype)

def _get_check_program_cmd(fn):
    def wrap(name, pconfig, config, default):
        is_ok = lambda f: os.path.isfile(f) and os.access(f, os.X_OK)
        bcbio_system = config.get("bcbio_system", None)
        if bcbio_system:
            system_bcbio_path = os.path.join(os.path.dirname(bcbio_system),
                                             os.pardir, "anaconda", "bin", name)
            if is_ok(system_bcbio_path):
                return system_bcbio_path
        # support bioconda installed programs
        if is_ok(os.path.join(os.path.dirname(sys.executable), name)):
            return (os.path.join(os.path.dirname(sys.executable), name))
        # find system bioconda installed programs if using private code install
        program = expand_path(fn(name, pconfig, config, default))
        if is_ok(program):
            return program
        # search the PATH now
        for adir in os.environ['PATH'].split(":"):
            if is_ok(os.path.join(adir, program)):
                return os.path.join(adir, program)
        raise CmdNotFound(" ".join(map(repr, (fn.__name__ if six.PY3 else fn.func_name, name, pconfig, default))))
    return wrap

@_get_check_program_cmd
def _get_program_cmd(name, pconfig, config, default):
    """Retrieve commandline of a program.
    """
    if pconfig is None:
        return name
    elif isinstance(pconfig, six.string_types):
        return pconfig
    elif "cmd" in pconfig:
        return pconfig["cmd"]
    elif default is not None:
        return default
    else:
        return name

def _get_program_dir(name, config):
    """Retrieve directory for a program (local installs/java jars).
    """
    if config is None:
        raise ValueError("Could not find directory in config for %s" % name)
    elif isinstance(config, six.string_types):
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

def is_std_config_arg(x):
    return isinstance(x, dict) and "algorithm" in x and "resources" in x and "files" not in x

def is_nested_config_arg(x):
    return isinstance(x, dict) and "config" in x and is_std_config_arg(x["config"])

def get_algorithm_config(xs):
    """Flexibly extract algorithm configuration for a sample from any function arguments.
    """
    if isinstance(xs, dict):
        xs = [xs]
    for x in xs:
        if is_std_config_arg(x):
            return x["algorithm"]
        elif is_nested_config_arg(x):
            return x["config"]["algorithm"]
        elif isinstance(x, (list, tuple)) and is_nested_config_arg(x[0]):
            return x[0]["config"]["algorithm"]
    raise ValueError("Did not find algorithm configuration in items: {0}"
                     .format(pprint.pformat(xs)))

def get_dataarg(args):
    """Retrieve the world 'data' argument from a set of input parameters.
    """
    for i, arg in enumerate(args):
        if is_nested_config_arg(arg):
            return i, arg
        elif is_std_config_arg(arg):
            return i, {"config": arg}
        elif isinstance(arg, (list, tuple)) and is_nested_config_arg(arg[0]):
            return i, arg[0]
    raise ValueError("Did not find configuration or data object in arguments: %s" % args)

def add_cores_to_config(args, cores_per_job, parallel=None):
    """Add information about available cores for a job to configuration.
    Ugly hack to update core information in a configuration dictionary.
    """
    def _update_cores(config):
        config["algorithm"]["num_cores"] = int(cores_per_job)
        if parallel:
            parallel.pop("view", None)
            config["parallel"] = parallel
        return config
    return _update_config(args, _update_cores)

def _update_config(args, update_fn, allow_missing=False):
    """Update configuration, nested in argument list, with the provided update function.
    """
    new_i = None
    for i, arg in enumerate(args):
        if (is_std_config_arg(arg) or is_nested_config_arg(arg) or
              (isinstance(arg, (list, tuple)) and is_nested_config_arg(arg[0]))):
            new_i = i
            break
    if new_i is None:
        if allow_missing:
            return args
        else:
            raise ValueError("Could not find configuration in args: %s" % str(args))

    new_arg = args[new_i]
    if is_nested_config_arg(new_arg):
        new_arg["config"] = update_fn(copy.deepcopy(new_arg["config"]))
    elif is_std_config_arg(new_arg):
        new_arg = update_fn(copy.deepcopy(new_arg))
    elif isinstance(arg, (list, tuple)) and is_nested_config_arg(new_arg[0]):
        new_arg_first = new_arg[0]
        new_arg_first["config"] = update_fn(copy.deepcopy(new_arg_first["config"]))
        new_arg = [new_arg_first] + new_arg[1:]
    else:
        raise ValueError("Unexpected configuration dictionary: %s" % new_arg)
    args = list(args)[:]
    args[new_i] = new_arg
    return args

def convert_to_bytes(mem_str):
    """Convert a memory specification, potentially with M or G, into bytes.
    """
    if str(mem_str)[-1].upper().endswith("G"):
        return int(round(float(mem_str[:-1]) * 1024 * 1024))
    elif str(mem_str)[-1].upper().endswith("M"):
        return int(round(float(mem_str[:-1]) * 1024))
    else:
        return int(round(float(mem_str)))

def adjust_cores_to_mb_target(target_mb, mem_str, cores):
    """Scale core usage to match a Mb/core target.

    Useful for memory dependent programs where we don't have control
    over memory usage so need to scale cores.
    """
    cur_mb = convert_to_bytes(mem_str) / 1024.0
    scale = target_mb / cur_mb
    if scale >= 1:
        return cores
    else:
        return max(1, int(math.ceil(scale * cores)))

def adjust_memory(val, magnitude, direction="increase", out_modifier="", maximum=None):
    """Adjust memory based on number of cores utilized.
    """
    modifier = val[-1:]
    amount = float(val[:-1])
    if direction == "decrease":
        new_amount = amount / float(magnitude)
        # dealing with a specifier like 1G, need to scale to Mb
        if new_amount < 1 or (out_modifier.upper().startswith("M") and modifier.upper().startswith("G")):
            if modifier.upper().startswith("G"):
                new_amount = (amount * 1024) / magnitude
                modifier = "M" + modifier[1:]
            else:
                raise ValueError("Unexpected decrease in memory: %s by %s" % (val, magnitude))
        amount = int(new_amount)
    elif direction == "increase" and magnitude > 1:
        # for increases with multiple cores, leave small percentage of
        # memory for system to maintain process running resource and
        # avoid OOM killers
        adjuster = 0.91
        amount = int(math.ceil(amount * (adjuster * magnitude)))
    if out_modifier.upper().startswith("G") and modifier.upper().startswith("M"):
        modifier = out_modifier
        amount = int(math.floor(amount / 1024.0))
    if out_modifier.upper().startswith("M") and modifier.upper().startswith("G"):
        modifier = out_modifier
        modifier = int(amount * 1024)
    if maximum:
        max_modifier = maximum[-1]
        max_amount = float(maximum[:-1])
        if modifier.upper() == "G" and max_modifier.upper() == "M":
            max_amount = max_amount / 1024.0
        elif modifier.upper() == "M" and max_modifier.upper() == "G":
            max_amount = max_amount * 1024.0
        amount = min([amount, max_amount])
    return "{amount}{modifier}".format(amount=int(math.floor(amount)), modifier=modifier)

def adjust_opts(in_opts, config):
    """Establish JVM opts, adjusting memory for the context if needed.

    This allows using less or more memory for highly parallel or multicore
    supporting processes, respectively.
    """
    memory_adjust = config["algorithm"].get("memory_adjust", {})
    out_opts = []
    for opt in in_opts:
        if opt.startswith("-Xmx") or (opt.startswith("-Xms") and memory_adjust.get("direction") == "decrease"):
            arg = opt[:4]
            opt = "{arg}{val}".format(arg=arg,
                                      val=adjust_memory(opt[4:],
                                                        memory_adjust.get("magnitude", 1),
                                                        memory_adjust.get("direction"),
                                                        maximum=memory_adjust.get("maximum")))
        out_opts.append(opt)
    return out_opts

# specific program usage

def use_vqsr(algs, call_file=None):
    """Processing uses GATK's Variant Quality Score Recalibration.
    """
    from bcbio.variation import vcfutils
    vqsr_callers = set(["gatk", "gatk-haplotype"])
    vqsr_sample_thresh = 50
    vqsr_supported = collections.defaultdict(int)
    coverage_intervals = set([])
    for alg in algs:
        callers = alg.get("variantcaller")
        if isinstance(callers, six.string_types):
            callers = [callers]
        if not callers:  # no variant calling, no VQSR
            continue
        if "vqsr" in (alg.get("tools_off") or []):  # VQSR turned off
            continue
        for c in callers:
            if c in vqsr_callers:
                if "vqsr" in (alg.get("tools_on") or []):  # VQSR turned on:
                    vqsr_supported[c] += 1
                    coverage_intervals.add("genome")
                # Do not try VQSR for gVCF inputs
                elif call_file and vcfutils.is_gvcf_file(call_file):
                    pass
                else:
                    coverage_intervals.add(alg.get("coverage_interval", "exome").lower())
                    vqsr_supported[c] += 1
    if len(vqsr_supported) > 0:
        num_samples = max(vqsr_supported.values())
        if "genome" in coverage_intervals or num_samples >= vqsr_sample_thresh:
            return True
    return False

def use_snpeff(algs):
    """Processing uses snpEff. Avoids memory requirements if not used.
    """
    return any(alg.get("effects", "snpeff") == "snpeff" and alg.get("variantcaller") for alg in algs)

def use_bcbio_variation_recall(algs):
    """Processing uses bcbio-variation-recall. Avoids core requirement if not used.
    """
    for alg in algs:
        jointcaller = alg.get("jointcaller", [])
        if not isinstance(jointcaller, (tuple, list)):
            jointcaller = [jointcaller]
        for caller in jointcaller:
            if caller not in set(["gatk-haplotype-joint", None, False]):
                return True
    return False

## functions for navigating through the standard galaxy directory of files

def get_rRNA_interval(genome_dir):
    return os.path.join(genome_dir, "rnaseq", "rRNA.interval_list")

def get_transcript_refflat(genome_dir):
    return os.path.join(genome_dir, "rnaseq", "ref-transcripts.refFlat")

def get_rRNA_sequence(genome_dir):
    return os.path.join(genome_dir, "rnaseq", "rRNA.fa")

def program_installed(program, data):
    """
    returns True if the path to a program can be found
    """
    try:
        path = get_program(program, data)
    except CmdNotFound:
        return False
    return True
