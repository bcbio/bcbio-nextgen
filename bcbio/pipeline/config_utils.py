"""Loads configurations from .yaml files and expands environment variables.
"""
import os
import glob
import yaml

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
    resources = config.get("resources", {}).get(name, {})
    if "jvm_opts" not in resources:
        java_memory = config["algorithm"].get("java_memory", None)
        if java_memory:
            resources["jvm_opts"] = ["-Xms%s" % java_memory, "-Xmx%s" % java_memory]
    return resources

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
    jars = glob.glob(os.path.join(dname, "%s*.jar" % base_name))
    if len(jars) == 1:
        return jars[0]
    else:
        raise ValueError("Could not find java jar %s in %s: %s" % (
            base_name, dname, jars))
