"""Loads configurations from .yaml files and expands environment variables.

The configuration yaml has the structure

galaxy_config:
program:
	program1:
	program2:
algorithm:
	setting1:
	setting2:
log_dir:
store_dir:
store_host:
analysis:
	config_file:
	towig_script:
distributed:
	rabbitmq_vhost:
custom_algorithms:
	setting1:
	setting2:

galaxy_config, program and analysis supports
environment variables.
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

def get_program(name, config, ptype="cmd"):
    """Retrieve program information from the configuration.

    This handles back compatible location specification in input
    YAML. The preferred location for program information is in
    `resources` but the older `program` tag is also supported.
    """
    try:
        pconfig = config.get("resources", {})[name]
    except KeyError:
        pconfig = config.get("program", {}).get("name", None)
    if ptype == "cmd":
        return _get_program_cmd(name, pconfig)
    elif ptype == "dir":
        return _get_program_dir(name, pconfig)
    else:
        raise ValueError("Don't understand program type: %s" % ptype)

def _get_program_cmd(name, config):
    """Retrieve commandline of a program.
    """
    if config.has_key("cmd"):
        return config["cmd"]
    elif isinstance(config, basestring):
        return config
    else:
        return name

def _get_program_dir(name, config):
    """Retrieve directory for a program (local installs/java jars).
    """
    if config.has_key("dir"):
        return config["dir"]
    elif isinstance(config, basestring):
        return config
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
