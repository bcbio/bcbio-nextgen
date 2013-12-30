"""Expose bcbio-nextgen installation and upgrade capabilities.
"""
import json
import os

import tornado.web
import yaml

import bcbio.install

def add_defaults(dct, defaults, overwrite=False):
    for k, v in defaults.iteritems():
        if overwrite or k not in dct:
            dct[k] = v
    return dct

def json_to_args(x, defaults=None, mandatory=None):
    dictargs = json.loads(x)
    if defaults:
        dictargs = add_defaults(dictargs, defaults)
    if mandatory:
        dictargs = add_defaults(dictargs, mandatory, overwrite=True)
    class Args: pass
    args = Args()
    for k, v in dictargs.iteritems():
        setattr(args, k, v)
    return args

def _install_bcbio_system(datadir):
    """Install limited bcbio_system.yaml file for setting core and memory usage.

    Adds any non-specific programs to the exposed bcbio_system.yaml file.
    """
    expose = set(["memory", "cores", "jvm_opts"])
    base_file = os.path.join(datadir, "config", "bcbio_system.yaml")
    expose_file = os.path.join(datadir, "galaxy", "bcbio_system.yaml")
    with open(base_file) as in_handle:
        config = yaml.load(in_handle)
    if os.path.exists(expose_file):
        with open(expose_file) as in_handle:
            expose_config = yaml.load(in_handle)
    else:
        expose_config = {"resources": {}}
    for pname, vals in config["resources"].iteritems():
        expose_vals = {}
        for k, v in vals.iteritems():
            if k in expose:
                expose_vals[k] = v
        if len(expose_vals) > 0 and pname not in expose_config["resources"]:
            expose_config["resources"][pname] = expose_vals
    with open(expose_file, "w") as out_handle:
        yaml.dump(expose_config, out_handle, default_flow_style=False, allow_unicode=False)
    return expose_file

def get_handler(args):
    """Enable upgrade of data in place on external system from inside container.
    """
    defaults = {"genomes": [], "aligners": [], "install_data": False,
                "toolplus": []}
    mandatory = {"sudo": False, "upgrade": "skip", "tools": False, "tooldir": None,
                 "tooldist": "minimal", "isolate": True, "distribution": ""}
    class InstallHandler(tornado.web.RequestHandler):
        def get(self):
            iargs = json_to_args(self.get_argument("args", "{}"), defaults,
                                 mandatory)
            iargs = bcbio.install.upgrade_bcbio(iargs)
            _install_bcbio_system(iargs.datadir)
    return InstallHandler
