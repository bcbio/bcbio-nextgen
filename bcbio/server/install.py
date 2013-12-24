"""Expose bcbio-nextgen installation and upgrade capabilities.
"""
import json
import tornado.web
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

def get_handler(args):
    """Enable upgrade of data in place on external system.
    """
    defaults = {"genomes": [], "aligners": [], "install_data": False}
    mandatory = {"sudo": False, "upgrade": "skip", "tools": False, "tooldir": None,
                 "toolplus": [], "tooldist": "minimal", "isolate": True, "distribution": ""}
    class InstallHandler(tornado.web.RequestHandler):
        def get(self):
            iargs = json_to_args(self.get_argument("args", "{}"), defaults,
                                 mandatory)
            bcbio.install.upgrade_bcbio(iargs)
    return InstallHandler
