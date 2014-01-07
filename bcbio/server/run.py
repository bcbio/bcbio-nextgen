"""Provide ability to run bcbio-nextgen workflows.
"""
from functools import wraps
import os
import StringIO
from threading import Thread
import uuid

import tornado.gen
import tornado.web
import yaml

from bcbio import utils
from bcbio.pipeline import config_utils

def run_async(func):
    """Run function in an asychronous background thread.

     http://stackoverflow.com/questions/13051591/tornado-blocking-asynchronous-requests
    """
    @wraps(func)
    def async_func(*args, **kwargs):
        func_hl = Thread(target=func, args=args, kwargs=kwargs)
        func_hl.start()
        return func_hl
    return async_func

@run_async
def run_bcbio_nextgen(**kwargs):
    from bcbio.pipeline.main import run_main
    callback = kwargs.pop("callback")
    app = kwargs.pop("app")
    run_id = str(uuid.uuid1())
    app.runmonitor.set_status(run_id, "running")
    callback(run_id)
    try:
        with utils.chdir(kwargs["work_dir"]):
            run_main(**kwargs)
    except:
        app.runmonitor.set_status(run_id, "failed")
        raise
    else:
        app.runmonitor.set_status(run_id, "finished")
    finally:
        print("Run ended: %s" % run_id)

def _merge_system_configs(host_config, container_config, work_dir):
    """Create a merged system configuration from external and internal specification.
    """
    out_file = os.path.join(work_dir, "web-bcbio_system.yaml")
    out, _ = config_utils.load_system_config(container_config)
    for k, v in host_config.iteritems():
        if k in set(["galaxy_config"]):
            out[k] = v
        elif k == "resources":
            for pname, resources in v.iteritems():
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
    with open(out_file, "w") as out_handle:
        yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file

def get_handler(args):
    class RunHandler(tornado.web.RequestHandler):
        @tornado.web.asynchronous
        @tornado.gen.coroutine
        def get(self):
            rargs = yaml.safe_load(StringIO.StringIO(str(self.get_argument("args", "{}"))))
            system_config = args.config or "bcbio_system.yaml"
            if "system_config" in rargs:
                system_config = _merge_system_configs(rargs["system_config"],
                                                      system_config, rargs["work_dir"])
            if "sample_config" in rargs:
                sample_config = os.path.join(rargs["work_dir"], "web-sample_config.yaml")
                with open(sample_config, "w") as out_handle:
                    yaml.dump(rargs["sample_config"], out_handle, default_flow_style=False, allow_unicode=False)
            else:
                sample_config = rargs.get("run_config")
            kwargs = {"work_dir": rargs["work_dir"],
                      "config_file": system_config,
                      "run_info_yaml": sample_config,
                      "fc_dir": rargs.get("fc_dir"),
                      "numcores": int(rargs.get("numcores", 1)),
                      "scheduler": rargs.get("scheduler"),
                      "queue": rargs.get("queue"),
                      "resources": rargs.get("resources", ""),
                      "timeout": int(rargs.get("timeout", 15)),
                      "retries": rargs.get("retries"),
                      "app": self.application}
            run_id = yield tornado.gen.Task(run_bcbio_nextgen, **kwargs)
            self.write(run_id)
            self.finish()
    return RunHandler

class StatusHandler(tornado.web.RequestHandler):
    def get(self):
        run_id = self.get_argument("run_id", None)
        if run_id is None:
            status = "not-running"
        else:
            status = self.application.runmonitor.get_status(run_id)
        self.write(status)
        self.finish()
