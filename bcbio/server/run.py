"""Provide ability to run bcbio-nextgen workflows.
"""
import collections
import os
import StringIO
import sys
import uuid

import tornado.gen
import tornado.web
import yaml

from bcbio import utils
from bcbio.distributed import clargs
from bcbio.server import background

def run_bcbio_nextgen(**kwargs):
    callback = kwargs.pop("callback", None)
    app = kwargs.pop("app")
    args = [x for x in [kwargs["config_file"], kwargs["fc_dir"], kwargs["run_info_yaml"]] if x]
    run_id = str(uuid.uuid1())
    def set_done(status, stdout, stderr, has_timed_out):
        app.runmonitor.set_status(run_id, "finished" if status == 0 else "failed")
    if utils.get_in(kwargs, ("parallel", "type")) == "local":
        _run_local(kwargs["workdir"], args, utils.get_in(kwargs, ("parallel", "cores")), set_done)
    else:
        # XXX Need to work on ways to prepare batch scripts for bcbio submission
        # when analysis server talks to an HPC cluster
        raise ValueError("Do not yet support automated execution of this parallel config: %s" % parallel)
    app.runmonitor.set_status(run_id, "running")
    if callback:
        callback(run_id)
    else:
        return run_id

def _run_local(workdir, args, cores, callback):
    cmd = [os.path.join(os.path.dirname(sys.executable), "bcbio_nextgen.py")] + args + \
          ["-n", cores]
    with utils.chdir(workdir):
        p = background.Subprocess(callback, timeout=-1, args=[str(x) for x in cmd])
        p.start()

def _rargs_to_parallel_args(rargs, args):
    Args = collections.namedtuple("Args", "numcores scheduler queue resources timeout retries tag")
    return Args(int(rargs.get("numcores", args.cores)), rargs.get("scheduler"),
                rargs.get("queue"), rargs.get("resources", ""),
                int(rargs.get("timeout", 15)), rargs.get("retries"),
                rargs.get("tag"))

def get_handler(args):
    class RunHandler(tornado.web.RequestHandler):
        @tornado.web.asynchronous
        @tornado.gen.coroutine
        def get(self):
            rargs = yaml.safe_load(StringIO.StringIO(str(self.get_argument("args", "{}"))))
            system_config = args.config or "bcbio_system.yaml"
            if "system_config" in rargs:
                system_config = os.path.join(rargs["work_dir"], "web-system_config.yaml")
                with open(system_config, "w") as out_handle:
                    yaml.safe_dump(rargs["system_config"], out_handle, default_flow_style=False, allow_unicode=False)
            if "sample_config" in rargs:
                sample_config = os.path.join(rargs["work_dir"], "web-sample_config.yaml")
                with open(sample_config, "w") as out_handle:
                    yaml.safe_dump(rargs["sample_config"], out_handle, default_flow_style=False, allow_unicode=False)
            else:
                sample_config = rargs.get("run_config")
            kwargs = {"workdir": rargs["work_dir"],
                      "config_file": system_config,
                      "run_info_yaml": sample_config,
                      "fc_dir": rargs.get("fc_dir"),
                      "parallel": clargs.to_parallel(_rargs_to_parallel_args(rargs, args)),
                      "app": self.application}
            run_id = yield tornado.gen.Task(run_bcbio_nextgen, **kwargs)
            self.write(run_id)
            self.finish()
    return RunHandler

class StatusHandler(tornado.web.RequestHandler):
    def get(self):
        run_id = self.get_argument("run_id", None)
        if run_id is None:
            status = "server-up"
        else:
            status = self.application.runmonitor.get_status(run_id)
        self.write(status)
        self.finish()
