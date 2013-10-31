"""Provide ability to run bcbio-nextgen workflows.
"""
import os
from functools import wraps
from threading import Thread

import tornado.gen
import tornado.web

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
    callback = kwargs["callback"]
    del kwargs["callback"]
    print kwargs
    callback(kwargs["work_dir"])
    run_main(**kwargs)

def get_handler(args):
    class RunHandler(tornado.web.RequestHandler):
        @tornado.web.asynchronous
        @tornado.gen.coroutine
        def get(self):
            kwargs = {"work_dir": str(os.path.abspath(self.get_argument("work_dir"))),
                      "config_file": args.config or "bcbio_system.yaml",
                      "fc_dir": self.get_argument("fc_dir", None),
                      "run_info_yaml": self.get_argument("run_config", None),
                      "numcores": self.get_argument("numcores", None),
                      "scheduler": self.get_argument("scheduler", None),
                      "queue": self.get_argument("queue", None),
                      "resources": self.get_argument("resources", ""),
                      "timeout": int(self.get_argument("timeout", 15)),
                      "retries": self.get_argument("retrieve", None)}
            response = yield tornado.gen.Task(run_bcbio_nextgen, **kwargs)
            self.write(response)
            self.finish()
    return RunHandler
