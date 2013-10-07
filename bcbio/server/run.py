"""Provide ability to run bcbio-nextgen workflows.
"""
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
def run_bcbio_nextgen(system_config, fc_dir, run_yaml, callback):
    print system_config, fc_dir, run_yaml
    callback("run_id")

def get_handler(args):
    class RunHandler(tornado.web.RequestHandler):
        @tornado.web.asynchronous
        @tornado.gen.coroutine
        def get(self):
            run_yaml = self.get_argument("run_yaml", None)
            fc_dir = self.get_argument("fc_dir", None)
            system_config = args.config
            response = yield tornado.gen.Task(run_bcbio_nextgen, system_config, fc_dir, run_yaml)
            self.write(response)
            self.finish()

    return RunHandler
