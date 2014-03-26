"""Top level functionality for running a bcbio-nextgen web server allowing remote jobs.
"""
import tornado.web
import tornado.ioloop

from bcbio.server import run

def start(args):
    """Run server with provided command line arguments.
    """
    application = tornado.web.Application([(r"/run", run.get_handler(args)),
                                           (r"/status", run.StatusHandler)])
    application.runmonitor = RunMonitor()
    application.listen(args.port)
    tornado.ioloop.IOLoop.instance().start()

class RunMonitor:
    """Track current runs and provide status.
    """
    def __init__(self):
        self._running = {}

    def set_status(self, run_id, status):
        self._running[run_id] = status

    def get_status(self, run_id):
        return self._running.get(run_id, "not-running")

def add_subparser(subparsers):
    """Add command line arguments as server subparser.
    """
    parser = subparsers.add_parser("server", help="Run a bcbio-nextgen server allowing remote job execution.")
    parser.add_argument("-c", "--config", help=("Global YAML configuration file specifying system details."
                                                "Defaults to installed bcbio_system.yaml"))
    parser.add_argument("-p", "--port", help="Port to listen on (default 8080)",
                        default=8080, type=int)
    parser.add_argument("-n", "--cores", help="Cores to use when processing locally when not requested (default 1)",
                        default=1, type=int)
    parser.add_argument("-d", "--biodata_dir", help="Directory with biological data",
                        default="/mnt/biodata", type=str)
    return parser
