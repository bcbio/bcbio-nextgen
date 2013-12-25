"""Top level functionality for running a bcbio-nextgen web server allowing remote jobs.
"""
import tornado.web
import tornado.ioloop

from bcbio.server import run, install

def start(args):
    """Run server with provided command line arguments.
    """
    application = tornado.web.Application([(r"/run", run.get_handler(args)),
                                           (r"/install", install.get_handler(args))])
    application.listen(args.port)
    tornado.ioloop.IOLoop.instance().start()

def add_subparser(subparsers):
    """Add command line arguments as server subparser.
    """
    parser = subparsers.add_parser("server", help="Run a bcbio-nextgen server allowing remote job execution.")
    parser.add_argument("-c", "--config", help=("Global YAML configuration file specifying system details."
                                                "Defaults to installed bcbio_system.yaml"))
    parser.add_argument("-p", "--port", help="Port to listen on",
                        default=8080, type=int)
    parser.add_argument("-d", "--biodata_dir", help="Directory with biological data",
                        default="/mnt/biodata", type=str)
    return parser
