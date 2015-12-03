"""Utility functionality for logging.
"""
import multiprocessing
import os
import socket
import sys
import time

import logbook
import logbook.queues

from bcbio import utils
from bcbio.log import logbook_zmqpush

LOG_NAME = "bcbio-nextgen"

def get_log_dir(config):
    d = config.get("log_dir", "log")
    return d

logger = logbook.Logger(LOG_NAME)
logger_cl = logbook.Logger(LOG_NAME + "-commands")
logger_stdout = logbook.Logger(LOG_NAME + "-stdout")
mpq = multiprocessing.Queue(-1)

def _is_cl(record, _):
    return record.channel == LOG_NAME + "-commands"

def _is_stdout(record, _):
    return record.channel == LOG_NAME + "-stdout"

def _not_cl(record, handler):
    return not _is_cl(record, handler) and not _is_stdout(record, handler)

class CloseableNestedSetup(logbook.NestedSetup):
    def close(self):
        for obj in self.objects:
            if hasattr(obj, "close"):
                obj.close()

class IOSafeMultiProcessingSubscriber(logbook.queues.MultiProcessingSubscriber):
    """Recover from zeromq interrupted system call IOErrors that stop logging.
    """
    def recv(self, timeout=None):
        try:
            return super(IOSafeMultiProcessingSubscriber, self).recv(timeout)
        except IOError, e:
            if "Interrupted system call" in str(e):
                return None
            else:
                raise

def _create_log_handler(config, add_hostname=False, direct_hostname=False):
    logbook.set_datetime_format("utc")
    handlers = [logbook.NullHandler()]
    format_str = "".join(["[{record.time:%Y-%m-%dT%H:%MZ}] " if config.get("include_time", True) else "",
                          "{record.extra[source]}: " if add_hostname else "",
                          "%s: " % (socket.gethostname)() if direct_hostname else "",
                          "{record.message}"])

    log_dir = get_log_dir(config)
    if log_dir:
        if not os.path.exists(log_dir):
            utils.safe_makedir(log_dir)
            # Wait to propagate, Otherwise see logging errors on distributed filesystems.
            time.sleep(5)
        handlers.append(logbook.FileHandler(os.path.join(log_dir, "%s.log" % LOG_NAME),
                                            format_string=format_str, level="INFO",
                                            filter=_not_cl))
        handlers.append(logbook.FileHandler(os.path.join(log_dir, "%s-debug.log" % LOG_NAME),
                                            format_string=format_str, level="DEBUG", bubble=True,
                                            filter=_not_cl))
        handlers.append(logbook.FileHandler(os.path.join(log_dir, "%s-commands.log" % LOG_NAME),
                                            format_string=format_str, level="DEBUG",
                                            filter=_is_cl))
    handlers.append(logbook.StreamHandler(sys.stdout, format_string="{record.message}",
                                          level="DEBUG", filter=_is_stdout))

    email = config.get("email", config.get("resources", {}).get("log", {}).get("email"))
    if email:
        email_str = u'''Subject: [bcbio-nextgen] {record.extra[run]} \n\n {record.message}'''
        handlers.append(logbook.MailHandler(email, [email],
                                            format_string=email_str,
                                            level='INFO', bubble=True))

    handlers.append(logbook.StreamHandler(sys.stderr, format_string=format_str, bubble=True,
                                          filter=_not_cl))
    return CloseableNestedSetup(handlers)

def create_base_logger(config=None, parallel=None):
    """Setup base logging configuration, also handling remote logging.

    Correctly sets up for local, multiprocessing and distributed runs.
    Creates subscribers for non-local runs that will be references from
    local logging.
    """
    if parallel is None: parallel = {}
    parallel_type = parallel.get("type", "local")
    cores = parallel.get("cores", 1)
    if parallel_type == "ipython":
        ips = [ip for ip in socket.gethostbyname_ex(socket.gethostname())[2]
               if not ip.startswith("127.0.0")]
        if not ips:
            sys.stderr.write("Cannot resolve a local IP address that isn't 127.0.0. "
                             "Your machines might not have a local IP address "
                             "assigned or are not able to resolve it.\n")
            sys.exit(1)
        uri = "tcp://%s" % ips[0]
        subscriber = logbook_zmqpush.ZeroMQPullSubscriber()
        mport = subscriber.socket.bind_to_random_port(uri)
        wport_uri = "%s:%s" % (uri, mport)
        parallel["log_queue"] = wport_uri
        subscriber.dispatch_in_background(_create_log_handler(config, True))
    elif cores > 1:
        subscriber = IOSafeMultiProcessingSubscriber(mpq)
        subscriber.dispatch_in_background(_create_log_handler(config))
    else:
        # Do not need to setup anything for local logging
        pass
    return parallel

def setup_local_logging(config=None, parallel=None):
    """Setup logging for a local context, directing messages to appropriate base loggers.

    Handles local, multiprocessing and distributed setup, connecting
    to handlers created by the base logger.
    """
    if config is None: config = {}
    if parallel is None: parallel = {}
    parallel_type = parallel.get("type", "local")
    cores = parallel.get("cores", 1)
    wrapper = parallel.get("wrapper", None)
    if parallel_type == "ipython":
        handler = logbook_zmqpush.ZeroMQPushHandler(parallel["log_queue"])
    elif cores > 1:
        handler = logbook.queues.MultiProcessingHandler(mpq)
    else:
        handler = _create_log_handler(config, direct_hostname=wrapper is not None)
    handler.push_thread()
    return handler
