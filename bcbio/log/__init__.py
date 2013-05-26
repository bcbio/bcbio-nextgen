"""Utility functionality for logging.
"""
import multiprocessing
import os
import socket
import sys

import logbook, logbook.queues

from bcbio import utils

LOG_NAME = "bcbio-nextgen"

def _get_log_dir(config):
    d = config.get("log_dir",
                   config.get("resources", {}).get("log", {}).get("dir", "log"))
    return d

logger = logbook.Logger(LOG_NAME)
mpq = multiprocessing.Queue(-1)

def _add_host(record):
    record.extra["hostname"] = socket.gethostname()
    return record

def _create_log_handler(config):
    handlers = [logbook.Processor(_add_host),
                logbook.NullHandler()]
    format_str = ("[{record.time:%Y-%m-%d %H:%M}] "
                  "{record.extra[hostname]}: {record.message}")

    log_dir = _get_log_dir(config)
    if log_dir:
        utils.safe_makedir(log_dir)
        handlers.append(logbook.FileHandler(os.path.join(log_dir, "%s.log" % LOG_NAME),
                                            format_string=format_str, level="INFO"))
        handlers.append(logbook.FileHandler(os.path.join(log_dir, "%s-debug.log" % LOG_NAME),
                                            format_string=format_str, level="DEBUG", bubble=True))

    email = config.get("email", config.get("resources", {}).get("log", {}).get("email"))
    if email:
        email_str = u'''Subject: [bcbio-nextgen] {record.extra[run]} \n\n {record.message}'''
        handlers.append(logbook.MailHandler(email, [email],
                                            format_string=email_str,
                                            level='INFO', bubble = True))

    handlers.append(logbook.StreamHandler(sys.stderr, format_string=format_str, bubble=True))
    return logbook.NestedSetup(handlers)

def create_base_logger(config, parallel=None):
    """Setup base logging configuration, also handling remote logging.

    Correctly sets up for local, multiprocessing and distributed runs.
    Creates subscribers for non-local runs that will be references from
    local logging.
    """
    if parallel is None: parallel = {}
    parallel_type = parallel.get("type", "local")
    cores = parallel.get("cores", 1)
    if parallel_type == "ipython":
        parallel["log_queue"] = "tcp://%s:5000" % socket.gethostbyname(socket.gethostname())
        subscriber = logbook.queues.ZeroMQSubscriber(parallel["log_queue"])
        subscriber.dispatch_in_background(_create_log_handler(config))
    elif cores > 1:
        subscriber = logbook.queues.MultiProcessingSubscriber(mpq)
        subscriber.dispatch_in_background(_create_log_handler(config))
    else:
        # Do not need to setup anything for local logging
        pass
    return parallel

def setup_local_logging(config, parallel=None):
    """Setup logging for a local context, directing messages to appropriate base loggers.

    Handles local, multiprocessing and distributed setup, connecting
    to handlers created by the base logger.
    """
    if parallel is None: parallel = {}
    parallel_type = parallel.get("type", "local")
    cores = parallel.get("cores", 1)
    if parallel_type == "ipython":
        handler = logbook.queues.ZeroMQHandler(parallel["log_queue"])
    elif cores > 1:
        handler = logbook.queues.MultiProcessingHandler(mpq)
    else:
        handler = _create_log_handler(config)
    handler.push_thread()
