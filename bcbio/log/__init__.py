"""Utility functionality for logging.
"""
import multiprocessing
import os
import socket
import sys

import logging

from bcbio import utils

LOG_NAME = "bcbio-nextgen"

logger = logging.getLogger(LOG_NAME)

def _get_log_dir(config):
    d = config.get("log_dir",
                   config.get("resources", {}).get("log", {}).get("dir", "log"))
    return d

def setup_logging(config):
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        formatter = logging.Formatter('[%(asctime)s] %(message)s')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        log_dir = _get_log_dir(config)
        if log_dir:
            logfile = os.path.join(utils.safe_makedir(log_dir),
                                   "{0}.log".format(LOG_NAME))
            handler = logging.FileHandler(logfile)
            handler.setFormatter(formatter)
            logger.addHandler(handler)

import logbook
logger2 = logbook.Logger(LOG_NAME)
mpq = multiprocessing.Queue(-1)

def _add_host(record):
    record.hostname = socket.gethostname()
    return record

def _create_log_handler(config):
    handlers = [logbook.Processor(_add_host),
                logbook.NullHandler()]

    log_dir = _get_log_dir(config)
    if log_dir:
        utils.safe_makedir(log_dir)
        handlers.append(logbook.FileHandler(os.path.join(log_dir, "%s.log" % LOG_NAME)),
                                           level="INFO")
        handlers.append(logbook.FileHandler(os.path.join(log_dir, "%s-debug.log" % LOG_NAME)),
                                           level="DEBUG", bubble=True)

    email = config.get("email", config.get("resources", {}).get("log", {}).get("email"))
    if email:
        handlers.append(logbook.MailHandler(email, [email],
                                            format_string=u'''Subject: [bcbio-nextgen] {record.extra[run]} \n\n {record.message}''',
                                            level='INFO', bubble = True))

    handlers.append(logbook.StreamHandler(sys.stderr, bubble=True))
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
        parallel["log_queue"] = "tcp://%s:5000" % socket.gethostname()
        subscriber = logbook.queues.ZeroMQSubscriber(parallel(["log_queue"]))
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
