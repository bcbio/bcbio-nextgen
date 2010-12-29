"""Utility functionality for logging.
"""
import os

import logbook

def create_log_handler(config, log_name):
    log_dir = config["log_dir"]
    if log_dir:
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        handler = logbook.FileHandler(os.path.join(log_dir, "%s.log" % log_name))
    else:
        handler = logbook.StreamHandler()
    return handler
