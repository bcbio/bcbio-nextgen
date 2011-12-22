"""Utility functionality for logging.
"""
import os
import sys

import logbook

from bcbio import utils

LOG_NAME = "nextgen_pipeline"

logger = logbook.Logger(LOG_NAME)

def create_log_handler(config):
    log_dir = config.get("log_dir", None)
    email = config.get("email", None)
    
    if log_dir:
        utils.safe_makedir(log_dir)
        handler = logbook.FileHandler(os.path.join(log_dir, "%s.log" % LOG_NAME))
    else:
        handler = logbook.StreamHandler(sys.stdout)
        
    if email:
        handler = logbook.MailHandler(email, [email], 
                                      format_string=u'''Subject: [BCBB pipeline] {record.extra[run]} \n\n {record.message}''',
                                      level='INFO', bubble = True)
    return handler
