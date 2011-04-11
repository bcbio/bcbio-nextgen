"""Utility functionality for logging.
"""
import os
import sys

import logbook

def create_log_handler(config, log_name):
    log_dir = config.get("log_dir", None)
    email = config.get("email_notify", None)
    
    if log_dir:
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        handler = logbook.FileHandler(os.path.join(log_dir, "%s.log" % log_name))
    else:
        handler = logbook.StreamHandler(sys.stdout)
        
    if email:
        handler = logbook.MailHandler("ngsmaster@example.com", [email],
                                      format_string=logbook.handlers.MAIL_FORMAT_STRING,
                                      level='INFO', bubble = True)
    return handler