"""Provide logging and diagnostics of running pipelines.

This wraps the BioLite diagnostics database format to provide
tracking of command lines, run times and inspection into run progress.
The goal is to allow traceability and reproducibility of pipelines.

https://bitbucket.org/caseywdunn/biolite

Also interfaces with Galaxy's history export format, to provide
a set of metadata that could be imported into object stores.
"""
import os
import uuid

import toolz as tz
try:
    import biolite, biolite.config, biolite.database
except ImportError:
    biolite = None

from bcbio import utils

def start_cmd(cmd, descr, data):
    """Retain details about starting a command, returning a command identifier.
    """
    if data and "provenance" in data:
        entity_id = tz.get_in(["provenance", "entity"], data)

def end_cmd(cmd_id, succeeded=True):
    """Mark a command as finished with success or failure.
    """
    pass

def initialize(dirs):
    """Initialize the biolite database to load provenance information.
    """
    if biolite and dirs.get("work"):
        base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
        p_db = os.path.join(base_dir, "biolite.db")
        biolite.config.resources["database"] = p_db
        biolite.database.connect()

def store_entity(data):
    fc_name = tz.get_in(["upload", "fn_name"], data)
    run_id = ("%s.%s" % (tz.get_in(["upload", "fc_date"], data), fc_name)
              if fc_name else str(uuid.uuid1()))
    return str(uuid.uuid1())

def track_parallel(items, sub_type):
    """Create entity identifiers to trace the given items in sub-commands.

    Helps handle nesting in parallel program execution:

    run id => sub-section id => parallel ids
    """
    out = []
    for i, args in enumerate(items):
        item_i, item = _get_provitem_from_args(args)
        if item:
            sub_entity = "%s.%s.%s" % (item["provenance"]["entity"], sub_type, i)
            item["provenance"]["entity"] = sub_entity
            args = list(args)
            args[item_i] = item
        out.append(args)
    # TODO: store mapping of entity to sub identifiers
    return out

def _has_provenance(x):
    return isinstance(x, dict) and "provenance" in x

def _get_provitem_from_args(xs):
    """Retrieve processed item from list of input arguments.
    """
    for i, x in enumerate(xs):
        if _has_provenance(x):
            return i, x
    return -1, None
