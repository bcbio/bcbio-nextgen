"""Provide logging and diagnostics of running pipelines.

This wraps the BioLite diagnostics database format to provide
tracking of command lines, run times and inspection into run progress.
The goal is to allow traceability and reproducibility of pipelines.

https://bitbucket.org/caseywdunn/biolite
"""

from bcbio.log import logger

def log_cmd(descr, entity, cmd):
    """Log and save details when running a program.
    """
    logger.info(descr)
    # TODO: storage of command line mapped to entity

def track_parallel(items, sub_type):
    """Create entity identifiers to trace the given items in sub-commands.

    Helps handle nesting in parallel program execution:

    run id => sub-section id => parallel ids
    """
    out = []
    for i, args in enumerate(items):
        item_i, item = get_item_from_args(args)
        if item:
            sub_entity = "%s.%s.%s" % (item["provenance"]["entity"], sub_type, i)
            item["provenance"]["entity"] = sub_entity
            args[item_i] = item
        out.append(args)
    # TODO: store mapping of entity to sub identifiers
    return out

def get_item_from_args(xs):
    """Retrieve processed item from list of input arguments.
    """
    for i, x in enumerate(xs):
        if isinstance(x, dict) and x.has_key("provenance"):
            return i, x
    return -1, None
