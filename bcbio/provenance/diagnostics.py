"""Provide logging and diagnostics of running pipelines.

This wraps the BioLite diagnostics database format to provide
tracking of command lines, run times and inspection into run progress.
The goal is to allow traceability and reproducibility of pipelines.

https://bitbucket.org/caseywdunn/biolite
"""

def start_cmd(cmd, descr, entity):
    """Retain details about starting a command, returning a command identifier.
    """
    pass

def end_cmd(cmd_id, succeeded=True):
    """Mark a command as finished with success or failure.
    """
    pass

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
            args = list(args)
            args[item_i] = item
        out.append(args)
    # TODO: store mapping of entity to sub identifiers
    return out

def _has_provenance(x):
    return isinstance(x, dict) and x.has_key("provenance")

def get_item_from_args(xs):
    """Retrieve processed item from list of input arguments.
    """
    for i, x in enumerate(xs):
        if _has_provenance(x):
            return i, x
    return -1, None
