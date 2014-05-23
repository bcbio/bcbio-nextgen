"""Parsing of command line arguments into parallel inputs.
"""

def to_parallel(args, module="bcbio.distributed"):
    """Convert input arguments into a parallel dictionary for passing to processing.
    """
    ptype, cores = _get_cores_and_type(args.numcores, getattr(args, "paralleltype", None),
                                       args.scheduler)
    parallel = {"type": ptype, "cores": cores,
                "scheduler": args.scheduler, "queue": args.queue,
                "tag": args.tag, "module": module,
                "resources": args.resources, "timeout": args.timeout,
                "retries": args.retries,
                "run_local": args.queue == "localrun"}
    return parallel

def _get_cores_and_type(numcores, paralleltype, scheduler):
    """Return core and parallelization approach from command line providing sane defaults.
    """
    if scheduler is not None:
        paralleltype = "ipython"
    if paralleltype is None:
        paralleltype = "local"
    if not numcores or int(numcores) < 1:
        numcores = 1
    return paralleltype, int(numcores)
