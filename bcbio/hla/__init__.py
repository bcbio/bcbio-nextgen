"""Main entry point for HLA typing on samples.
"""
import toolz as tz

from bcbio import utils
from bcbio.hla import bwakit, optitype

_CALLERS = {"bwakit": bwakit.run,
            "optitype": optitype.run}

def call_hla(data):
    data = utils.to_single_data(data)
    hlacaller = tz.get_in(["config", "algorithm", "hlacaller"], data)
    if hlacaller and "hla" in data and tz.get_in(["hla", "fastq"], data):
        data = _CALLERS[hlacaller](data)
    elif "hla" not in data:
        data["hla"] = {}
    return [[data]]

def run(samples, run_parallel):
    """Run HLA detection on the input samples.
    """
    to_process = []
    extras = []
    for data in (xs[0] for xs in samples):
        hlacaller = tz.get_in(["config", "algorithm", "hlacaller"], data)
        if hlacaller:
            to_process.append(data)
        else:
            extras.append([data])
    processed = run_parallel("call_hla", ([x] for x in to_process))
    return extras + processed
