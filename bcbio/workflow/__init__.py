"""Prepare pre-defined workflows that handle creation of input sample configs.
"""
from bcbio.workflow import xprize, stormseq

workflows = {"xprize": xprize,
             "stormseq": stormseq}

def setup(name, inputs):
    workflow = workflows[name]
    args = workflow.parse_args(inputs)
    return workflow.setup(args)
