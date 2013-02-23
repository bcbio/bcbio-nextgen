"""Prepare pre-defined workflows that handle creation of input sample configs.
"""
from bcbio.workflow import xprize

workflows = {"xprize": xprize}

def setup(name, inputs):
    workflow = workflows[name]
    args = workflow.parse_args(inputs)
    return workflow.setup(args)
