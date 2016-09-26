#!/bin/env python
"""Exploratory code to convert bcbio generated CWL into WDL.

Uses cwltool parser to build objects, calling out to current cwl2wdl
implementation:

https://github.com/adamstruck/cwl2wdl
"""
import pprint
import sys

from cwl2wdl import generators
from cwl2wdl import base_classes as cwl2wdl_classes
import cwltool.load_tool
import cwltool.workflow
import wdl.parser

def main(wf_file, json_file):
    main_wf = cwltool.load_tool.load_tool(wf_file, cwltool.workflow.defaultMakeTool)
    main_wf_dict = _wf_to_dict(main_wf)
    wdl_doc = generators.WdlWorkflowGenerator(main_wf_dict).generate_wdl()
    print(wdl.parser.parse(wdl_doc))

def _wf_to_dict(wf):
    """Parse a workflow into cwl2wdl style dictionary.
    """
    print wf.tool.keys()
    out = {"name": wf, "inputs": [], "outputs": [], "steps": [], "requirements": []}
    for step in wf.steps:
        if isinstance(step.embedded_tool, cwltool.workflow.Workflow):
            task = _wf_to_dict(step.embedded_tool)
        else:
            task = _tool_to_dict(step.embedded_tool)
    return cwl2wdl_classes.Workflow(out)

def _tool_to_dict(tool):
    """Parse a tool definition into a cwl2wdl style dictionary.
    """
    print(tool)
    print(dir(tool))
    pprint.pprint(tool.tool)
    print(tool.requirements)
    print(tool.hints)
    out = {"name": tool.split("#")[-1], "baseCommand": "", "arguments": [], "inputs": [],
           "outputs": [],
           "requirements": [], "stdin": None, "stdout": None}
    print(out)
    return cwl2wdl_classes.Task(out)

if __name__ == "__main__":
    main(*sys.argv[1:])
