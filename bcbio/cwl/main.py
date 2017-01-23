"""Main entry point for running preparation of CWL inputs.
"""
from bcbio.pipeline import run_info
from bcbio.cwl import create

def run(args):
    """Run a CWL preparation pipeline.
    """
    dirs, config, run_info_yaml = run_info.prep_system(args.sample_config, args.systemconfig)
    integrations = args.integrations if hasattr(args, "integrations") else {}
    world = run_info.organize(dirs, config, run_info_yaml, add_provenance=False, integrations=integrations)
    create.from_world(world, run_info_yaml, integrations=integrations)
