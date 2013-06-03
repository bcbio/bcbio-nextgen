"""Create bcbio_sample.yaml files from standard templates and lists of input files.

Provides an automated way to generate a full set of analysis files from an input
YAML template. Default templates are provided for common approaches which can be tweaked
as needed.
"""
import contextlib
import copy
import datetime
import itertools
import os
import urllib2

import yaml

from bcbio import utils
from bcbio.bam import fastq
from bcbio.variation.cortex import get_sample_name
from bcbio.workflow.xprize import HelpArgParser

def parse_args(inputs):
    parser = HelpArgParser(
        description="Create a bcbio_sample.yaml file from a standard template and inputs")
    parser.add_argument("template", help="Template name or path to template YAML file")
    parser.add_argument("project_name", help="Name of the current project")
    parser.add_argument("input_files", nargs="+", help="Input read files, in BAM or fastq format")
    return parser.parse_args(inputs)

# ## Prepare sequence data inputs

def _prep_bam_input(f, base):
    if not os.path.exists(f):
        raise ValueError("Could not find input file: %s" % f)
    cur = copy.deepcopy(base)
    cur["files"] = os.path.abspath(f)
    cur["description"] = get_sample_name(f)
    cur["algorithm"]["aligner"] = False
    return cur

def _prep_fastq_input(fs, base):
    for f in fs:
        if not os.path.exists(f):
            raise ValueError("Could not find input file: %s" % f)
    cur = copy.deepcopy(base)
    cur["files"] = [os.path.abspath(f) for f in fs]
    d = os.path.commonprefix([os.path.splitext(os.path.basename(f))[0] for f in fs])
    to_strip = ("_R", "_", "fastq", ".")
    while d.endswith(to_strip):
        for x in to_strip:
            if d.endswith(x):
                d = d[:len(d) - len(x)]
                break
    cur["description"] = d
    return cur

def _prep_items_from_base(base, in_files):
    """Prepare a set of configuration items for input files.
    """
    details = []
    fq_exts = [".fq", ".fastq", ".txt", ".gz"]
    gz_exts = tuple(["%s.gz" % ext for ext in fq_exts if ext != ".gz"])
    for ext, files in itertools.groupby(in_files, lambda x: os.path.splitext(x)[-1].lower()):
        if ext == ".bam":
            for f in files:
                details.append(_prep_bam_input(f, base))
        elif ext in fq_exts:
            files = list(files)
            if ext == ".gz": assert all(f.endswith(gz_exts) for f in files), (files, gz_exts)
            for fs in fastq.combine_pairs(files):
                details.append(_prep_fastq_input(fs, base))
        else:
            raise ValueError("File type not yet implemented: %s" % ext)
    return details

# ## Read and write configuration files

def _read_template(template):
    """Read template file into a dictionary to use as base for all samples.

    Handles well-known template names, pulled from GitHub repository and local
    files.
    """
    if os.path.isfile(template):
        with open(template) as in_handle:
            return yaml.load(in_handle)
    else:
        base_url = "https://raw.github.com/chapmanb/bcbio-nextgen/master/config/templates/%s.yaml"
        try:
            with contextlib.closing(urllib2.urlopen(base_url % template)) as in_handle:
                return yaml.load(in_handle)
        except urllib2.HTTPError:
            raise ValueError("Could not find template '%s' locally or in standard templates on GitHub"
                             % template)


def _write_config_file(items, template, project_name, out_dir):
    """Write configuration file, adding required top level attributes.
    """
    config_dir = utils.safe_makedir(os.path.join(out_dir, "config"))
    out_config_file = os.path.join(config_dir, "%s.yaml" % project_name)
    out = {"fc_date": datetime.datetime.now().strftime("%y%m%d"),
           "fc_name": project_name,
           "upload" : {"dir": "../final"},
           "details": items}
    for k, v in template.iteritems():
        if k not in ["details"]:
            out[k] = v
    with open(out_config_file, "w") as out_handle:
        yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_config_file

def setup(args):
    template = _read_template(args.template)
    base_item = template["details"][0]
    items = _prep_items_from_base(base_item, args.input_files)

    project_name = args.project_name.replace(" ", "_")
    out_dir = os.path.join(os.getcwd(), project_name)
    out_config_file = _write_config_file(items, template, project_name, out_dir)
    work_dir = utils.safe_makedir(os.path.join(out_dir, "work"))
    print "Configuration file created at: %s" % out_config_file
    print "Edit to finalize and run with:"
    print "  cd %s" % work_dir
    print "  bcbio_nextgen.py /path/to/bcbio_system.yaml %s" % out_config_file
