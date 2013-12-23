"""Create bcbio_sample.yaml files from standard templates and lists of input files.

Provides an automated way to generate a full set of analysis files from an input
YAML template. Default templates are provided for common approaches which can be tweaked
as needed.
"""
import contextlib
import copy
import csv
import datetime
import glob
import itertools
import os
import shutil
import urllib2

import yaml

from bcbio import utils
from bcbio.bam import fastq, sample_name
from bcbio.pipeline import run_info
from bcbio.workflow.xprize import HelpArgParser

def parse_args(inputs):
    parser = HelpArgParser(
        description="Create a bcbio_sample.yaml file from a standard template and inputs")
    parser.add_argument("template", help=("Template name or path to template YAML file. "
                                          "Built in choices: freebayes-variant, gatk-variant, tumor-paired, "
                                          "noalign-variant, illumina-rnaseq, illumina-chipseq"))
    parser.add_argument("metadata", help="CSV file with project metadata. Name of file used as project name.")
    parser.add_argument("input_files", nargs="*", help="Input read files, in BAM or fastq format")
    return parser.parse_args(inputs)

# ## Prepare sequence data inputs

def _prep_bam_input(f, i, base):
    if not os.path.exists(f):
        raise ValueError("Could not find input file: %s" % f)
    cur = copy.deepcopy(base)
    cur["files"] = [os.path.abspath(f)]
    cur["description"] = sample_name(f) or os.path.splitext(os.path.basename(f))[0]
    return cur

def _prep_fastq_input(fs, base):
    for f in fs:
        if not os.path.exists(f):
            raise ValueError("Could not find input file: %s" % f)
    cur = copy.deepcopy(base)
    cur["files"] = [os.path.abspath(f) for f in fs]
    d = os.path.commonprefix([utils.splitext_plus(os.path.basename(f))[0] for f in fs])
    cur["description"] = fastq.rstrip_extra(d)
    return cur

def _prep_items_from_base(base, in_files):
    """Prepare a set of configuration items for input files.
    """
    details = []
    known_exts = {".bam": "bam", ".fq": "fastq",
                  ".fastq": "fastq", ".txt": "fastq",
                  ".fastq.gz": "fastq", ".fq.gz": "fastq",
                  ".txt.gz": "fastq", ".gz": "fastq"}
    in_files = _expand_dirs(in_files, known_exts)
    in_files = _expand_wildcards(in_files)

    for i, (ext, files) in enumerate(itertools.groupby(
            in_files, lambda x: known_exts.get(utils.splitext_plus(x)[-1].lower()))):
        if ext == "bam":
            for f in files:
                details.append(_prep_bam_input(f, i, base))
        elif ext == "fastq":
            files = list(files)
            for fs in fastq.combine_pairs(files):
                details.append(_prep_fastq_input(fs, base))
        else:
            raise ValueError("Unexpected input file types: %s" % str(files))
    return details

def _expand_dirs(in_files, known_exts):
    def _is_dir(in_file):
        return os.path.isdir(os.path.expanduser(in_file))
    files, dirs = utils.partition(_is_dir, in_files)
    for dir in dirs:
        for ext in known_exts.keys():
            wildcard = os.path.join(os.path.expanduser(dir), "*"+ext)
            files = itertools.chain(glob.glob(wildcard), files)
    return list(files)

def _expand_wildcards(in_files):
    def _has_wildcard(in_file):
        return "*" in in_file

    files, wildcards = utils.partition(_has_wildcard, in_files)
    for wc in wildcards:
        abs_path = os.path.expanduser(wc)
        files = itertools.chain(glob.glob(abs_path), files)
    return list(files)

# ## Read and write configuration files

def _read_template(template):
    """Read template file into a dictionary to use as base for all samples.

    Handles well-known template names, pulled from GitHub repository and local
    files.
    """
    if os.path.isfile(template):
        with open(template) as in_handle:
            txt_config = in_handle.read()
        with open(template) as in_handle:
            config = yaml.load(in_handle)
    else:
        base_url = "https://raw.github.com/chapmanb/bcbio-nextgen/master/config/templates/%s.yaml"
        try:
            with contextlib.closing(urllib2.urlopen(base_url % template)) as in_handle:
                txt_config = in_handle.read()
            with contextlib.closing(urllib2.urlopen(base_url % template)) as in_handle:
                config = yaml.load(in_handle)
        except urllib2.HTTPError:
            raise ValueError("Could not find template '%s' locally or in standard templates on GitHub"
                             % template)
    return config, txt_config

def _write_template_config(template_txt, project_name, out_dir):
    config_dir = utils.safe_makedir(os.path.join(out_dir, "config"))
    out_config_file = os.path.join(config_dir, "%s-template.yaml" % project_name)
    with open(out_config_file, "w") as out_handle:
        out_handle.write(template_txt)
    return out_config_file

def _write_config_file(items, template, project_name, out_dir):
    """Write configuration file, adding required top level attributes.
    """
    config_dir = utils.safe_makedir(os.path.join(out_dir, "config"))
    out_config_file = os.path.join(config_dir, "%s.yaml" % project_name)
    out = {"fc_date": datetime.datetime.now().strftime("%Y-%m-%d"),
           "fc_name": project_name,
           "upload": {"dir": "../final"},
           "details": items}
    for k, v in template.iteritems():
        if k not in ["details"]:
            out[k] = v
    if os.path.exists(out_config_file):
        shutil.move(out_config_file,
                    out_config_file + ".bak%s" % datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))
    with open(out_config_file, "w") as out_handle:
        yaml.dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_config_file

def _safe_name(x):
    for prob in [" ", "."]:
        x = x.replace(prob, "_")
    return x

def _parse_metadata(in_file):
    """Reads metadata from a simple CSV structured input file.

    samplename,batch,phenotype
    ERR256785,batch1,normal
    """
    metadata = {}
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        while 1:
            header = reader.next()
            if not header[0].startswith("#"):
                break
        keys = [x.strip() for x in header[1:]]
        for sinfo in (x for x in reader if not x[0].startswith("#")):
            sample = sinfo[0].strip()
            metadata[sample] = dict(zip(keys, (x.strip() for x in sinfo[1:])))
    return metadata

def _pname_and_metadata(in_file):
    """Retrieve metadata and project name from the input metadata CSV file.

    Uses the input file name for the project name and

    For back compatibility, accepts the project name as an input, providing no metadata.
    """
    if not os.path.isfile(in_file):
        return _safe_name(in_file), {}
    else:
        return (_safe_name(os.path.splitext(os.path.basename(in_file))[0]),
                _parse_metadata(in_file))

def _add_metadata(item, metadata):
    """Add metadata information from CSV file to current item.

    Retrieves metadata based on 'description' parsed from input CSV file.
    Adds to object and handles special keys:
    - `description`: A new description for the item. Used to relabel items
       based on the pre-determined description from fastq name or BAM read groups.
    - Keys matching supported names in the algorithm section map
      to key/value pairs there instead of metadata.
    """
    item_md = metadata.get(item["description"],
                           metadata.get(os.path.basename(item["files"][0]), {}))
    TOP_LEVEL = set(["description"])
    if len(item_md) > 0:
        if "metadata" not in item:
            item["metadata"] = {}
        for k, v in item_md.iteritems():
            if v:
                if k in TOP_LEVEL:
                    item[k] = v
                elif k in run_info.ALGORITHM_KEYS:
                    # Handle keys that pass integer or boolean values
                    try:
                        v = int(v)
                    except ValueError:
                        if v.lower() == "true":
                            v = True
                        elif v.lower() == "false":
                            v = False
                    item["algorithm"][k] = v
                else:
                    item["metadata"][k] = v
    elif len(metadata) > 0:
        print "Metadata not found for sample %s, %s" % (item["description"],
                                                        os.path.basename(item["files"][0]))
    return item

def setup(args):
    template, template_txt = _read_template(args.template)
    base_item = template["details"][0]
    project_name, metadata = _pname_and_metadata(args.metadata)
    items = [_add_metadata(item, metadata)
             for item in _prep_items_from_base(base_item, args.input_files)]

    out_dir = os.path.join(os.getcwd(), project_name)
    work_dir = utils.safe_makedir(os.path.join(out_dir, "work"))
    if len(items) == 0:
        out_config_file = _write_template_config(template_txt, project_name, out_dir)
        print "Template configuration file created at: %s" % out_config_file
        print "Edit to finalize custom options, then prepare full sample config with:"
        print "  bcbio_nextgen.py -w template %s %s sample1.bam sample2.fq" % \
            (out_config_file, project_name)
    else:
        out_config_file = _write_config_file(items, template, project_name, out_dir)
        print "Configuration file created at: %s" % out_config_file
        print "Edit to finalize and run with:"
        print "  cd %s" % work_dir
        print "  bcbio_nextgen.py ../config/%s" % os.path.basename(out_config_file)
