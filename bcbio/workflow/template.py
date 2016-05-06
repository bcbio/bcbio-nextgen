"""Create bcbio_sample.yaml files from standard templates and lists of input files.

Provides an automated way to generate a full set of analysis files from an inpu
YAML template. Default templates are provided for common approaches which can be tweaked
as needed.
"""
import collections
import contextlib
import copy
import csv
import datetime
import glob
import itertools
import os
import shutil
import urllib2

import toolz as tz
import yaml

from bcbio import utils
from bcbio.bam import fastq, sample_name
from bcbio.distributed import objectstore
from bcbio.upload import s3
from bcbio.pipeline import config_utils, run_info
from bcbio.workflow.xprize import HelpArgParser

def parse_args(inputs):
    parser = HelpArgParser(
        description="Create a bcbio_sample.yaml file from a standard template and inputs")
    parser = setup_args(parser)
    return parser.parse_args(inputs)

def setup_args(parser):
    parser.add_argument("template", help=("Template name or path to template YAML file. "
                                          "Built in choices: freebayes-variant, gatk-variant, tumor-paired, "
                                          "noalign-variant, illumina-rnaseq, illumina-chipseq"))
    parser.add_argument("metadata", help="CSV file with project metadata. Name of file used as project name.")
    parser.add_argument("input_files", nargs="*", help="Input read files, in BAM or fastq format")
    parser.add_argument("--only-metadata", help="Ignore samples not present in metadata CSV file",
                        action="store_true", default=False)
    return parser

# ## Prepare sequence data inputs

def _prep_bam_input(f, i, base):
    if not os.path.exists(f) and not objectstore.is_remote(f):
        raise ValueError("Could not find input file: %s" % f)
    cur = copy.deepcopy(base)
    if objectstore.is_remote(f):
        cur["files"] = [f]
        cur["description"] = os.path.splitext(os.path.basename(f))[0]
    else:
        cur["files"] = [os.path.abspath(f)]
        cur["description"] = ((sample_name(f) if f.endswith(".bam") else None)
                              or os.path.splitext(os.path.basename(f))[0])
    return cur

def _prep_fastq_input(fs, base):
    for f in fs:
        if not os.path.exists(f) and not objectstore.is_remote(f):
            raise ValueError("Could not find input file: %s" % f)
    cur = copy.deepcopy(base)
    cur["files"] = [os.path.abspath(f) if not objectstore.is_remote(f) else f for f in fs]
    d = os.path.commonprefix([utils.splitext_plus(os.path.basename(f))[0] for f in fs])
    cur["description"] = fastq.rstrip_extra(d)
    return cur

KNOWN_EXTS = {".bam": "bam", ".cram": "bam", ".fq": "fastq",
              ".fastq": "fastq", ".txt": "fastq",
              ".fastq.gz": "fastq", ".fq.gz": "fastq",
              ".txt.gz": "fastq", ".gz": "fastq",
              ".fastq.bz2": "fastq", ".fq.bz2": "fastq",
              ".txt.bz2": "fastq", ".bz2": "fastq"}

def _prep_items_from_base(base, in_files):
    """Prepare a set of configuration items for input files.
    """
    details = []
    in_files = _expand_dirs(in_files, KNOWN_EXTS)
    in_files = _expand_wildcards(in_files)

    for i, (ext, files) in enumerate(itertools.groupby(
            in_files, lambda x: KNOWN_EXTS.get(utils.splitext_plus(x)[-1].lower()))):
        if ext == "bam":
            for f in files:
                details.append(_prep_bam_input(f, i, base))
        elif ext in ["fastq", "fq", "fasta"]:
            files = list(files)
            for fs in fastq.combine_pairs(files):
                details.append(_prep_fastq_input(fs, base))
        else:
            print("Ignoring unexpected input file types %s: %s" % (ext, list(files)))
    return details

def _expand_file(x):
    return os.path.abspath(os.path.normpath(os.path.expanduser(os.path.expandvars(x))))

def _expand_dirs(in_files, known_exts):
    def _is_dir(in_file):
        return os.path.isdir(os.path.expanduser(in_file))
    files, dirs = utils.partition(_is_dir, in_files)
    for dir in dirs:
        for ext in known_exts.keys():
            wildcard = os.path.join(os.path.expanduser(dir), "*" + ext)
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

def name_to_config(template):
    """Read template file into a dictionary to use as base for all samples.

    Handles well-known template names, pulled from GitHub repository and local
    files.
    """
    if objectstore.is_remote(template):
        with objectstore.open(template) as in_handle:
            config = yaml.load(in_handle)
        with objectstore.open(template) as in_handle:
            txt_config = in_handle.read()
    elif os.path.isfile(template):
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
        except (urllib2.HTTPError, urllib2.URLError):
            raise ValueError("Could not find template '%s' locally or in standard templates on GitHub"
                             % template)
    return config, txt_config

def _write_template_config(template_txt, project_name, out_dir):
    config_dir = utils.safe_makedir(os.path.join(out_dir, "config"))
    out_config_file = os.path.join(config_dir, "%s-template.yaml" % project_name)
    with open(out_config_file, "w") as out_handle:
        out_handle.write(template_txt)
    return out_config_file

def _write_config_file(items, global_vars, template, project_name, out_dir,
                       remotes):
    """Write configuration file, adding required top level attributes.
    """
    config_dir = utils.safe_makedir(os.path.join(out_dir, "config"))
    out_config_file = os.path.join(config_dir, "%s.yaml" % project_name)
    out = {"fc_date": datetime.datetime.now().strftime("%Y-%m-%d"),
           "fc_name": project_name,
           "upload": {"dir": "../final"},
           "details": items}
    if remotes.get("base"):
        r_base = objectstore.parse_remote(remotes.get("base"))
        out["upload"]["method"] = r_base.store
        out["upload"]["bucket"] = r_base.bucket
        out["upload"]["folder"] = os.path.join(r_base.key, "final") if r_base.key else "final"
        if r_base.region:
            out["upload"]["region"] = r_base.region
    if global_vars:
        out["globals"] = global_vars
    for k, v in template.iteritems():
        if k not in ["details"]:
            out[k] = v
    if os.path.exists(out_config_file):
        shutil.move(out_config_file,
                    out_config_file + ".bak%s" % datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))
    with open(out_config_file, "w") as out_handle:
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_config_file

def _safe_name(x):
    for prob in [" ", "."]:
        x = x.replace(prob, "_")
    return x

def _set_global_vars(metadata):
    """Identify files used multiple times in metadata and replace with global variables
    """
    fnames = collections.defaultdict(list)
    for sample in metadata.keys():
        for k, v in metadata[sample].items():
            if isinstance(v, basestring) and os.path.isfile(v):
                v = _expand_file(v)
                metadata[sample][k] = v
                fnames[v].append(k)
    global_vars = {}
    # Skip global vars -- more confusing than useful
    # loc_counts = collections.defaultdict(int)
    # global_var_sub = {}
    # for fname, locs in fnames.items():
    #     if len(locs) > 1:
    #         loc_counts[locs[0]] += 1
    #         name = "%s%s" % (locs[0], loc_counts[locs[0]])
    #         global_var_sub[fname] = name
    #         global_vars[name] = fname
    # for sample in metadata.keys():
    #     for k, v in metadata[sample].items():
    #         if isinstance(v, basestring) and v in global_var_sub:
    #             metadata[sample][k] = global_var_sub[v]
    return metadata, global_vars

def _parse_metadata(in_handle):
    """Reads metadata from a simple CSV structured input file.

    samplename,batch,phenotype
    ERR256785,batch1,normal
    """
    metadata = {}
    reader = csv.reader(in_handle)
    while 1:
        header = reader.next()
        if not header[0].startswith("#"):
            break
    keys = [x.strip() for x in header[1:]]
    for sinfo in (x for x in reader if x and not x[0].startswith("#")):
        sinfo = [_strip_and_convert_lists(x) for x in sinfo]
        sample = sinfo[0]
        # sanity check to avoid duplicate rows
        if sample in metadata:
            raise ValueError("Sample %s present multiple times in metadata file.\n"
                             "If you need to specify multiple attributes as a list "
                             "use a semi-colon to separate them on a single line.\n"
                             "https://bcbio-nextgen.readthedocs.org/en/latest/"
                             "contents/configuration.html#automated-sample-configuration\n"
                             "Duplicate line is %s" % (sample, sinfo))
        metadata[sample] = dict(zip(keys, sinfo[1:]))
    metadata, global_vars = _set_global_vars(metadata)
    return metadata, global_vars

def _strip_and_convert_lists(field):
    field = field.strip()
    if "," in field:
        field = [x.strip() for x in field.split(",")]
    return field

def _pname_and_metadata(in_file):
    """Retrieve metadata and project name from the input metadata CSV file.

    Uses the input file name for the project name and for back compatibility,
    accepts the project name as an input, providing no metadata.
    """
    if os.path.isfile(in_file):
        with open(in_file) as in_handle:
            md, global_vars = _parse_metadata(in_handle)
        base = os.path.splitext(os.path.basename(in_file))[0]
        md_file = in_file
    elif objectstore.is_remote(in_file):
        with objectstore.open(in_file) as in_handle:
            md, global_vars = _parse_metadata(in_handle)
        base = os.path.splitext(os.path.basename(in_file))[0]
        md_file = None
    else:
        if in_file.endswith(".csv"):
            raise ValueError("Did not find input metadata file: %s" % in_file)
        base, md, global_vars = _safe_name(os.path.splitext(os.path.basename(in_file))[0]), {}, {}
        md_file = None
    return _safe_name(base), md, global_vars, md_file

def _handle_special_yaml_cases(v):
    """Handle values that pass integer, boolean or list values.
    """
    if ";" in v:
        v = v.split(";")
    elif isinstance(v, list):
        v = v
    else:
        try:
            v = int(v)
        except ValueError:
            if v.lower() == "true":
                v = True
            elif v.lower() == "false":
                    v = False
    return v

def _add_ped_metadata(name, metadata):
    """Add standard PED file attributes into metadata if not present.

    http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
    """
    ignore = set(["-9", "undefined", "unknown", "."])
    def _ped_mapping(x, valmap):
        try:
            x = int(x)
        except ValueError:
            x = -1
        for k, v in valmap.items():
            if k == x:
                return v
        return None
    def _ped_to_gender(x):
        return _ped_mapping(x, {1: "male", 2: "female"})
    def _ped_to_phenotype(x):
        known_phenotypes = set(["unaffected", "affected", "tumor", "normal"])
        if x in known_phenotypes:
            return x
        else:
            return _ped_mapping(x, {1: "unaffected", 2: "affected"})
    def _ped_to_batch(x):
        if x not in ignore and x != "0":
            return x
    with open(metadata["ped"]) as in_handle:
        for line in in_handle:
            parts = line.split("\t")[:6]
            if parts[1] == str(name):
                for index, key, convert_fn in [(4, "sex", _ped_to_gender), (0, "batch", _ped_to_batch),
                                               (5, "phenotype", _ped_to_phenotype)]:
                    val = convert_fn(parts[index])
                    if val is not None and key not in metadata:
                        metadata[key] = val
                break
    return metadata

def _add_metadata(item, metadata, remotes, only_metadata=False):
    """Add metadata information from CSV file to current item.

    Retrieves metadata based on 'description' parsed from input CSV file.
    Adds to object and handles special keys:
    - `description`: A new description for the item. Used to relabel items
       based on the pre-determined description from fastq name or BAM read groups.
    - Keys matching supported names in the algorithm section map
      to key/value pairs there instead of metadata.
    """
    item_md = metadata.get(item["description"],
                           metadata.get(os.path.basename(item["files"][0]),
                                        metadata.get(utils.splitext_plus(os.path.basename(item["files"][0]))[0], {})))
    if remotes.get("region"):
        item["algorithm"]["variant_regions"] = remotes["region"]
    TOP_LEVEL = set(["description", "genome_build", "lane", "vrn_files", "files", "analysis"])
    keep_sample = True
    if len(item_md) > 0:
        if "metadata" not in item:
            item["metadata"] = {}
        for k, v in item_md.iteritems():
            if v:
                if k in TOP_LEVEL:
                    item[k] = v
                elif k in run_info.ALGORITHM_KEYS:
                    v = _handle_special_yaml_cases(v)
                    item["algorithm"][k] = v
                else:
                    v = _handle_special_yaml_cases(v)
                    item["metadata"][k] = v
    elif len(metadata) > 0:
        warn = "Dropped sample" if only_metadata else "Added minimal sample information"
        print "WARNING: %s: metadata not found for %s, %s" % (warn, item["description"],
                                                              os.path.basename(item["files"][0]))
        keep_sample = not only_metadata
    if tz.get_in(["metadata", "ped"], item):
        item["metadata"] = _add_ped_metadata(item["description"], item["metadata"])
    return item if keep_sample else None

def _retrieve_remote(fnames):
    """Retrieve remote inputs found in the same bucket as the template or metadata files.
    """
    for fname in fnames:
        if objectstore.is_remote(fname):
            inputs = []
            regions = []
            remote_base = os.path.dirname(fname)
            for rfname in objectstore.list(remote_base):
                if rfname.endswith(tuple(KNOWN_EXTS.keys())):
                    inputs.append(rfname)
                elif rfname.endswith((".bed", ".bed.gz")):
                    regions.append(rfname)
            return {"base": remote_base,
                    "inputs": inputs,
                    "region": regions[0] if len(regions) == 1 else None}
    return {}

def _convert_to_relpaths(data, work_dir):
    """Convert absolute paths in the input data to relative paths to the work directory.
    """
    work_dir = os.path.abspath(work_dir)
    data["files"] = [os.path.relpath(f, work_dir) for f in data["files"]]
    for topk in ["metadata", "algorithm"]:
        for k, v in data[topk].items():
            if isinstance(v, basestring) and os.path.isfile(v) and os.path.isabs(v):
                data[topk][k] = os.path.relpath(v, work_dir)
    return data

def _check_all_metadata_found(metadata, items):
    """Print warning if samples in CSV file are missing in folder"""
    for name in metadata:
        seen = False
        for sample in items:
            if sample['files'][0].find(name) > -1:
                seen = True
        if not seen:
            print "WARNING: sample not found %s" % name

def setup(args):
    template, template_txt = name_to_config(args.template)
    base_item = template["details"][0]
    project_name, metadata, global_vars, md_file = _pname_and_metadata(args.metadata)
    remotes = _retrieve_remote([args.metadata, args.template])
    inputs = args.input_files + remotes.get("inputs", [])
    if hasattr(args, "systemconfig") and args.systemconfig and hasattr(args, "integrations"):
        config, _ = config_utils.load_system_config(args.systemconfig)
        for iname, retriever in args.integrations.items():
            if iname in config:
                inputs += retriever.get_files(metadata, config[iname])
    raw_items = [_add_metadata(item, metadata, remotes, args.only_metadata)
                 for item in _prep_items_from_base(base_item, inputs)]
    items = [x for x in raw_items if x]
    _check_all_metadata_found(metadata, items)
    out_dir = os.path.join(os.getcwd(), project_name)
    work_dir = utils.safe_makedir(os.path.join(out_dir, "work"))
    if hasattr(args, "relpaths") and args.relpaths:
        items = [_convert_to_relpaths(x, work_dir) for x in items]
    out_config_file = _write_template_config(template_txt, project_name, out_dir)
    if md_file:
        shutil.copyfile(md_file, os.path.join(out_dir, "config", os.path.basename(md_file)))
    if len(items) == 0:
        print
        print "Template configuration file created at: %s" % out_config_file
        print "Edit to finalize custom options, then prepare full sample config with:"
        print "  bcbio_nextgen.py -w template %s %s sample1.bam sample2.fq" % \
            (out_config_file, project_name)
    else:
        out_config_file = _write_config_file(items, global_vars, template, project_name, out_dir,
                                             remotes)
        print
        print "Configuration file created at: %s" % out_config_file
        print "Edit to finalize and run with:"
        print "  cd %s" % work_dir
        print "  bcbio_nextgen.py ../config/%s" % os.path.basename(out_config_file)
        if remotes.get("base"):
            remote_path = os.path.join(remotes["base"], os.path.basename(out_config_file))
            s3.upload_file_boto(out_config_file, remote_path)
            print "Also uploaded to AWS S3 in %s" % remotes["base"]
            print "Run directly with bcbio_vm.py run %s" % remote_path
