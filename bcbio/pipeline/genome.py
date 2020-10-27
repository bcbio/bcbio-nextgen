"""Read genome build configurations from Galaxy *.loc and bcbio-nextgen resource files.
"""
from six.moves import configparser
import glob
import os
import sys
from xml.etree import ElementTree

import six
import toolz as tz
import yaml

from bcbio import utils
from bcbio.cwl import cwlutils
from bcbio.distributed import objectstore
from bcbio.log import logger
from bcbio.ngsalign import star
from bcbio.pipeline import alignment
from bcbio.provenance import do
from bcbio.rnaseq import gtf

# ## bcbio-nextgen genome resource files

def get_resources(genome, ref_file, data):
    """Retrieve genome information from a genome-references.yaml file.
    """
    base_dir = os.path.normpath(os.path.dirname(ref_file))
    resource_file = os.path.join(base_dir, "%s-resources.yaml" % genome.replace("-test", ""))
    if not os.path.exists(resource_file):
        raise IOError("Did not find resource file for %s: %s\n"
                      "To update bcbio_nextgen.py with genome resources for standard builds, run:\n"
                      "bcbio_nextgen.py upgrade -u skip"
                      % (genome, resource_file))
    with open(resource_file) as in_handle:
        resources = yaml.safe_load(in_handle)

    def resource_file_path(x):
        if isinstance(x, six.string_types) and os.path.exists(os.path.join(base_dir, x)):
            return os.path.normpath(os.path.join(base_dir, x))
        return x
    cleaned = utils.dictapply(resources, resource_file_path)
    return ensure_annotations(cleaned, data)

def add_required_resources(resources):
    """Add default or empty values for required resources referenced in CWL
    """
    required = [["variation", "cosmic"], ["variation", "clinvar"], ["variation", "dbsnp"],
                ["variation", "lcr"], ["variation", "polyx"],
                ["variation", "encode_blacklist"], ["variation", "gc_profile"],
                ["variation", "germline_het_pon"],
                ["variation", "train_hapmap"], ["variation", "train_indels"],
                ["variation", "editing"], ["variation", "exac"], ["variation", "esp"],
                ["variation", "gnomad_exome"],
                ["variation", "1000g"], ["aliases", "human"]]
    for key in required:
        if not tz.get_in(key, resources):
            resources = tz.update_in(resources, key, lambda x: None)
    return resources

def ensure_annotations(resources, data):
    """Prepare any potentially missing annotations for downstream processing in a local directory.
    """
    transcript_gff = tz.get_in(["rnaseq", "transcripts"], resources)
    if transcript_gff and utils.file_exists(transcript_gff):
        out_dir = os.path.join(tz.get_in(["dirs", "work"], data),
                               "inputs", "data", "annotations")
        resources["rnaseq"]["gene_bed"] = gtf.gtf_to_bed(transcript_gff, out_dir)
    return resources

# ## Utilities

def abs_file_paths(xs, base_dir=None, ignore_keys=None, fileonly_keys=None, cur_key=None,
                   do_download=True):
    """Normalize any file paths found in a subdirectory of configuration input.

    base_dir -- directory to normalize relative paths to
    ignore_keys -- algorithm key names to ignore normalize for (keywords, not files/directories)
    fileonly_keys -- algorithm key names to only expand files (not directories)
    cur_key -- current key when calling recursively
    """
    ignore_keys = set([]) if ignore_keys is None else set(ignore_keys)
    fileonly_keys = set([]) if fileonly_keys is None else set(fileonly_keys)
    if base_dir is None:
        base_dir = os.getcwd()
    orig_dir = os.getcwd()
    os.chdir(base_dir)
    input_dir = os.path.join(base_dir, "inputs")
    if isinstance(xs, dict):
        out = {}
        for k, v in xs.items():
            if k not in ignore_keys and v and isinstance(v, six.string_types):
                if v.lower() == "none":
                    out[k] = None
                else:
                    out[k] = abs_file_paths(v, base_dir, ignore_keys, fileonly_keys, k, do_download=do_download)
            elif isinstance(v, (list, tuple)):
                out[k] = [abs_file_paths(x, base_dir, ignore_keys, fileonly_keys, k, do_download=do_download)
                          for x in v]
            else:
                out[k] = v
    elif isinstance(xs, six.string_types):
        if os.path.exists(xs) or (do_download and objectstore.is_remote(xs)):
            dl = objectstore.download(xs, input_dir)
            if dl and cur_key not in ignore_keys and not (cur_key in fileonly_keys and not os.path.isfile(dl)):
                out = os.path.normpath(os.path.join(base_dir, dl))
            else:
                out = xs
        else:
            out = xs
    else:
        out = xs
    os.chdir(orig_dir)
    return out

# ## Galaxy integration -- *.loc files

def _get_galaxy_loc_file(name, galaxy_dt, ref_dir, galaxy_base):
    """Retrieve Galaxy *.loc file for the given reference/aligner name.

    First tries to find an aligner specific *.loc file. If not defined
    or does not exist, then we need to try and remap it from the
    default reference file
    """
    if "file" in galaxy_dt and os.path.exists(os.path.join(galaxy_base, galaxy_dt["file"])):
        loc_file = os.path.join(galaxy_base, galaxy_dt["file"])
        need_remap = False
    elif alignment.TOOLS[name].galaxy_loc_file is None:
        loc_file = os.path.join(ref_dir, alignment.BASE_LOCATION_FILE)
        need_remap = True
    else:
        loc_file = os.path.join(ref_dir, alignment.TOOLS[name].galaxy_loc_file)
        need_remap = False
    if not os.path.exists(loc_file):
        loc_file = os.path.join(ref_dir, alignment.BASE_LOCATION_FILE)
        need_remap = True
    return loc_file, need_remap

def _galaxy_loc_iter(loc_file, galaxy_dt, need_remap=False):
    """Iterator returning genome build and references from Galaxy *.loc file.
    """
    if "column" in galaxy_dt:
        dbkey_i = galaxy_dt["column"].index("dbkey")
        path_i = galaxy_dt["column"].index("path")
    else:
        dbkey_i = None
    if os.path.exists(loc_file):
        with open(loc_file) as in_handle:
            for line in in_handle:
                if line.strip() and not line.startswith("#"):
                    parts = [x.strip() for x in line.strip().split("\t")]
                    # Detect and report spaces instead of tabs
                    if len(parts) == 1:
                        parts = [x.strip() for x in line.strip().split(" ") if x.strip()]
                        if len(parts) > 1:
                            raise IOError("Galaxy location file uses spaces instead of "
                                          "tabs to separate fields: %s" % loc_file)
                    if dbkey_i is not None and not need_remap:
                        dbkey = parts[dbkey_i]
                        cur_ref = parts[path_i]
                    else:
                        if parts[0] == "index":
                            parts = parts[1:]
                        dbkey = parts[0]
                        cur_ref = parts[-1]
                    yield (dbkey, cur_ref)

def _get_ref_from_galaxy_loc(name, genome_build, loc_file, galaxy_dt, need_remap,
                             galaxy_config, data):
    """Retrieve reference genome file from Galaxy *.loc file.

    Reads from tool_data_table_conf.xml information for the index if it
    exists, otherwise uses heuristics to find line based on most common setups.
    """
    refs = [ref for dbkey, ref in _galaxy_loc_iter(loc_file, galaxy_dt, need_remap)
            if dbkey == genome_build]
    remap_fn = alignment.TOOLS[name].remap_index_fn
    need_remap = remap_fn is not None
    if len(refs) == 0:
        raise ValueError("Did not find genome build %s in bcbio installation: %s" %
                         (genome_build, os.path.normpath(loc_file)))
    else:
        cur_ref = refs[-1]
    # Find genome directory and check for packed wf tarballs
    cur_ref_norm = os.path.normpath(utils.add_full_path(cur_ref, galaxy_config["tool_data_path"]))
    base_dir_i = cur_ref_norm.find("/%s/" % genome_build)
    base_dir = os.path.join(cur_ref_norm[:base_dir_i], genome_build)
    for tarball in glob.glob(os.path.join(base_dir, "*-wf.tar.gz")):
        cwlutils.unpack_tarballs(tarball, {"dirs": {"work": base_dir}}, use_subdir=False)
    if need_remap:
        assert remap_fn is not None, "%s requires remapping function from base location file" % name
        cur_ref = os.path.normpath(utils.add_full_path(cur_ref, galaxy_config["tool_data_path"]))
        cur_ref = remap_fn(os.path.abspath(cur_ref))
    return cur_ref

def _get_galaxy_tool_info(galaxy_base):
    """Retrieve Galaxy tool-data information from defaults or galaxy config file.
    """
    ini_file = os.path.join(galaxy_base, "universe_wsgi.ini")
    info = {"tool_data_table_config_path": os.path.join(galaxy_base, "tool_data_table_conf.xml"),
            "tool_data_path": os.path.join(galaxy_base, "tool-data")}
    config = configparser.ConfigParser()
    config.read(ini_file)
    if "app:main" in config.sections():
        for option in config.options("app:main"):
            if option in info:
                info[option] = os.path.join(galaxy_base, config.get("app:main", option))
    return info

def _get_galaxy_data_table(name, dt_config_file):
    """Parse data table config file for details on tool *.loc location and columns.
    """
    out = {}
    if os.path.exists(dt_config_file):
        tdtc = ElementTree.parse(dt_config_file)
        for t in tdtc.getiterator("table"):
            if t.attrib.get("name", "") in [name, "%s_indexes" % name]:
                out["column"] = [x.strip() for x in t.find("columns").text.split(",")]
                out["file"] = t.find("file").attrib.get("path", "")
    return out

def get_refs(genome_build, aligner, galaxy_base, data):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    out = {}
    name_remap = {"samtools": "fasta"}
    if genome_build:
        galaxy_config = _get_galaxy_tool_info(galaxy_base)
        for name in [x for x in ("samtools", aligner) if x]:
            galaxy_dt = _get_galaxy_data_table(name, galaxy_config["tool_data_table_config_path"])
            loc_file, need_remap = _get_galaxy_loc_file(name, galaxy_dt, galaxy_config["tool_data_path"],
                                                        galaxy_base)
            cur_ref = _get_ref_from_galaxy_loc(name, genome_build, loc_file, galaxy_dt, need_remap,
                                               galaxy_config, data)
            base = os.path.normpath(utils.add_full_path(cur_ref, galaxy_config["tool_data_path"]))
            # Expand directories unless we are an aligner like minimap2 that uses the seq directory
            if os.path.isdir(base) and not (need_remap and os.path.basename(base) == "seq"):
                indexes = sorted(glob.glob(os.path.join(base, "*")))
            elif name != "samtools":
                indexes = sorted(glob.glob("%s*" % utils.splitext_plus(base)[0]))
            else:
                indexes = []
            name = name_remap.get(name, name)
            out[name] = {}
            if os.path.exists(base) and os.path.isfile(base):
                out[name]["base"] = base
            if indexes:
                out[name]["indexes"] = indexes
            # For references, add compressed inputs and indexes if they exist
            if name == "fasta" and "base" in out[name] and os.path.exists(out[name]["base"] + ".gz"):
                indexes = [out[name]["base"] + ".gz.fai", out[name]["base"] + ".gz.gzi",
                           utils.splitext_plus(out[name]["base"])[0] + ".dict"]
                out[name + "gz"] = {"base": out[name]["base"] + ".gz",
                                    "indexes": [x for x in indexes if os.path.exists(x)]}
        # add additional indices relative to the base
        if tz.get_in(["fasta", "base"], out):
            ref_dir, ref_filebase = os.path.split(out["fasta"]["base"])
            rtg_dir = os.path.normpath(os.path.join(ref_dir, os.path.pardir, "rtg",
                                                    "%s.sdf" % (os.path.splitext(ref_filebase)[0])))
            out["rtg"] = {"base": os.path.join(rtg_dir, "mainIndex"),
                          "indexes": [x for x in glob.glob(os.path.join(rtg_dir, "*"))
                                      if not x.endswith("/mainIndex")]}
            twobit = os.path.normpath(os.path.join(ref_dir, os.path.pardir, "ucsc",
                                                   "%s.2bit" % (os.path.splitext(ref_filebase)[0])))
            if os.path.exists(twobit):
                out["twobit"] = twobit
    return out

def get_builds(galaxy_base):
    """Retrieve configured genome builds and reference files, using Galaxy configuration files.

    Allows multiple dbkey specifications in the same file, using the most recently added.
    """
    name = "samtools"
    galaxy_config = _get_galaxy_tool_info(galaxy_base)
    galaxy_dt = _get_galaxy_data_table(name, galaxy_config["tool_data_table_config_path"])
    loc_file, need_remap = _get_galaxy_loc_file(name, galaxy_dt, galaxy_config["tool_data_path"],
                                                galaxy_base)
    assert not need_remap, "Should not need to remap reference files"
    fnames = {}
    for dbkey, fname in _galaxy_loc_iter(loc_file, galaxy_dt):
        fnames[dbkey] = fname
    out = []
    for dbkey in sorted(fnames.keys()):
        out.append((dbkey, fnames[dbkey]))
    return out
