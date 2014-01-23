"""Read genome build configurations from Galaxy *.loc and bcbio-nextgen resource files.
"""
import ConfigParser
import os
from xml.etree import ElementTree

import yaml

from bcbio import utils
from bcbio.pipeline import alignment

# ## bcbio-nextgen genome resource files

def get_resources(genome, ref_file):
    """Retrieve genome information from a genome-references.yaml file.
    """
    base_dir = os.path.normpath(os.path.dirname(ref_file))
    resource_file = os.path.join(base_dir, "%s-resources.yaml" % genome)
    if not os.path.exists(resource_file):
        raise IOError("Did not find resource file for %s: %s\n"
                      "To update bcbio_nextgen.py with genome resources for standard builds, run:\n"
                      "bcbio_nextgen.py upgrade -u skip"
                      % (genome, resource_file))
    with open(resource_file) as in_handle:
        resources = yaml.load(in_handle)

    def resource_file_path(x):
        if isinstance(x, basestring) and os.path.exists(os.path.join(base_dir, x)):
            return os.path.normpath(os.path.join(base_dir, x))
        return x

    return utils.dictapply(resources, resource_file_path)

# ## Utilities


def abs_file_paths(xs, base_dir=None, ignore_keys=None):
    """Normalize any file paths found in a subdirectory of configuration input.
    """
    ignore_keys = set([]) if ignore_keys is None else set(ignore_keys)
    if not isinstance(xs, dict):
        return xs
    if base_dir is None:
        base_dir = os.getcwd()
    orig_dir = os.getcwd()
    os.chdir(base_dir)
    out = {}
    for k, v in xs.iteritems():
        if k not in ignore_keys and v and isinstance(v, basestring) and os.path.exists(v):
            out[k] = os.path.normpath(os.path.join(base_dir, v))
        else:
            out[k] = v
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
    with open(loc_file) as in_handle:
        for line in in_handle:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
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

def _get_ref_from_galaxy_loc(name, genome_build, loc_file, galaxy_dt, need_remap):
    """Retrieve reference genome file from Galaxy *.loc file.

    Reads from tool_data_table_conf.xml information for the index if it
    exists, otherwise uses heuristics to find line based on most common setups.
    """
    refs = [ref for dbkey, ref in _galaxy_loc_iter(loc_file, galaxy_dt, need_remap)
            if dbkey == genome_build]
    if len(refs) == 0:
        raise IndexError("Genome %s not found in %s" % (genome_build, loc_file))
    elif len(refs) > 1:
        raise IndexError("Genome %s found multiple times in %s" % (genome_build, loc_file))
    else:
        cur_ref = refs[0]
    if need_remap:
        remap_fn = alignment.TOOLS[name].remap_index_fn
        assert remap_fn is not None, "%s requires remapping function from base location file" % name
        cur_ref = remap_fn(cur_ref)
    return cur_ref

def _get_galaxy_tool_info(galaxy_base):
    """Retrieve Galaxy tool-data information from defaults or galaxy config file.
    """
    ini_file = os.path.join(galaxy_base, "universe_wsgi.ini")
    info = {"tool_data_table_config_path": os.path.join(galaxy_base, "tool_data_table_conf.xml"),
            "tool_data_path": os.path.join(galaxy_base, "tool-data")}
    config = ConfigParser.ConfigParser()
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

def get_refs(genome_build, aligner, galaxy_base):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    if not genome_build:
        return (None, None)
    galaxy_config = _get_galaxy_tool_info(galaxy_base)
    out_info = []
    for name in [aligner, "samtools"]:
        if not name:
            out_info.append(None)
            continue
        galaxy_dt = _get_galaxy_data_table(name, galaxy_config["tool_data_table_config_path"])
        loc_file, need_remap = _get_galaxy_loc_file(name, galaxy_dt, galaxy_config["tool_data_path"],
                                                    galaxy_base)
        cur_ref = _get_ref_from_galaxy_loc(name, genome_build, loc_file, galaxy_dt, need_remap)
        out_info.append(utils.add_full_path(cur_ref, galaxy_config["tool_data_path"]))

    if len(out_info) != 2:
        raise ValueError("Did not find genome reference for %s %s" %
                         (genome_build, aligner))
    else:
        return tuple(out_info)

def get_builds(galaxy_base):
    """Retrieve configured genome builds and reference files, using Galaxy configuration files.
    """
    name = "samtools"
    galaxy_config = _get_galaxy_tool_info(galaxy_base)
    galaxy_dt = _get_galaxy_data_table(name, galaxy_config["tool_data_table_config_path"])
    loc_file, need_remap = _get_galaxy_loc_file(name, galaxy_dt, galaxy_config["tool_data_path"],
                                                galaxy_base)
    assert not need_remap, "Should not need to remap reference files"
    return _galaxy_loc_iter(loc_file, galaxy_dt)

