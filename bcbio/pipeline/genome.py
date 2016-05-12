"""Read genome build configurations from Galaxy *.loc and bcbio-nextgen resource files.
"""
import ConfigParser
import glob
import os
import sys
from xml.etree import ElementTree

import toolz as tz
import yaml

from bcbio import utils
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
        resources = yaml.load(in_handle)

    def resource_file_path(x):
        if isinstance(x, basestring) and os.path.exists(os.path.join(base_dir, x)):
            return os.path.normpath(os.path.join(base_dir, x))
        return x

    cleaned = utils.dictapply(resources, resource_file_path)
    return _ensure_annotations(cleaned, data)

def _ensure_annotations(resources, data):
    """Prepare any potentially missing annotations for downstream processing in a local directory.
    """
    out_dir = utils.safe_makedir(os.path.join(tz.get_in(["dirs", "work"], data),
                                              "inputs", "data", "annotations"))
    transcript_gff = tz.get_in(["rnaseq", "transcripts"], resources)
    if transcript_gff and utils.file_exists(transcript_gff):
        resources["rnaseq"]["gene_bed"] = gtf.gtf_to_bed(transcript_gff, out_dir)
    return resources

# ## Utilities

def abs_file_paths(xs, base_dir=None, ignore_keys=None):
    """Normalize any file paths found in a subdirectory of configuration input.
    """
    ignore_keys = set([]) if ignore_keys is None else set(ignore_keys)
    if base_dir is None:
        base_dir = os.getcwd()
    orig_dir = os.getcwd()
    os.chdir(base_dir)
    input_dir = os.path.join(base_dir, "inputs")
    if isinstance(xs, dict):
        out = {}
        for k, v in xs.iteritems():
            if k not in ignore_keys and v and isinstance(v, basestring):
                if v.lower() == "none":
                    out[k] = None
                elif os.path.exists(v) or objectstore.is_remote(v):
                    dl = objectstore.download(v, input_dir)
                    if dl:
                        out[k] = os.path.normpath(os.path.join(base_dir, dl))
                    else:
                        out[k] = v
                else:
                    out[k] = v
            else:
                out[k] = v
    elif isinstance(xs, basestring):
        if os.path.exists(xs) or objectstore.is_remote(xs):
            dl = objectstore.download(xs, input_dir)
            if dl:
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
        logger.info("Downloading %s %s from AWS" % (genome_build, name))
        cur_ref = download_prepped_genome(genome_build, data, name, need_remap)
    # allow multiple references in a file and use the most recently added
    else:
        cur_ref = refs[-1]
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
            if os.path.isdir(base):
                indexes = glob.glob(os.path.join(base, "*"))
            elif name != "samtools":
                indexes = glob.glob("%s*" % utils.splitext_plus(base)[0])
            else:
                indexes = []
            out[name_remap.get(name, name)] = {}
            if os.path.exists(base) and os.path.isfile(base):
                out[name_remap.get(name, name)]["base"] = base
            if indexes:
                out[name_remap.get(name, name)]["indexes"] = indexes
        # add additional indices relative to the base
        if tz.get_in(["fasta", "base"], out):
            ref_dir, ref_filebase = os.path.split(out["fasta"]["base"])
            out["rtg"] = os.path.normpath(os.path.join(ref_dir, os.path.pardir, "rtg",
                                                       "%s.sdf" % (os.path.splitext(ref_filebase)[0])))
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

# ## Retrieve pre-prepared genomes

REMAP_NAMES = {"tophat2": ["bowtie2"],
               "samtools": ["rtg", "seq"]}
INPLACE_INDEX = {"star": star.index}

def download_prepped_genome(genome_build, data, name, need_remap, out_dir=None):
    """Get a pre-prepared genome from S3, unpacking it locally.

    Supports runs on AWS where we can retrieve the resources on demand. Upgrades
    GEMINI in place if installed inside a Docker container with the biological data.
    GEMINI install requires write permissions to standard data directories -- works
    on AWS but not generalizable elsewhere.
    """
    from bcbio.variation import population
    from bcbio import install
    if not out_dir:
        out_dir = utils.safe_makedir(os.path.join(tz.get_in(["dirs", "work"], data),
                                                  "inputs", "data", "genomes"))
    for target in REMAP_NAMES.get(name, [name]):
        ref_dir = os.path.join(out_dir, genome_build, target)
        if not os.path.exists(ref_dir):
            if target in INPLACE_INDEX:
                ref_file = glob.glob(os.path.normpath(os.path.join(ref_dir, os.pardir, "seq", "*.fa")))[0]
                INPLACE_INDEX[target](ref_file, ref_dir, data)
            else:
                # XXX Currently only supports genomes from S3 us-east-1 bucket.
                # Need to assess how slow this is from multiple regions and generalize to non-AWS.
                fname = objectstore.BIODATA_INFO["s3"].format(build=genome_build, target=target)
                try:
                    objectstore.connect(fname)
                except:
                    raise ValueError("Could not find reference genome file %s %s" % (genome_build, name))
                with utils.chdir(out_dir):
                    cmd = objectstore.cl_input(fname, unpack=False, anonpipe=False) + " | pigz -d -c | tar -xvp"
                    do.run(cmd.format(**locals()), "Download pre-prepared genome data: %s" % genome_build)
    ref_file = glob.glob(os.path.normpath(os.path.join(ref_dir, os.pardir, "seq", "*.fa")))[0]
    if data.get("genome_build"):
        gresources = get_resources(data["genome_build"], ref_file, data)
        if data.get("files") and population.do_db_build([data], need_bam=False, gresources=gresources):
            # symlink base GEMINI directory to work directory, avoiding write/space issues
            out_gemini_dir = utils.safe_makedir(os.path.join(os.path.dirname(ref_dir), "gemini_data"))
            orig_gemini_dir = install.get_gemini_dir()
            # Remove empty initial directory created by installer
            if os.path.isdir(orig_gemini_dir) and len(os.listdir(orig_gemini_dir)) == 0:
                if os.path.islink(orig_gemini_dir):
                    os.remove(orig_gemini_dir)
                else:
                    os.rmdir(orig_gemini_dir)
            if not os.path.exists(orig_gemini_dir):
                os.symlink(out_gemini_dir, orig_gemini_dir)
            cmd = [os.path.join(os.path.dirname(sys.executable), "gemini"), "update", "--dataonly"]
            do.run(cmd, "Download GEMINI data")
    genome_dir = os.path.join(out_dir, genome_build)
    genome_build = genome_build.replace("-test", "")
    if need_remap or name == "samtools":
        return os.path.join(genome_dir, "seq", "%s.fa" % genome_build)
    else:
        ref_dir = os.path.join(genome_dir, REMAP_NAMES.get(name, [name])[-1])
        base_name = os.path.commonprefix(os.listdir(ref_dir))
        while base_name.endswith("."):
            base_name = base_name[:-1]
        return os.path.join(ref_dir, base_name)
