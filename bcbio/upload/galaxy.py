"""Move files to local Galaxy upload directory and add to Galaxy Data Libraries.

Required configurable variables in upload:
  dir
"""
import collections
import os
import shutil
import time

from bcbio import utils
from bcbio.log import logger
from bcbio.upload import filesystem
from bcbio.pipeline import qcsummary

# Avoid bioblend import errors, raising at time of use
try:
    import bioblend
    from bioblend.galaxy import GalaxyInstance
    import simplejson
except ImportError:
    GalaxyInstance, bioblend, simplejson = None, None, None

def update_file(finfo, sample_info, config):
    """Update file in Galaxy data libraries.
    """
    if GalaxyInstance is None:
        raise ImportError("Could not import bioblend.galaxy")
    if "dir" not in config:
        raise ValueError("Galaxy upload requires `dir` parameter in config specifying the "
                         "shared filesystem path to move files to.")
    if "outputs" in config:
        _galaxy_tool_copy(finfo, config["outputs"])
    else:
        _galaxy_library_upload(finfo, sample_info, config)

def _galaxy_tool_copy(finfo, outputs):
    """Copy information directly to pre-defined outputs from a Galaxy tool.

    XXX Needs generalization
    """
    tool_map = {"align": "bam", "variants": "vcf.gz"}
    for galaxy_key, finfo_type in tool_map.items():
        if galaxy_key in outputs and finfo.get("type") == finfo_type:
            shutil.copy(finfo["path"], outputs[galaxy_key])

def _galaxy_library_upload(finfo, sample_info, config):
    """Upload results to galaxy library.
    """
    folder_name = "%s_%s" % (config["fc_date"], config["fc_name"])
    storage_dir = utils.safe_makedir(os.path.join(config["dir"], folder_name))
    if finfo.get("type") == "directory":
        storage_file = None
        if finfo.get("ext") == "qc":
            pdf_file = qcsummary.prep_pdf(finfo["path"], config)
            if pdf_file:
                finfo["path"] = pdf_file
                finfo["type"] = "pdf"
                storage_file = filesystem.copy_finfo(finfo, storage_dir, pass_uptodate=True)
    else:
        storage_file = filesystem.copy_finfo(finfo, storage_dir, pass_uptodate=True)
    if "galaxy_url" in config and "galaxy_api_key" in config:
        galaxy_url = config["galaxy_url"]
        if not galaxy_url.endswith("/"):
            galaxy_url += "/"
        gi = GalaxyInstance(galaxy_url, config["galaxy_api_key"])
    else:
        raise ValueError("Galaxy upload requires `galaxy_url` and `galaxy_api_key` in config")
    if storage_file and sample_info and not finfo.get("index", False) and not finfo.get("plus", False):
        _to_datalibrary_safe(storage_file, gi, folder_name, sample_info, config)

def _to_datalibrary_safe(fname, gi, folder_name, sample_info, config):
    """Upload with retries for intermittent JSON failures.
    """
    num_tries = 0
    max_tries = 5
    while 1:
        try:
            _to_datalibrary(fname, gi, folder_name, sample_info, config)
            break
        except (simplejson.scanner.JSONDecodeError, bioblend.galaxy.client.ConnectionError) as e:
            num_tries += 1
            if num_tries > max_tries:
                raise
            print "Retrying upload, failed with:", str(e)
            time.sleep(5)

def _to_datalibrary(fname, gi, folder_name, sample_info, config):
    """Upload a file to a Galaxy data library in a project specific folder.
    """
    library = _get_library(gi, sample_info, config)
    libitems = gi.libraries.show_library(library.id, contents=True)
    folder = _get_folder(gi, folder_name, library, libitems)
    _file_to_folder(gi, fname, sample_info, libitems, library, folder)

def _file_to_folder(gi, fname, sample_info, libitems, library, folder):
    """Check if file exists on Galaxy, if not upload to specified folder.
    """
    full_name = os.path.join(folder["name"], os.path.basename(fname))

    # Handle VCF: Galaxy reports VCF files without the gzip extension
    file_type = "vcf_bgzip" if full_name.endswith(".vcf.gz") else "auto"
    if full_name.endswith(".vcf.gz"):
        full_name = full_name.replace(".vcf.gz", ".vcf")

    for item in libitems:
        if item["name"] == full_name:
            return item
    logger.info("Uploading to Galaxy library '%s': %s" % (library.name, full_name))
    return gi.libraries.upload_from_galaxy_filesystem(str(library.id), fname, folder_id=str(folder["id"]),
                                                      link_data_only="link_to_files",
                                                      dbkey=sample_info["genome_build"],
                                                      file_type=file_type,
                                                      roles=str(library.roles) if library.roles else None)

def _get_folder(gi, folder_name, library, libitems):
    """Retrieve or create a folder inside the library with the specified name.
    """
    for item in libitems:
        if item["type"] == "folder" and item["name"] == "/%s" % folder_name:
            return item
    return gi.libraries.create_folder(library.id, folder_name)[0]

GalaxyLibrary = collections.namedtuple("GalaxyLibrary", ["id", "name", "roles"])

def _get_library(gi, sample_info, config):
    """Retrieve the appropriate data library for the current user.
    """
    galaxy_lib = sample_info.get("galaxy_library",
                                 config.get("galaxy_library"))
    role = sample_info.get("galaxy_role",
                           config.get("galaxy_role"))
    if galaxy_lib:
        return _get_library_from_name(gi, galaxy_lib, role, sample_info, create=True)
    elif config.get("private_libs") or config.get("lab_association") or config.get("researcher"):
        return _library_from_nglims(gi, sample_info, config)
    else:
        raise ValueError("No Galaxy library specified for sample: %s" %
                         sample_info["description"])

def _get_library_from_name(gi, name, role, sample_info, create=False):
    for lib in gi.libraries.get_libraries():
        if lib["name"].lower() == name.lower() and not lib.get("deleted", False):
            return GalaxyLibrary(lib["id"], lib["name"], role)
    if create and name:
        logger.info("Creating Galaxy library: '%s'" % name)
        lib = gi.libraries.create_library(name)
        librole = str(gi.users.get_current_user()["id"] if not role else role)
        try:
            gi.libraries.set_library_permissions(str(lib["id"]), librole, librole, librole, librole)
        # XXX Returns error on Galaxy side but seems to work -- ugly
        except:
            pass
        return GalaxyLibrary(lib["id"], lib["name"], role)
    else:
        raise ValueError("Could not find Galaxy library matching '%s' for sample %s" %
                         (name, sample_info["description"]))

def _library_from_nglims(gi, sample_info, config):
    """Retrieve upload library from nglims specified user libraries.
    """
    names = [config.get(x, "").strip() for x in ["lab_association", "researcher"]
             if config.get(x)]
    for name in names:
        for ext in ["sequencing", "lab"]:
            check_name = "%s %s" % (name.split()[0], ext)
            try:
                return _get_library_from_name(gi, check_name, None, sample_info)
            except ValueError:
                pass
    check_names = set([x.lower() for x in names])
    for libname, role in config["private_libs"]:
        # Try to find library for lab or rsearcher
        if libname.lower() in check_names:
            return _get_library_from_name(gi, libname, role, sample_info)
    # default to first private library if available
    if len(config.get("private_libs", [])) > 0:
        libname, role = config["private_libs"][0]
        return _get_library_from_name(gi, libname, role, sample_info)
    # otherwise use the lab association or researcher name
    elif len(names) > 0:
        return _get_library_from_name(gi, names[0], None, sample_info, create=True)
    else:
        raise ValueError("Could not find Galaxy library for sample %s" % sample_info["description"])
