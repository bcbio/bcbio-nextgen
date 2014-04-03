"""Move files to local Galaxy upload directory and add to Galaxy Data Libraries.

Required configurable variables in upload:
  dir
"""
import collections
import os

from bcbio import utils
from bcbio.log import logger
from bcbio.upload import filesystem

# Avoid bioblend import errors, raising at time of use
try:
    from bioblend.galaxy import GalaxyInstance
except ImportError:
    GalaxyInstance = None

def update_file(finfo, sample_info, config):
    """Update file in Galaxy data libraries.
    """
    if GalaxyInstance is None:
        raise ImportError("Could not import bioblend.galaxy")
    if "dir" not in config:
        raise ValueError("Galaxy upload requires `dir` parameter in config specifying the "
                         "shared filesystem path to move files to.")
    folder_name = "%s_%s" % (config["fc_date"], config["fc_name"])
    storage_dir = utils.safe_makedir(os.path.join(config["dir"], folder_name))
    storage_file = (filesystem.copy_finfo(finfo, storage_dir, pass_uptodate=True)
                    if finfo.get("type") != "directory" else None)
    if "galaxy_url" in config and "galaxy_api_key" in config:
        galaxy_url = config["galaxy_url"]
        if not galaxy_url.endswith("/"):
            galaxy_url += "/"
        gi = GalaxyInstance(galaxy_url, config["galaxy_api_key"])
    else:
        raise ValueError("Galaxy upload requires `galaxy_url` and `galaxy_api_key` in config")
    if storage_file and sample_info and not finfo.get("index", False):
        _to_datalibrary(storage_file, gi, folder_name, sample_info, config)

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
    for item in libitems:
        if item["name"] == full_name:
            return item
    logger.info("Uploading to Galaxy library '%s': %s" % (library.name, full_name))
    return gi.libraries.upload_from_galaxy_filesystem(str(library.id), fname, folder_id=str(folder["id"]),
                                                      link_data_only="link_to_files",
                                                      dbkey=sample_info["genome_build"],
                                                      roles=str(library.roles))

def _get_folder(gi, folder_name, library, libitems):
    """Retrieve or create a folder inside the library with the right now.
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
        return _get_library_from_name(gi, galaxy_lib, role, sample_info)
    elif config.get("private_libs") or config.get("lab_association") or config.get("researcher"):
        return _library_from_nglims(gi, sample_info, config)
    else:
        raise ValueError("No Galaxy library specified for sample: %s" %
                         sample_info["description"])

def _get_library_from_name(gi, name, role, sample_info, create=False):
    for lib in gi.libraries.get_libraries():
        if lib["name"].lower().find(name.lower()) >= 0:
            return GalaxyLibrary(lib["id"], lib["name"], role)
    if create and name:
        logger.info("Creating Galaxy library: '%s'" % name)
        lib = gi.libraries.create_library(name)[0]
        return GalaxyLibrary(lib["id"], lib["name"], role)
    else:
        raise ValueError("Could not find Galaxy library matching '%s' for sample %s" %
                         (name, sample_info["description"]))

def _library_from_nglims(gi, sample_info, config):
    """Retrieve upload library from nglims specified user libraries.
    """
    names = [config.get(x, "").strip() for x in ["lab_association", "researcher"]
             if config.get(x)]
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
