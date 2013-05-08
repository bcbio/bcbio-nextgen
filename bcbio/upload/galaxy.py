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
    folder_name = "%s_%s" % (config["fc_name"], config["fc_date"])
    storage_dir = utils.safe_makedir(os.path.join(config["dir"], folder_name))
    storage_file = filesystem.copy_finfo(finfo, storage_dir)
    if config.has_key("galaxy_url") and config.has_key("galaxy_api_key"):
        gi = GalaxyInstance(config["galaxy_url"], config["galaxy_api_key"])
    else:
        raise ValueError("Galaxy upload requires `galaxy_url` and `galaxy_api_key` in config")
    if storage_file and sample_info:
        _to_datalibrary(storage_file, gi, folder_name, sample_info, config)

def _to_datalibrary(fname, gi, folder_name, sample_info, config):
    """Upload a file to a Galaxy data library in a project specific folder.
    """
    library = _get_library(gi, sample_info, config)
    libitems =  gi.libraries.show_library(library.id, contents=True)
    folder = _get_folder(gi, folder_name, library, libitems)
    _file_to_folder(gi, fname, sample_info, libitems, library, folder)

def _file_to_folder(gi, fname, sample_info, libitems, library, folder):
    """Check if file exists on Galaxy, if not upload to specified folder.
    """
    full_name = os.path.join(folder["name"], os.path.basename(fname))
    for item in libitems:
        if item["name"] == full_name:
            return item
    logger.info("Uploading to Galaxy: %s" % full_name)
    return gi.libraries.upload_from_galaxy_filesystem(library.id, fname, folder_id=folder["id"],
                                                      link_data_only="link_to_files",
                                                      dbkey=sample_info["genome_build"],
                                                      roles=library.roles)

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
    elif sample_info.get("private_libs"):
        return _library_from_nglims(gi, sample_info)
    else:
        raise ValueError("No Galaxy library specified for sample: %s" %
                         sample_info["description"])

def _get_library_from_name(gi, name, role, sample_info):
    for lib in gi.libraries.get_libraries():
        if lib["name"].lower().find(name.lower()) >= 0:
            return GalaxyLibrary(lib["id"], lib["name"], role)
    else:
        raise ValueError("Could not find Galaxy library matching %s for sample %s" %
                         (name, sample_info["description"]))

def _library_from_nglims(gi, sample_info):
    """Retrieve upload library from nglims specified user libraries.
    """
    check_names = [sample_info[x].lower()
                   for x in ["lab_association", "researcher"]
                   if sample_info[x]]
    for libname, role in sample_info["private_libs"]:
        # Try to find library for lab or rsearcher
        if libname.lower() in check_names:
            return _get_library_from_name(gi, libname, role, sample_info)
    # default to first private library if available
    if len(sample_info["private_libs"]) > 0:
        libname, role = sample_info["private_libs"][0]
        return _get_library_from_name(gi, libname, role, sample_info)
    # otherwise use the lab association or researcher name
    else:
        libname = check_names[0]
        return _get_library_from_name(gi, libname, None, sample_info)
