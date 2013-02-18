"""Move files to local Galaxy upload directory and add to Galaxy Data Libraries.

Required configurable variables in upload:
  dir
  galaxy_url
  galaxy_api_key
  galaxy_library
"""
import collections
import os

from bioblend.galaxy import GalaxyInstance

from bcbio import utils
from bcbio.upload import filesystem

def update_file(finfo, sample_info, config):
    """Update file in Galaxy data libraries.
    """
    folder_name = "%s_%s" % (config["fc_name"], config["fc_date"])
    storage_dir = utils.safe_makedir(os.path.join(config["dir"], folder_name))
    storage_file = filesystem.copy_finfo(finfo, storage_dir)
    if config.has_key("galaxy_url") and config.has_key("galaxy_api_key"):
        gi = GalaxyInstance(config["galaxy_url"], config["galaxy_api_key"])
    else:
        raise ValueError("Galaxy upload requires `galaxy_url` and `galaxy_api_key` in config")
    if storage_file:
        _to_datalibrary(storage_file, gi, folder_name, sample_info, config)

def _to_datalibrary(fname, gi, folder_name, sample_info, config):
    library = _get_library(gi, sample_info, config)
    folder = _get_folder(gi, folder_name, library)
    gi.libraries.upload_from_galaxy_filesystem(library.id, fname, folder_id=folder["id"],
                                               link_data_only="link_to_files",
                                               dbkey=sample_info["genome_build"])

def _get_folder(gi, folder_name, library):
    """Retrieve or create a folder inside the library with the right now.
    """
    for item in gi.libraries.show_library(library.id, contents=True):
        if item["type"] == "folder" and item["name"] == "/%s" % folder_name:
            return item
    return gi.libraries.create_folder(library.id, folder_name)[0]

def _get_library(gi, sample_info, config):
    """Retrieve the appropriate data library for the current user.
    """
    Library = collections.namedtuple("Library", ["id", "name", "role"])
    galaxy_lib = sample_info.get("galaxy_library",
                                 config.get("galaxy_library"))
    if galaxy_lib:
        for lib in gi.libraries.get_libraries():
            if lib["name"].lower().find(galaxy_lib.lower()) >= 0:
                return Library(lib["id"], lib["name"], "")
        else:
            raise ValueError("Could not find Galaxy library matching %s for sample %s" %
                             (galaxy_lib, sample_info["description"]))
    else:
        raise ValueError("No Galaxy library specified for sample: %s" %
                         sample_info["description"])
