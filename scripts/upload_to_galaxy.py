#!/usr/bin/env python
"""Upload a set of next-gen sequencing data files to a data library in Galaxy.

Usage:
    upload_to_galaxy.py <config file> <flowcell directory> <analysis output dir>
                        [<YAML run information>]

The optional <YAML run information> file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/run_info.yaml'

The configuration file is in YAML format with the following key/value pairs:

galaxy_url: Base URL of Galaxy for uploading.
galaxy_api_key: Developer's API key.
galaxy_config: Path to Galaxy's universe_wsgi.ini file. This is required so
we know where to organize directories for upload based on library_import_dir.
"""
import sys
import os
import glob
import shutil
import ConfigParser
import time

import yaml

from bcbio.solexa.flowcell import get_fastq_dir
from bcbio.pipeline.run_info import get_run_info
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio import utils

def main(config_file, fc_dir, analysis_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    galaxy_api = (GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
                  if config.has_key("galaxy_api_key") else None)
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)

    base_folder_name = "%s_%s" % (fc_date, fc_name)
    run_details = lims_run_details(run_info, base_folder_name)
    for (library_name, access_role, dbkey, lane, bc_id, name, desc,
            local_name) in run_details:
        library_id = (get_galaxy_library(library_name, galaxy_api)
                      if library_name else None)
        upload_files = list(select_upload_files(local_name, bc_id, fc_dir,
                                                analysis_dir, config))
        if len(upload_files) > 0:
            print lane, bc_id, name, desc, library_name
            print "Creating storage directory"
            if library_id:
                folder, cur_galaxy_files = get_galaxy_folder(library_id,
                               base_folder_name, name, desc, galaxy_api)
            else:
                cur_galaxy_files = []
            store_dir = move_to_storage(lane, bc_id, base_folder_name, upload_files,
                                        cur_galaxy_files, config, config_file)
            if store_dir and library_id:
                print "Uploading directory of files to Galaxy"
                print galaxy_api.upload_directory(library_id, folder['id'],
                                                  store_dir, dbkey, access_role)
    if galaxy_api and not run_info_yaml:
        add_run_summary_metrics(analysis_dir, galaxy_api)

# LIMS specific code for retrieving information on what to upload from
# the Galaxy NGLIMs.
# Also includes function for selecting files to upload from flow cell and
# analysis directories.
# These should be edited to match a local workflow if adjusting this.

def lims_run_details(run_info, base_folder_name):
    """Retrieve run infomation on a flow cell from Next Gen LIMS.
    """
    for lane_items in run_info["details"]:
        for lane_info in lane_items:
            if not run_info["run_id"] or lane_info.has_key("researcher"):
                if lane_info.get("private_libs", None) is not None:
                    libname, role = _get_galaxy_libname(lane_info["private_libs"],
                                                        lane_info["lab_association"],
                                                        lane_info["researcher"])
                elif lane_info.has_key("galaxy_library"):
                    libname = lane_info["galaxy_library"]
                    role = lane_info["galaxy_role"]
                else:
                    libname, role = (None, None)
                remote_folder = str(lane_info.get("name", lane_info["lane"]))
                description = ": ".join([lane_info[n] for n in ["researcher", "description"]
                                         if lane_info.has_key(n)])
                local_name = "%s_%s" % (lane_info["lane"], base_folder_name)
                if lane_info["barcode_id"] is not None:
                    remote_folder += "_%s" % lane_info["barcode_id"]
                    local_name += "_%s" % lane_info["barcode_id"]
                yield (libname, role, lane_info["genome_build"],
                       lane_info["lane"], lane_info["barcode_id"],
                       remote_folder, description, local_name)

def _get_galaxy_libname(private_libs, lab_association, researcher):
    """Retrieve most appropriate Galaxy data library.

    Gives preference to private data libraries associated with the user. If not
    found will create a user specific library.
    """
    print private_libs, lab_association
    # simple case -- one private library. Upload there
    if len(private_libs) == 1:
        return private_libs[0]
    # no private libraries -- use the lab association or researcher name
    elif len(private_libs) == 0:
        if not lab_association:
            return researcher, ""
        else:
            return lab_association, ""
    # multiple libraries -- find the one that matches the lab association
    else:
        check_libs = [l.lower() for (l, _) in private_libs]
        try:
            i = check_libs.index(lab_association.lower())
            return private_libs[i]
        # can't find the lab association, give us the first library
        except (IndexError, ValueError):
            return private_libs[0]

def select_upload_files(base, bc_id, fc_dir, analysis_dir, config):
    """Select fastq, bam alignment and summary files for upload to Galaxy.
    """
    base_glob = _dir_glob(base, analysis_dir)
    # Configurable upload of fastq files -- BAM provide same information, compacted
    if config["algorithm"].get("upload_fastq", True):
        # look for fastq files in a barcode directory or the main fastq directory
        bc_base = base.rsplit("_", 1)[0] if bc_id else base
        bc_dir = os.path.join(analysis_dir, "%s_barcode" % bc_base)
        fastq_glob = "%s_*fastq.txt" % base
        found_fastq = False
        for fname in glob.glob(os.path.join(bc_dir, fastq_glob)):
            found_fastq = True
            yield (fname, os.path.basename(fname))
        if not found_fastq:
            fastq_dir = get_fastq_dir(fc_dir)
            for fname in glob.glob(os.path.join(fastq_dir, fastq_glob)):
                yield (fname, os.path.basename(fname))
    for summary_file in base_glob("summary.pdf"):
        yield (summary_file, _name_with_ext(summary_file, "-summary.pdf"))
    for wig_file in base_glob(".bigwig"):
        yield (wig_file, _name_with_ext(wig_file, "-coverage.bigwig"))
    # upload BAM files, preferring recalibrated and realigned files
    found_recal = False
    for bam_file in base_glob("gatkrecal-realign-dup.bam"):
        found_recal = True
        yield (bam_file, _name_with_ext(bam_file, "-gatkrecal-realign.bam"))
    if not found_recal:
        for bam_file in base_glob("gatkrecal.bam"):
            found_recal = True
            yield (bam_file, _name_with_ext(bam_file, "-gatkrecal.bam"))
    if not found_recal:
        for bam_file in base_glob("sort-dup.bam"):
            yield (bam_file, _name_with_ext(bam_file, ".bam"))
    # Genotype files produced by SNP calling
    for snp_file in base_glob("variants-combined-annotated.vcf"):
        yield (snp_file, _name_with_ext(bam_file, "-variants.vcf"))
    # Effect information on SNPs
    for snp_file in base_glob("variants-combined-effects.tsv"):
        yield (snp_file, _name_with_ext(bam_file, "-variants-effects.tsv"))

def _dir_glob(base, work_dir):
    # Allowed characters that can trail the base. This prevents picking up
    # NAME_10 when globbing for NAME_1
    trailers = "[-_.]"
    def _safe_glob(ext):
        return glob.glob(os.path.join(work_dir, "%s%s*%s" % (base, trailers, ext)))
    return _safe_glob

def _name_with_ext(orig_file, ext):
    """Return a normalized filename without internal processing names.
    """
    base = os.path.basename(orig_file).split("-")[0]
    for extra in ["_trim"]:
        if base.endswith(extra):
            base = base[:-len(extra)]
    return "%s%s" % (base, ext)

def add_run_summary_metrics(analysis_dir, galaxy_api):
    """Upload YAML file of run information to Galaxy though the NGLims API.
    """
    run_file = os.path.join(analysis_dir, "run_summary.yaml")
    if os.path.exists(run_file):
        with open(run_file) as in_handle:
            run_summary = yaml.load(in_handle)
        galaxy_api.sqn_run_summary(run_summary)

# General functionality for interacting with Galaxy via the Library API

def get_galaxy_folder(library_id, folder_name, lane, description, galaxy_api):
    """Return or create a folder within the given library.

    Creates or retrieves a top level directory for a run, and then creates
    a lane specific directory within this run.
    """
    items = galaxy_api.library_contents(library_id)
    root = _folders_by_name('/', items)[0]
    run_folders = _safe_get_folders("/%s" % folder_name, items,
                                    library_id, root["id"], folder_name, "",
                                    galaxy_api)
    lane_folders = _safe_get_folders("/%s/%s" % (folder_name, lane), items,
                                     library_id, run_folders[0]['id'],
                                     str(lane), description, galaxy_api)
    cur_files = [f for f in items if f['type'] == 'file'
                 and f['name'].startswith("/%s/%s" % (folder_name, lane))]
    return lane_folders[0], cur_files

def _safe_get_folders(base_name, items, library_id, base_folder_id, name,
                      description, galaxy_api):
    """Retrieve folders for a run or lane, retrying in the case of network errors.
    """
    max_tries = 5
    num_tries = 0
    while 1:
        try:
            folders = _folders_by_name(base_name, items)
            if len(folders) == 0:
                folders = galaxy_api.create_folder(library_id, base_folder_id,
                                                   name, description)
            break
        except ValueError:
            if num_tries > max_tries:
                raise
            time.sleep(2)
            num_tries += 1
    return folders

def _folders_by_name(name, items):
    return [f for f in items if f['type'] == 'folder' and
                                f['name'] == name]

def move_to_storage(lane, bc_id, fc_dir, select_files, cur_galaxy_files,
                    config, config_file):
    """Create directory for long term storage before linking to Galaxy.
    """
    galaxy_config_file = utils.add_full_path(config["galaxy_config"],
                                             os.path.dirname(config_file))
    galaxy_conf = ConfigParser.SafeConfigParser({'here' : ''})
    galaxy_conf.read(galaxy_config_file)
    try:
        lib_import_dir = galaxy_conf.get("app:main", "library_import_dir")
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        raise ValueError("Galaxy config %s needs library_import_dir to be set."
                         % galaxy_config_file)
    storage_dir = _get_storage_dir(fc_dir, lane, bc_id, os.path.join(lib_import_dir,
                                   "storage"))
    existing_files = [os.path.basename(f['name']) for f in cur_galaxy_files]
    need_upload = False
    for orig_file, new_file in select_files:
        if new_file in existing_files:
            need_upload = False
            break
        else:
            new_file = os.path.join(storage_dir, new_file)
            if not os.path.exists(new_file):
                shutil.copy(orig_file, new_file)
            need_upload = True
    return (storage_dir if need_upload else None)

def _get_storage_dir(cur_folder, lane, bc_id, storage_base):
    base = "%s_%s" % (lane, bc_id) if bc_id else str(lane)
    store_dir = os.path.join(storage_base, cur_folder, base)
    utils.safe_makedir(store_dir)
    return store_dir

def get_galaxy_library(lab_association, galaxy_api):
    ret_info = None
    for lib_info in galaxy_api.get_libraries():
        if lib_info["name"].find(lab_association) >= 0:
            ret_info = lib_info
            break
    # need to add a new library
    if ret_info is None:
        ret_info = galaxy_api.create_library(lab_association)[0]
    return ret_info["id"]

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    main(*sys.argv[1:])
