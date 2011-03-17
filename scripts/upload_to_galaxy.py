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
import urllib
import urllib2
import json

import yaml

from bcbio.solexa.flowcell import get_flowcell_info, get_fastq_dir
from bcbio.galaxy.api import GalaxyApiAccess

def main(config_file, fc_dir, analysis_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    fc_name, fc_date = get_flowcell_info(fc_dir)
    if run_info_yaml:
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        run_info = dict(details=run_details, run_id="")
        galaxy_api = None
    else:
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        run_info = galaxy_api.run_details(fc_name)

    base_folder_name = "%s_%s" % (fc_date, fc_name)
    run_details = lims_run_details(run_info, fc_name, base_folder_name)
    for (library_name, access_role, dbkey, lane, bc_id, name, desc,
            local_name) in run_details:
        library_id = (get_galaxy_library(library_name, galaxy_api)
                      if library_name else None)
        upload_files = list(select_upload_files(local_name, bc_id, fc_dir,
            analysis_dir))
        if len(upload_files) > 0:
            print lane, bc_id, name, desc, library_name
            print "Creating storage directory"
            if library_id:
                folder, cur_galaxy_files = get_galaxy_folder(library_id,
                               base_folder_name, name, desc, galaxy_api)
            else:
                cur_galaxy_files = []
            store_dir = move_to_storage(lane, bc_id, base_folder_name, upload_files,
                    cur_galaxy_files, config)
            if store_dir and library_id:
                print "Uploading directory of files to Galaxy"
                print galaxy_api.upload_directory(library_id, folder['id'],
                                                  store_dir, dbkey, access_role)
    if galaxy_api:
        add_run_summary_metrics(analysis_dir, galaxy_api)

# LIMS specific code for retrieving information on what to upload from
# the Galaxy NGLIMs.
# Also includes function for selecting files to upload from flow cell and
# analysis directories.
# These should be editing to match a local workflow if adjusting this.

def lims_run_details(run_info, fc_name, base_folder_name):
    """Retrieve run infomation on a flow cell from Next Gen LIMS.
    """
    for lane_info in (l for l in run_info["details"] if l.has_key("researcher")
                      or not run_info["run_id"]):
        if lane_info.get("private_libs", None) is not None:
            libname, role = _get_galaxy_libname(lane_info["private_libs"],
                                                lane_info["lab_association"],
                                                lane_info["researcher"])
        else:
            libname, role = (None, None)
        for barcode in lane_info.get("multiplex", [None]):
            remote_folder = lane_info.get("name", "")
            description = "%s: %s" % (lane_info.get("researcher", ""),
                    lane_info["description"])
            local_name = "%s_%s" % (lane_info["lane"], base_folder_name)
            if barcode:
                remote_folder += "_%s" % barcode["barcode_id"]
                description += ": %s" % barcode["name"]
                local_name += "_%s" % barcode["barcode_id"]
            yield (libname, role, lane_info["genome_build"],
                    lane_info["lane"], barcode["barcode_id"] if barcode else "",
                    remote_folder, description, local_name)

def _get_galaxy_libname(private_libs, lab_association, researcher):
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
        except IndexError:
            return private_libs[0]

def select_upload_files(base, bc_id, fc_dir, analysis_dir):
    """Select fastq, bam alignment and summary files for upload to Galaxy.
    """
    fastq_dir = analysis_dir if bc_id else get_fastq_dir(fc_dir)
    for fname in glob.glob(os.path.join(fastq_dir, "%s_*fastq.txt" % base)):
        yield (fname, os.path.basename(fname))
    for summary_file in glob.glob(os.path.join(analysis_dir,
            "%s-*summary.pdf" % base)):
        yield (summary_file, _name_with_ext(summary_file, "-summary.pdf"))
    for bam_file in glob.glob(os.path.join(analysis_dir,
            "%s-*sort-dup.bam" % base)):
        yield (bam_file, _name_with_ext(bam_file, ".bam"))
    for wig_file in glob.glob(os.path.join(analysis_dir,
            "%s-*sort.bigwig" % base)):
        yield (wig_file, _name_with_ext(wig_file, "-coverage.bigwig"))
    # upload any recalibrated BAM files used for SNP calling
    found_recal = False
    for bam_file in glob.glob(os.path.join(analysis_dir,
            "%s-*gatkrecal-realign-sort.bam" % base)):
        found_recal = True
        yield (bam_file, _name_with_ext(bam_file, "-gatkrecal-realign.bam"))
    if not found_recal:
        for bam_file in glob.glob(os.path.join(analysis_dir,
                "%s-*gatkrecal.bam" % base)):
            yield (bam_file, _name_with_ext(bam_file, "-gatkrecal.bam"))
    # Genotype files produced by SNP calling
    for snp_file in glob.glob(os.path.join(analysis_dir,
            "%s-*snp-filter.vcf" % base)):
        yield (snp_file, _name_with_ext(bam_file, "-snp-filter.vcf"))
    # Effect information on SNPs
    for snp_file in glob.glob(os.path.join(analysis_dir,
            "%s-*snp-filter-effects.tsv" % base)):
        yield (snp_file, _name_with_ext(bam_file, "-snp-effects.tsv"))

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
    run_folders = _folders_by_name("/%s" % folder_name, items)
    if len(run_folders) == 0:
        run_folders = galaxy_api.create_folder(library_id,
                root['id'], folder_name)
    lane_folders = _folders_by_name("/%s/%s" % (folder_name, lane), items)
    if len(lane_folders) == 0:
        lane_folders = galaxy_api.create_folder(library_id,
                run_folders[0]['id'], str(lane), description)
        cur_files = []
    else:
        cur_files = [f for f in items if f['type'] == 'file'
                and f['name'].startswith("/%s/%s" % (folder_name, lane))]
    return lane_folders[0], cur_files

def _folders_by_name(name, items):
    return [f for f in items if f['type'] == 'folder' and
                                f['name'] == name]

def move_to_storage(lane, bc_id, fc_dir, select_files, cur_galaxy_files, config):
    """Create directory for long term storage before linking to Galaxy.
    """
    galaxy_conf = ConfigParser.SafeConfigParser({'here' : ''})
    galaxy_conf.read(config["galaxy_config"])
    try:
        lib_import_dir = galaxy_conf.get("app:main", "library_import_dir")
    except ConfigParser.NoOptionError:
        raise ValueError("Galaxy config %s needs library_import_dir to be set."
                % config["galaxy_config"])
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
    if not os.path.exists(store_dir):
        os.makedirs(store_dir)
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
