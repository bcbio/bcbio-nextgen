"""Support integration with Illumina sequencer machines.
"""
from __future__ import print_function
import glob
import json
import os
import operator
import subprocess
from xml.etree.ElementTree import ElementTree

import requests
import yaml

from bcbio import utils
from bcbio.log import setup_local_logging
from bcbio.illumina import demultiplex, samplesheet, transfer
from bcbio.galaxy import nglims
from functools import reduce

# ## bcbio-nextgen integration

def check_and_postprocess(args):
    """Check for newly dumped sequencer output, post-processing and transferring.
    """
    with open(args.process_config) as in_handle:
        config = yaml.safe_load(in_handle)
    setup_local_logging(config)
    for dname in _find_unprocessed(config):
        lane_details = nglims.get_runinfo(config["galaxy_url"], config["galaxy_apikey"], dname,
                                          utils.get_in(config, ("process", "storedir")))
        if isinstance(lane_details, dict) and "error" in lane_details:
            print("Flowcell not found in Galaxy: %s" % lane_details)
        else:
            lane_details = _tweak_lane(lane_details, dname)
            fcid_ss = samplesheet.from_flowcell(dname, lane_details)
            _update_reported(config["msg_db"], dname)
            fastq_dir = demultiplex.run_bcl2fastq(dname, fcid_ss, config)
            bcbio_config, ready_fastq_dir = nglims.prep_samples_and_config(dname, lane_details, fastq_dir, config)
            transfer.copy_flowcell(dname, ready_fastq_dir, bcbio_config, config)
            _start_processing(dname, bcbio_config, config)

def _tweak_lane(lane_details, dname):
    """Potentially tweak lane information to handle custom processing, reading a lane_config.yaml file.
    """
    tweak_config_file = os.path.join(dname, "lane_config.yaml")
    if os.path.exists(tweak_config_file):
        with open(tweak_config_file) as in_handle:
            tweak_config = yaml.safe_load(in_handle)
        if tweak_config.get("uniquify_lanes"):
            out = []
            for ld in lane_details:
                ld["name"] = "%s-%s" % (ld["name"], ld["lane"])
                out.append(ld)
            return out
    return lane_details

def _remap_dirname(local, remote):
    """Remap directory names from local to remote.
    """
    def do(x):
        return x.replace(local, remote, 1)
    return do

def _start_processing(dname, sample_file, config):
    """Initiate processing: on a remote server or locally on a cluster.
    """
    to_remote = _remap_dirname(dname, os.path.join(utils.get_in(config, ("process", "dir")),
                                                   os.path.basename(dname)))
    args = {"work_dir": to_remote(os.path.join(dname, "analysis")),
            "run_config": to_remote(sample_file),
            "fc_dir": to_remote(dname)}
    # call a remote server
    if utils.get_in(config, ("process", "server")):
        print("%s/run?args=%s" % (utils.get_in(config, ("process", "server")), json.dumps(args)))
        requests.get(url="%s/run" % utils.get_in(config, ("process", "server")),
                     params={"args": json.dumps(args)})
    # submit to a cluster scheduler
    elif "submit_cmd" in config["process"] and "bcbio_batch" in config["process"]:
        with utils.chdir(utils.safe_makedir(args["work_dir"])):
            batch_script = "submit_bcbio.sh"
            with open(batch_script, "w") as out_handle:
                out_handle.write(config["process"]["bcbio_batch"].format(fcdir=args["fc_dir"],
                                                                         run_config=args["run_config"]))
            submit_cmd = utils.get_in(config, ("process", "submit_cmd"))
            subprocess.check_call(submit_cmd.format(batch_script=batch_script), shell=True)
    else:
        raise ValueError("Unexpected processing approach: %s" % config["process"])

def add_subparser(subparsers):
    """Add command line arguments for post-processing sequencer results.
    """
    parser = subparsers.add_parser("sequencer", help="Post process results from a sequencer.")
    parser.add_argument("process_config", help="YAML file specifying sequencer details for post-processing.")
    return parser

# ## Dump directory processing

def _find_unprocessed(config):
    """Find any finished directories that have not been processed.
    """
    reported = _read_reported(config["msg_db"])
    for dname in _get_directories(config):
        if os.path.isdir(dname) and dname not in reported:
            if _is_finished_dumping(dname):
                yield dname

def _get_directories(config):
    for directory in config["dump_directories"]:
        for dname in sorted(glob.glob(os.path.join(directory, "*[Aa]*[Xx][XxYy23]"))):
            if os.path.isdir(dname):
                yield dname

def _is_finished_dumping(directory):
    """Determine if the sequencing directory has all files.

    The final checkpoint file will differ depending if we are a
    single or paired end run.
    """
    #if _is_finished_dumping_checkpoint(directory):
    #    return True
    # Check final output files; handles both HiSeq and GAII
    run_info = os.path.join(directory, "RunInfo.xml")
    hi_seq_checkpoint = "Basecalling_Netcopy_complete_Read%s.txt" % \
                        _expected_reads(run_info)
    to_check = ["Basecalling_Netcopy_complete_SINGLEREAD.txt",
                "Basecalling_Netcopy_complete_READ2.txt",
                hi_seq_checkpoint]
    return reduce(operator.or_,
                  [os.path.exists(os.path.join(directory, f)) for f in to_check])

def _is_finished_dumping_checkpoint(directory):
    """Recent versions of RTA (1.10 or better), write the complete file.

    This is the most straightforward source but as of 1.10 still does not
    work correctly as the file will be created at the end of Read 1 even
    if there are multiple reads.
    """
    check_file = os.path.join(directory, "Basecalling_Netcopy_complete.txt")
    check_v1, check_v2 = (1, 10)
    if os.path.exists(check_file):
        with open(check_file) as in_handle:
            line = in_handle.readline().strip()
        if line:
            version = line.split()[-1]
            v1, v2 = [float(v) for v in version.split(".")[:2]]
            if ((v1 > check_v1) or (v1 == check_v1 and v2 >= check_v2)):
                return True

def _expected_reads(run_info_file):
    """Parse the number of expected reads from the RunInfo.xml file.
    """
    reads = []
    if os.path.exists(run_info_file):
        tree = ElementTree()
        tree.parse(run_info_file)
        read_elem = tree.find("Run/Reads")
        reads = read_elem.findall("Read")
    return len(reads)

# ## Flat file of processed directories

def _read_reported(msg_db):
    """Retrieve a list of directories previous reported.
    """
    reported = []
    if os.path.exists(msg_db):
        with open(msg_db) as in_handle:
            for line in in_handle:
                reported.append(line.strip())
    return reported

def _update_reported(msg_db, new_dname):
    """Add a new directory to the database of reported messages.
    """
    with open(msg_db, "a") as out_handle:
        out_handle.write("%s\n" % new_dname)
