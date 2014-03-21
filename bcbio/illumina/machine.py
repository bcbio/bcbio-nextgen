"""Support integration with Illumina sequencer machines.
"""
import glob
import os
import operator
from xml.etree.ElementTree import ElementTree

import yaml
import logbook

from bcbio.log import setup_local_logging, logger
from bcbio.illumina import demultiplex, samplesheet
from bcbio.galaxy import nglims

# ## bcbio-nextgen integration

def check_and_postprocess(args):
    """Check for newly dumped sequencer output, post-processing and transferring.
    """
    with open(args.process_config) as in_handle:
        config = yaml.safe_load(in_handle)
    setup_local_logging()
    for dname in _find_unprocessed(config):
        runinfo = nglims.get_runinfo(config["galaxy_url"], config["galaxy_apikey"], dname)
        lane_details = nglims.flatten_lane_detail(runinfo)
        fcid_ss = samplesheet.from_flowcell(dname, lane_details)
        fastq_dir = demultiplex.run_bcl2fastq(dname, fcid_ss, config)
        #_update_reported(config["msg_db"], dname)

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
        for dname in sorted(glob.glob(os.path.join(directory, "*[Aa]*[Xx][Xx]"))):
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
