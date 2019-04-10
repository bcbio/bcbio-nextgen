"""Converts Illumina SampleSheet CSV files to the run_info.yaml input file.

This allows running the analysis pipeline without Galaxy, using CSV input
files from Illumina SampleSheet or Genesifter.
"""
import os
import csv
import itertools
import difflib
import glob
import io

import yaml

from bcbio.illumina import flowcell
from bcbio import utils

# ## Create samplesheets

def from_flowcell(run_folder, lane_details, out_dir=None):
    """Convert a flowcell into a samplesheet for demultiplexing.
    """
    fcid = os.path.basename(run_folder)
    if out_dir is None:
        out_dir = run_folder
    out_file = os.path.join(out_dir, "%s.csv" % fcid)
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["FCID", "Lane", "Sample_ID", "SampleRef", "Index",
                         "Description", "Control", "Recipe", "Operator", "SampleProject"])
        for ldetail in lane_details:
            writer.writerow(_lane_detail_to_ss(fcid, ldetail))
    return out_file

def _lane_detail_to_ss(fcid, ldetail):
    """Convert information about a lane into Illumina samplesheet output.
    """
    return [fcid, ldetail["lane"], ldetail["name"], ldetail["genome_build"],
            ldetail["bc_index"], ldetail["description"].encode("ascii", "ignore"), "N", "", "",
            ldetail["project_name"]]

# ## Use samplesheets to create YAML files

def _organize_lanes(info_iter, barcode_ids):
    """Organize flat lane information into nested YAML structure.
    """
    all_lanes = []
    for (fcid, lane, sampleref), info in itertools.groupby(info_iter, lambda x: (x[0], x[1], x[1])):
        info = list(info)
        cur_lane = dict(flowcell_id=fcid, lane=lane, genome_build=info[0][3], analysis="Standard")
        if not _has_barcode(info):
            cur_lane["description"] = info[0][1]
        else: # barcoded sample
            cur_lane["description"] = "Barcoded lane %s" % lane
            multiplex = []
            for (_, _, sample_id, _, bc_seq) in info:
                bc_type, bc_id = barcode_ids[bc_seq]
                multiplex.append(dict(barcode_type=bc_type,
                                      barcode_id=bc_id,
                                      sequence=bc_seq,
                                      name=sample_id))
            cur_lane["multiplex"] = multiplex
        all_lanes.append(cur_lane)
    return all_lanes

def _has_barcode(sample):
    if sample[0][4]:
        return True

def _generate_barcode_ids(info_iter):
    """Create unique barcode IDs assigned to sequences
    """
    bc_type = "SampleSheet"
    barcodes = list(set([x[-1] for x in info_iter]))
    barcodes.sort()
    barcode_ids = {}
    for i, bc in enumerate(barcodes):
        barcode_ids[bc] = (bc_type, i+1)
    return barcode_ids

def _read_input_csv(in_file):
    """Parse useful details from SampleSheet CSV file.
    """
    with io.open(in_file, newline=None) as in_handle:
        reader = csv.reader(in_handle)
        next(reader) # header
        for line in reader:
            if line: # empty lines
                (fc_id, lane, sample_id, genome, barcode) = line[:5]
                yield fc_id, lane, sample_id, genome, barcode

def _get_flowcell_id(in_file, require_single=True):
    """Retrieve the unique flowcell id represented in the SampleSheet.
    """
    fc_ids = set([x[0] for x in _read_input_csv(in_file)])
    if require_single and len(fc_ids) > 1:
        raise ValueError("There are several FCIDs in the same samplesheet file: %s" % in_file)
    else:
        return fc_ids

def csv2yaml(in_file, out_file=None):
    """Convert a CSV SampleSheet to YAML run_info format.
    """
    if out_file is None:
        out_file = "%s.yaml" % os.path.splitext(in_file)[0]
    barcode_ids = _generate_barcode_ids(_read_input_csv(in_file))
    lanes = _organize_lanes(_read_input_csv(in_file), barcode_ids)
    with open(out_file, "w") as out_handle:
        out_handle.write(yaml.safe_dump(lanes, default_flow_style=False))
    return out_file

def run_has_samplesheet(fc_dir, config, require_single=True):
    """Checks if there's a suitable SampleSheet.csv present for the run
    """
    fc_name, _ = flowcell.parse_dirname(fc_dir)
    sheet_dirs = config.get("samplesheet_directories", [])
    fcid_sheet = {}
    for ss_dir in (s for s in sheet_dirs if os.path.exists(s)):
        with utils.chdir(ss_dir):
            for ss in glob.glob("*.csv"):
                fc_ids = _get_flowcell_id(ss, require_single)
                for fcid in fc_ids:
                    if fcid:
                        fcid_sheet[fcid] = os.path.join(ss_dir, ss)
    # difflib handles human errors while entering data on the SampleSheet.
    # Only one best candidate is returned (if any). 0.85 cutoff allows for
    # maximum of 2 mismatches in fcid

    potential_fcids = difflib.get_close_matches(fc_name, fcid_sheet.keys(), 1, 0.85)
    if len(potential_fcids) > 0 and potential_fcids[0] in fcid_sheet:
        return fcid_sheet[potential_fcids[0]]
    else:
        return None
