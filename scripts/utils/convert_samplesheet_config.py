#!/usr/bin/env python
"""Convert Illumina SampleSheet CSV files to the run_info.yaml input file.

This allows running the analysis pipeline without Galaxy, using CSV input
files from Illumina SampleSheet or Genesifter.

Usage:
  convert_samplesheet_config.py <input csv>
"""
import os
import sys
import csv
import itertools

import yaml

def main(in_file):
    out_file = "%s.yaml" % os.path.splitext(in_file)[0]
    barcode_ids = generate_barcode_ids(read_input_csv(in_file))
    lanes = organize_lanes(read_input_csv(in_file), barcode_ids)
    with open(out_file, "w") as out_handle:
        yaml.dump(lanes, out_handle, default_flow_style=False)

def organize_lanes(info_iter, barcode_ids):
    """Organize flat lane information into nested YAML structure.
    """
    all_lanes = []
    for (lane, org), info in itertools.groupby(info_iter, lambda x: (x[0], x[2])):
        cur_lane = dict(lane=lane, genome_build=org, analysis="Standard")
        info = list(info)
        if len(info) == 1: # non-barcoded sample
            cur_lane["description"] = info[0][1]
        else: # barcoded sample
            cur_lane["description"] = "Barcoded %s" % lane
            multiplex = []
            for (_, sample_id, _, bc_seq) in info:
                bc_type, bc_id = barcode_ids[bc_seq]
                multiplex.append(dict(barcode_type=bc_type,
                                      barcode_id=bc_id,
                                      sequence=bc_seq,
                                      name=sample_id))
            cur_lane["multiplex"] = multiplex
        all_lanes.append(cur_lane)
    return all_lanes

def generate_barcode_ids(info_iter):
    """Create unique barcode IDs assigned to sequences
    """
    bc_type = "SampleSheet"
    barcodes = list(set([x[-1] for x in info_iter]))
    barcodes.sort()
    barcode_ids = {}
    for i, bc in enumerate(barcodes):
        barcode_ids[bc] = (bc_type, i+1)
    return barcode_ids

def read_input_csv(in_file):
    """Parse useful details from SampleSheet CSV file.
    """
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        reader.next() # header
        for line in reader:
            (_, lane, sample_id, genome, barcode) = line[:5]
            yield lane, sample_id, genome, barcode

if __name__ == "__main__":
    main(*sys.argv[1:])
