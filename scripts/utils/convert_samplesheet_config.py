#!/usr/bin/env python
"""Convert Illumina SampleSheet CSV files to the run_info.yaml input file.

This allows running the analysis pipeline without Galaxy, using CSV input
files from Illumina SampleSheet or Genesifter.

Usage:
  convert_samplesheet_config.py <input csv>
"""

import sys
import os
from bcbio.solexa.samplesheet import SampleSheet

if __name__ == "__main__":
    
    in_file = sys.argv[1]
    out_file = "%s.yaml" % os.path.splitext(in_file)[0]
    
    sh = SampleSheet(in_file)
    yaml_sheet = sh.csv2yaml(in_file)
    
    open(out_file, "w").write(yaml_sheet)
