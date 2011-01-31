#!/usr/bin/env python
"""Convert Illumina SampleSheet CSV files to the run_info.yaml input file.

This allows running the analysis pipeline without Galaxy, using CSV input
files from Illumina SampleSheet or Genesifter.

Usage:
  convert_samplesheet_config.py <input csv>
"""
import sys

from bcbio.solexa import samplesheet

if __name__ == "__main__":
    samplesheet.csv2yaml(sys.argv[1])
