"""Tests associated with detecting sequencing results dumped from a machine.
"""
import os

import yaml

import pytest
from bcbio.illumina import samplesheet


class TestSampleSheet(object):
    """Deal with Illumina SampleSheets and convert to YAML input.
    """
    @property
    def ss_file(self):
        return os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "data", "illumina_samplesheet.csv")

    @pytest.mark.speed1
    def test_toyaml(self):
        """Convert CSV Illumina SampleSheet to YAML.
        """
        out_file = samplesheet.csv2yaml(self.ss_file)
        assert os.path.exists(out_file)
        with open(out_file) as in_handle:
            info = yaml.safe_load(in_handle)
        assert info[0]['lane'] == '1'
        assert info[0]['multiplex'][0]['barcode_id'] == 5
        os.remove(out_file)

    @pytest.mark.speed1
    def test_checkforrun(self):
        """Check for the presence of runs in an Illumina SampleSheet.
        """
        fcdir = "fake/101007_80HM7ABXX"
        config = {"samplesheet_directories": [os.path.dirname(self.ss_file)]}
        ss = samplesheet.run_has_samplesheet(fcdir, config, False)
        assert ss is not None
        fcdir = "fake/101007_NOPEXX"
        ss = samplesheet.run_has_samplesheet(fcdir, config, False)
        assert ss is None
