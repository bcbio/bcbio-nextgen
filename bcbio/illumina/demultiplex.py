"""Demultiplex and fastq conversion from Illumina output directories.

Uses Illumina's bcl2fastq: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.ilmn
"""
import os
import subprocess

from bcbio import utils

def run_bcl2fastq(run_folder, ss_csv, config):
    """Run bcl2fastq for de-multiplexing and fastq generation.
    run_folder -- directory of Illumina outputs
    ss_csv -- Samplesheet CSV file describing samples.
    """
    bc_dir = os.path.join(run_folder, "Data", "Intensities", "BaseCalls")
    output_dir = os.path.join(run_folder, "fastq")

    if not os.path.exists(output_dir):
        subprocess.check_call(["configureBclToFastq.pl", "--no-eamss",
                               "--input-dir", bc_dir, "--output-dir", output_dir,
                               "--sample-sheet", ss_csv])
    with utils.chdir(output_dir):
        subprocess.check_call(["make", "-j", str(utils.get_in(config, ("algorithm", "num_cores"), 1))])
    return output_dir
