"""High level code for driving a next-gen analysis pipeline.

This structures processing steps into the following modules:

  - lane.py: Analyze a single fastq file.
    - fastq.py: Utilities to retrieve fastq files.
    - alignment.py: Align to a reference genome.

  - sample.py: Analyze a sample, which may consist of multiple lanes or
               barcoded samples on a lane.
    - merge.py: Merge multiple sample files in one processing run.
    - variation.py: Calculate SNP/indel variations for a sample.
    - qcsummary.py: Quality control, alignment metrics and summary information.
"""
