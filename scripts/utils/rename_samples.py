#!/usr/bin/env python
"""Rename sample name in a BAM file, eliminating spaces and colon characters.

Usage:
    rename_samples.py [<one or more> <input BAM files>]

Requires:
    pysam -- http://code.google.com/p/pysam/
      % sudo easy_install Cython
      % wget http://pysam.googlecode.com/files/pysam-0.3.tar.gz
      % tar -xzvpf pysam-0.3.tar.gz
      % cd pysam-0.3
      % python setup.py build && sudo python setup.py install
"""
import os
import sys

import pysam

def main(in_files):
    for in_file in in_files:
        out_file = "%s.rename" % in_file
        backup_file = "%s.orig" % in_file
        orig = pysam.Samfile(in_file, "rb")
        new_header = orig.header
        new_header["RG"][0]["SM"] = new_header["RG"][0]["SM"].split(": ")[-1]
        new = pysam.Samfile(out_file, "wb", header=new_header)
        for read in orig:
            new.write(read)
        orig.close()
        new.close()
        os.rename(in_file, backup_file)
        os.rename(out_file, in_file)

if __name__ == "__main__":
    main(sys.argv[1:])
