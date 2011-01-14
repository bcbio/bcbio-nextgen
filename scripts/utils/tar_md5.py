#!/usr/bin/env python
"""Script to tar a dataset for archival.

It generates a TAR archive while computing MD5 checksum for each file.

From VeriTAR:
http://www.g-loaded.eu/2007/12/01/veritar-verify-checksums-of-files-within-a-tar-archive/

ToDo: Stream the backup to archival system using ssh/rsync while tarring (more pipework)
"""
import os
import fnmatch
import sys
import tarfile
import hashlib
import pprint
from subprocess import *

# Subclassing attempt on TarFile, "with" not supported yet on 2.6
#class myTarFile(tarfile):
#    def __exit__(self, *args):
#        self.close()

# From O'reilly Python Cookbook 2nd edition
def all_files(root, patterns='*', single_level=False, yield_folders=False):
# Expand patterns from semicolon-separated string to list
    patterns = patterns.split(';')
    for path, subdirs, files in os.walk(root):
        if yield_folders:
            files.extend(subdirs)
        files.sort()
        for name in files:
            for pattern in patterns:
                if fnmatch.fnmatch(name, pattern):
                    yield os.path.join(path, name)
                    break
        if single_level:
            break

def make_hashed_tar(tar_filename, tree_to_backup):
    tar = tarfile.open("%s.tar" % tar_filename, "w")
    hashfile = open("%s.md5" % tar_filename, "wa")
    for item in list(tree_to_backup):
        tar.add(item)
        if os.path.isfile(item):
            hashfile.write(hashlib.md5(item).hexdigest()+"\t"+item)
            hashfile.write("\n")
    tar.close()
    hashfile.close()

    #hashfile = open("%s.md5" % filename, "w")
    #with myTarFile.open("%s.tar" % filename, "w"):
    #    for f in dirtotar:
    #        tar.add(f)
    #        if f.isfile():
    #            hashfile.write(hashlib.md5(f).hexdigest())

if __name__ == "__main__":
    dataset = sys.argv[1]
    make_hashed_tar(dataset, all_files(dataset))
