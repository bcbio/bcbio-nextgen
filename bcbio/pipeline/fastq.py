"""Pipeline utilities to retrieve
"""
import os
import glob
import subprocess

def get_fastq_files(directory, lane, fc_name, bc_name=None):
    """Retrieve fastq files for the given lane, ready to process.
    """
    if bc_name:
        glob_str = "%s_*%s_%s_*_fastq.txt" % (lane, fc_name, bc_name)
    else:
        glob_str = "%s_*%s*_fastq.txt" % (lane, fc_name)
    files = glob.glob(os.path.join(directory, glob_str))
    files.sort()
    if len(files) > 2 or len(files) == 0:
        raise ValueError("Did not find correct files for %s %s %s %s" %
                (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz"):
            cl = ["gunzip", fname]
            subprocess.check_call(cl)
            ready_files.append(os.path.splitext(fname)[0])
        else:
            ready_files.append(fname)
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)

