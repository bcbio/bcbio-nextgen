"""Report version information for data used in an analysis.
"""
import csv
import glob
import os

from bcbio import install, utils

def write_versions(dirs, items):
    """Write data versioning for genomes present in the configuration.
    """
    genomes = sorted(list(set([x["genome_build"] for x in items])))
    genome_dir = install.get_genomes_dir()
    out_file = _get_out_file(dirs)
    if out_file:
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["genome", "resource", "version"])
            for genome in genomes:
                in_version_files = glob.glob(os.path.join(genome_dir, "*", genome, "versions.csv"))
                if len(in_version_files) == 1:
                    with open(in_version_files[0]) as in_handle:
                        reader = csv.reader(in_handle)
                        for resource, version in reader:
                            writer.writerow([genome, resource, version])
    return out_file

def _get_out_file(dirs):
    if dirs and dirs.get("work"):
        base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
        return os.path.join(base_dir, "data_versions.csv")
