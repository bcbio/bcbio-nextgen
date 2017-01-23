"""Report version information for data used in an analysis.
"""
import csv
import os

from bcbio import install, utils

def write_versions(dirs, items):
    """Write data versioning for genomes present in the configuration.
    """
    genomes = sorted(list(set([x["genome_build"] for x in items])))
    out_file = _get_out_file(dirs)
    found_versions = False
    if genomes and out_file:
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["genome", "resource", "version"])
            for genome in genomes:
                genome_dir = install.get_genome_dir(genome, dirs.get("galaxy"), items[0])
                if genome_dir:
                    version_file = os.path.join(genome_dir, "versions.csv")
                    if os.path.exists(version_file):
                        found_versions = True
                        with open(version_file) as in_handle:
                            reader = csv.reader(in_handle)
                            for resource, version in reader:
                                writer.writerow([genome, resource, version])
    if found_versions:
        return out_file

def _get_out_file(dirs):
    if dirs and dirs.get("work"):
        base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
        return os.path.join(base_dir, "data_versions.csv")
