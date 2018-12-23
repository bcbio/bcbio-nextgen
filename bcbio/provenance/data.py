"""Report version information for data used in an analysis.
"""
import csv
import os

from bcbio import install, utils

def write_versions(dirs, items):
    """Write data versioning for genomes present in the configuration.
    """
    genomes = {}
    for d in items:
        genomes[d["genome_build"]] = d.get("reference", {}).get("versions")
    out_file = _get_out_file(dirs)
    found_versions = False
    if genomes and out_file:
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle)
            writer.writerow(["genome", "resource", "version"])
            for genome, version_file in genomes.items():
                if not version_file:
                    genome_dir = install.get_genome_dir(genome, dirs.get("galaxy"), items[0])
                    if genome_dir:
                        version_file = os.path.join(genome_dir, "versions.csv")
                if version_file and os.path.exists(version_file):
                    found_versions = True
                    with open(version_file) as in_handle:
                        reader = csv.reader(in_handle)
                        for parts in reader:
                            if len(parts) >= 2:
                                resource, version = parts[:2]
                                writer.writerow([genome, resource, version])
    if found_versions:
        return out_file

def _get_out_file(dirs):
    if dirs and dirs.get("work"):
        base_dir = utils.safe_makedir(os.path.join(dirs["work"], "provenance"))
        return os.path.join(base_dir, "data_versions.csv")
