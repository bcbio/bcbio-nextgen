"""
functions for finding, retrieving data from and updating the Galaxy .loc files
"""
import os
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists
import tempfile
import shutil

REF_FILES = {"bwa": "bwa_index.loc",
             "bowtie": "bowtie_indices.loc",
             "samtools": "sam_fa_indices.loc",
             "maq": "bowtie_indices.loc",
             "seq": "alignseq.loc",
             "bowtie2": "bowtie2_indices.loc"}

def samtools_formatter(build, loc, name=None):
    return "\t".join(["index", build, loc]) + "\n"

def ucsc_formatter(build, loc, name=None):
    return "\t".join(["seq", build, loc]) + "\n"

def generic_formatter(build, loc, name=None):
    name = build if not name else name
    return "\t".join([build, build, name, loc]) + "\n"

def get_locformatter(loc_type):
    if loc_type in ["samtools", "ucsc"]:
        return eval("%s_formatter" % loc_type)
    else:
        return generic_formatter

def get_loc_file(galaxy_base, loc_type):
    loc_file = REF_FILES.get(loc_type, None)
    if not loc_file:
        return None
    return os.path.abspath(os.path.join(galaxy_base, "tool-data", loc_file))

def get_loc_files(galaxy_base):
    """
    get dictionary of loc_type: loc_file, .loc files in the galaxy base
    for example: {"bwa": "/galaxy_base_path/tool-dir/bwa_index.loc"}
    """
    return {k: os.path.join(galaxy_base, "tool-data", v) for k, v in REF_FILES.items()}

def get_genome_refs(loc_file, loc_type):
    """
    get dictionary of genome: location for all genomes of type in a .loc file
    for example: {'hg19': '/genomedir/Hsapiens/hg19/seq/hg19.fa'}
    """
    if not file_exists(out_file):
        return None
    refs = {}
    with open(loc_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                parts = line.strip().split()
                if loc_type in ["bowtie2", "samtools", "alignseq"]:
                    refs[parts[1]] = parts[-1]
                else:
                    refs[parts[0]] = parts[-1]
    return refs

def update_loc_file(galaxy_base, loc_type, genome_build, ref_loc):
    ref_loc = os.path.abspath(ref_loc)
    loc_file = get_loc_file(galaxy_base, loc_type)
    if not loc_file:
        return None
    formatter = get_locformatter(loc_type)
    builds = []
    tmp_out = tempfile.NamedTemporaryFile(delete=False).name
    if file_exists(loc_file):
        with open(loc_file) as in_handle, open(tmp_out, "w") as out_handle:
            for line in in_handle:
                if line.startswith("#"):
                    out_handle.write(line)
                else:
                    parts = line.strip().split()
                    build = parts[1]
                    builds.append(build)
                    if build != genome_build:
                        out_handle.write(line)
                    else:
                        out_handle.write(formatter(genome_build, ref_loc))
        shutil.copyfile(tmp_out, loc_file)
    if genome_build not in builds or not file_exists(loc_file):
        with open(loc_file, "a") as out_handle:
            out_handle.write(formatter(genome_build, ref_loc))
    return loc_file
