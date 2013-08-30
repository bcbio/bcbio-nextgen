"""Pipeline code to run alignments and prepare BAM files.

This works as part of the lane/flowcell process step of the pipeline.
"""
from collections import namedtuple
import ConfigParser
import os
from xml.etree import ElementTree

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio import utils, broad
from bcbio.bam import cram
from bcbio.ngsalign import (bowtie, bwa, tophat, bowtie2, mosaik,
                            novoalign, star)
from bcbio.distributed.transaction import file_transaction

# Define a next-generation sequencing tool to plugin:
# align_fn -- runs an aligner and generates SAM output
# galaxy_loc_file -- name of a Galaxy location file to retrieve
#  the genome index location
# bam_align_fn -- runs an aligner on a BAM file
# remap_index_fn -- Function that will take the location provided
#  from galaxy_loc_file and find the actual location of the index file.
#  This is useful for indexes that don't have an associated location file
#  but are stored in the same directory structure.
NgsTool = namedtuple("NgsTool", ["align_fn", "pipe_align_fn", "bam_align_fn",
                                 "galaxy_loc_file", "remap_index_fn", "can_pipe"])

BASE_LOCATION_FILE = "sam_fa_indices.loc"
_tools = {
    "bowtie": NgsTool(bowtie.align, None, None, bowtie.galaxy_location_file, None, None),
    "bowtie2": NgsTool(bowtie2.align, None, None, bowtie2.galaxy_location_file, bowtie2.remap_index_fn,
                       None),
    "bwa": NgsTool(bwa.align, bwa.align_pipe, bwa.align_bam, bwa.galaxy_location_file, None,
                   bwa.can_pipe),
    "mosaik": NgsTool(mosaik.align, None, None, mosaik.galaxy_location_file, None,
                      None),
    "novoalign": NgsTool(novoalign.align, novoalign.align_pipe, novoalign.align_bam,
                         novoalign.galaxy_location_file, novoalign.remap_index_fn, novoalign.can_pipe),
    "tophat": NgsTool(tophat.align, None, None, bowtie2.galaxy_location_file, bowtie2.remap_index_fn,
                      None),
    "samtools": NgsTool(None, None, None, BASE_LOCATION_FILE, None, None),
    "star": NgsTool(star.align, None, None, None, star.remap_index_fn, None),
    "tophat2": NgsTool(tophat.align, None, None, bowtie2.galaxy_location_file, bowtie2.remap_index_fn,
                      None)}

metadata = {"support_bam": [k for k, v in _tools.iteritems() if v.bam_align_fn is not None]}

def align_to_sort_bam(fastq1, fastq2, names, genome_build, aligner,
                      dirs, config, dir_ext=""):
    """Align to the named genome build, returning a sorted BAM file.
    """
    align_dir = utils.safe_makedir(os.path.join(dirs["work"], "align", names["sample"], dir_ext))
    align_ref, sam_ref = get_genome_ref(genome_build, aligner, dirs["galaxy"])
    if fastq1.endswith(".bam"):
        out_bam = _align_from_bam(fastq1, aligner, align_ref, sam_ref, names, align_dir, config)
    elif _can_pipe(aligner, fastq1):
        out_bam = _align_from_fastq_pipe(fastq1, fastq2, aligner, align_ref, sam_ref, names,
                                         align_dir, config)
    else:
        out_bam = _align_from_fastq(fastq1, fastq2, aligner, align_ref, sam_ref, names,
                                    align_dir, config)
    return out_bam, sam_ref

def _can_pipe(aligner, fastq_file):
    """Check if current aligner support piping for a particular input fastq file.
    """
    if _tools[aligner].can_pipe and _tools[aligner].pipe_align_fn:
        return _tools[aligner].can_pipe(fastq_file)
    return False

def _align_from_fastq_pipe(fastq1, fastq2, aligner, align_ref, sam_ref, names, align_dir, config):
    """Align longer reads using new piped strategies that avoid disk IO.
    """
    align_fn = _tools[aligner].pipe_align_fn
    if align_fn is None:
        raise NotImplementedError("Do not yet support piped alignment with %s" % aligner)
    return align_fn(fastq1, fastq2, align_ref, names, align_dir, config)

def _align_from_bam(fastq1, aligner, align_ref, sam_ref, names, align_dir, config):
    qual_bin_method = config["algorithm"].get("quality_bin")
    if (qual_bin_method == "prealignment" or
         (isinstance(qual_bin_method, list) and "prealignment" in qual_bin_method)):
        out_dir = utils.safe_makedir(os.path.join(align_dir, "qualbin"))
        fastq1 = cram.illumina_qual_bin(fastq1, sam_ref, out_dir, config)
    align_fn = _tools[aligner].bam_align_fn
    if align_fn is None:
        raise NotImplementedError("Do not yet support BAM alignment with %s" % aligner)
    return align_fn(fastq1, align_ref, names, align_dir, config)

def _align_from_fastq(fastq1, fastq2, aligner, align_ref, sam_ref, names,
                      align_dir, config):
    """Align from fastq inputs, producing sorted BAM output.
    """
    align_fn = _tools[aligner].align_fn
    sam_file = align_fn(fastq1, fastq2, align_ref, names["lane"], align_dir, config,
                        names)
    if fastq2 is None and aligner in ["bwa", "bowtie2"]:
        fastq1 = _remove_read_number(fastq1, sam_file)
    sort_method = config["algorithm"].get("bam_sort", "coordinate")

    if sort_method == "queryname":
        return sam_to_querysort_bam(sam_file, config)
    else:
        return sam_to_sort_bam(sam_file, sam_ref, fastq1, fastq2, names, config)

def _remove_read_number(in_file, sam_file):
    """Work around problem with MergeBamAlignment with BWA and single end reads.

    Need to remove read number ends from Fastq to match BWA stripping of numbers.

    http://sourceforge.net/mailarchive/forum.php?thread_name=87bosvbbqz.fsf%
    40fastmail.fm&forum_name=samtools-help
    http://sourceforge.net/mailarchive/forum.php?thread_name=4EB03C42.2060405%
    40broadinstitute.org&forum_name=samtools-help
    """
    out_file = os.path.join(os.path.dirname(sam_file),
                            "%s-safe%s" % os.path.splitext(os.path.basename(in_file)))
    if not os.path.exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            with open(in_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for i, (name, seq, qual) in enumerate(FastqGeneralIterator(in_handle)):
                        if i == 0 and not name.endswith("/1"):
                            out_file = in_file
                            break
                        else:
                            name = name.rsplit("/", 1)[0]
                            out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
    return out_file

def sam_to_querysort_bam(sam_file, config):
    """Convert SAM file directly to a query sorted BAM without merging of FASTQ reads.

    This allows merging of multiple mappers which do not work with MergeBamAlignment.
    """
    runner = broad.runner_from_config(config)
    out_file = "{}-querysorted.bam".format(os.path.splitext(sam_file)[0])
    return runner.run_fn("picard_sort", sam_file, "queryname", out_file)

def sam_to_querysort_sam(sam_file, config):
    """Convert SAM file directly to a query sorted SAM without merging of FASTQ reads.

    This allows merging of multiple mappers which do not work with MergeBamAlignment.
    """
    runner = broad.runner_from_config(config)
    out_file = "{}-querysorted.sam".format(os.path.splitext(sam_file)[0])
    return runner.run_fn("picard_sort", sam_file, "queryname", out_file)

def sam_to_sort_bam(sam_file, ref_file, fastq1, fastq2, names, config):
    """Convert SAM file to merged and sorted BAM file.
    """
    picard = broad.runner_from_config(config)
    base_dir = os.path.dirname(sam_file)

    picard.run_fn("picard_index_ref", ref_file)
    out_fastq_bam = picard.run_fn("picard_fastq_to_bam", fastq1, fastq2, base_dir, names)
    out_bam = picard.run_fn("picard_sam_to_bam", sam_file, out_fastq_bam, ref_file,
                            fastq2 is not None)
    sort_bam = picard.run_fn("picard_sort", out_bam)

    utils.save_diskspace(sam_file, "SAM converted to BAM", config)
    utils.save_diskspace(out_fastq_bam, "Combined into output BAM %s" % out_bam, config)
    utils.save_diskspace(out_bam, "Sorted to %s" % sort_bam, config)
    # merge FASTQ files, only if barcoded samples in the work directory
    if (os.path.commonprefix([fastq1, sort_bam]) ==
             os.path.split(os.path.dirname(sort_bam))[0]
          and not config["algorithm"].get("upload_fastq", True)):
        utils.save_diskspace(fastq1, "Merged into output BAM %s" % out_bam, config)
        if fastq2:
            utils.save_diskspace(fastq2, "Merged into output BAM %s" % out_bam, config)
    return sort_bam

# ## Galaxy integration -- *.loc files

def _get_galaxy_loc_file(name, galaxy_dt, ref_dir, galaxy_base):
    """Retrieve Galaxy *.loc file for the given reference/aligner name.

    First tries to find an aligner specific *.loc file. If not defined
    or does not exist, then we need to try and remap it from the
    default reference file
    """
    if "file" in galaxy_dt and os.path.exists(os.path.join(galaxy_base, galaxy_dt["file"])):
        loc_file = os.path.join(galaxy_base, galaxy_dt["file"])
        need_remap = False
    elif _tools[name].galaxy_loc_file is None:
        loc_file = os.path.join(ref_dir, BASE_LOCATION_FILE)
        need_remap = True
    else:
        loc_file = os.path.join(ref_dir, _tools[name].galaxy_loc_file)
        need_remap = False
    if not os.path.exists(loc_file):
        loc_file = os.path.join(ref_dir, BASE_LOCATION_FILE)
        need_remap = True
    return loc_file, need_remap

def _get_ref_from_galaxy_loc(name, genome_build, loc_file, galaxy_dt, need_remap):
    """Retrieve reference genome file from Galaxy *.loc file.

    Reads from tool_data_table_conf.xml information for the index if it
    exists, otherwise uses heuristics to find line based on most common setups.
    """
    if "column" in galaxy_dt:
        dbkey_i = galaxy_dt["column"].index("dbkey")
        path_i = galaxy_dt["column"].index("path")
    else:
        dbkey_i = None
    cur_ref = None
    with open(loc_file) as in_handle:
        for line in in_handle:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                if len(parts) == 1: # spaces instead of tabs
                    parts = [x.strip() for x in line.strip().split("  ") if x.strip()]
                if dbkey_i is not None and not need_remap:
                    if parts[dbkey_i] == genome_build:
                        cur_ref = parts[path_i]
                        break
                else:
                    if parts[0] == "index":
                        parts = parts[1:]
                    if parts[0] == genome_build:
                        cur_ref = parts[-1]
                        break
    if cur_ref is None:
        raise IndexError("Genome %s not found in %s" % (genome_build, loc_file))
    if need_remap:
        remap_fn = _tools[name].remap_index_fn
        assert remap_fn is not None, "%s requires remapping function from base location file" % name
        cur_ref = remap_fn(cur_ref)
    return cur_ref

def _get_galaxy_tool_info(galaxy_base):
    """Retrieve Galaxy tool-data information from defaults or galaxy config file.
    """
    ini_file = os.path.join(galaxy_base, "universe_wsgi.ini")
    info = {"tool_data_table_config_path": os.path.join(galaxy_base, "tool_data_table_conf.xml"),
            "tool_data_path": os.path.join(galaxy_base, "tool-data")}
    config = ConfigParser.ConfigParser()
    config.read(ini_file)
    if "app:main" in config.sections():
        for option in config.options("app:main"):
            if option in info:
                info[option] = os.path.join(galaxy_base, config.get("app:main", option))
    return info

def _get_galaxy_data_table(name, dt_config_file):
    """Parse data table config file for details on tool *.loc location and columns.
    """
    out = {}
    if os.path.exists(dt_config_file):
        tdtc = ElementTree.parse(dt_config_file)
        for t in tdtc.getiterator("table"):
            if t.attrib.get("name", "") in [name, "%s_indexes" % name]:
                out["column"] = [x.strip() for x in t.find("columns").text.split(",")]
                out["file"] = t.find("file").attrib.get("path", "")
    return out

def get_genome_ref(genome_build, aligner, galaxy_base):
    """Retrieve the reference genome file location from galaxy configuration.
    """
    if not genome_build:
        return (None, None)
    galaxy_config = _get_galaxy_tool_info(galaxy_base)
    out_info = []
    for name in [aligner, "samtools"]:
        if not name:
            out_info.append(None)
            continue
        galaxy_dt = _get_galaxy_data_table(name, galaxy_config["tool_data_table_config_path"])
        loc_file, need_remap = _get_galaxy_loc_file(name, galaxy_dt, galaxy_config["tool_data_path"],
                                                    galaxy_base)
        cur_ref = _get_ref_from_galaxy_loc(name, genome_build, loc_file, galaxy_dt, need_remap)
        out_info.append(utils.add_full_path(cur_ref, galaxy_config["tool_data_path"]))

    if len(out_info) != 2:
        raise ValueError("Did not find genome reference for %s %s" %
                (genome_build, aligner))
    else:
        return tuple(out_info)
