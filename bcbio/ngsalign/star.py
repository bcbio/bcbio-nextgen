from os import path
import tempfile

from bcbio.pipeline import config_utils
from bcbio.utils import safe_makedir, file_exists, get_in
from bcbio.provenance import do
from bcbio.bam import fastq

CLEANUP_FILES = ["Aligned.out.sam", "Log.out", "Log.progress.out"]



def align(fastq_file, pair_file, ref_file, out_base, align_dir, data,
          names=None):
    genome_dir = create_genome(fastq_file, data)

    config = data["config"]
    out_prefix = path.join(align_dir, out_base)
    out_file = out_prefix + "Aligned.out.sam"
    if file_exists(out_file):
        return out_file
    star_path = config_utils.get_program("STAR", config)
    fastq = " ".join([fastq_file, pair_file])
    num_cores = config["algorithm"].get("num_cores", 1)

    safe_makedir(align_dir)
    cmd = ("{star_path} --genomeDir {genome_dir} --readFilesIn {fastq} "
           "--runThreadN {num_cores} --outFileNamePrefix {out_prefix} "
           "--outReadsUnmapped Fastx --outFilterMultimapNmax 10")
    fusion_mode = get_in(data, ("config", "algorithm", "fusion_mode"), False)
    if fusion_mode:
        cmd += " --chimSegmentMin 15 --chimJunctionOverhangMin 15"
    strandedness = get_in(data, ("config", "algorithm", "strandedness"),
                          "unstranded").lower()
    if strandedness == "unstranded":
        cmd += " --outSAMstrandField intronMotif"
    run_message = "Running STAR aligner on %s and %s." % (pair_file, ref_file)
    do.run(cmd.format(**locals()), run_message, None)
    return out_file

def _get_quality_format(config):
    qual_format = config["algorithm"].get("quality_format", None)
    if qual_format.lower() == "illumina":
        return "fastq-illumina"
    elif qual_format.lower() == "solexa":
        return "fastq-solexa"
    else:
        return "fastq-sanger"

def remap_index_fn(ref_file):
    """Map sequence references to equivalent star indexes
    """
    return path.join(path.dirname(path.dirname(ref_file)), "star")

def create_genome(fastq_file, data):
    """
    create a genome file with splicing information. this depends on read length
    so must be made custom for each mapping
    """
    config = data["config"]
    star_path = config_utils.get_program("star", config, default="STAR")
    tempdir = tempfile.mkdtemp("_star")
    quality_format = _get_quality_format(data["config"])
    overhang = fastq.estimate_read_length(fastq_file, quality_format, 100000) - 1
    fasta_file = data["sam_ref"]
    gtf_file = get_in(data, ("genome_resources", "rnaseq", "transcripts"), None)
    num_cores = config["algorithm"].get("num_cores", 1)

    mem_limit = 10000000000

    cmd = ("{star_path} --genomeDir {tempdir} --genomeFastaFiles {fasta_file} "
           " --runThreadN {num_cores} --runMode genomeGenerate "
           "--sjdbOverhang {overhang} --limitGenomeGenerateRAM {mem_limit}")
    if gtf_file:
        cmd + " --sjdbGTFfile {gtf_file}"
    run_message = ("Generating a STAR index for %s using %s as transcripts with an "
                   "overhang of %d") % (fasta_file, gtf_file, overhang)

    do.run(cmd.format(**locals()), run_message, None)
    return tempdir

def job_requirements(cores, memory):
    MIN_STAR_MEMORY = 30.0
    if not memory or cores * memory < MIN_STAR_MEMORY:
        memory = MIN_STAR_MEMORY / cores
    return cores, memory

align.job_requirements = job_requirements
