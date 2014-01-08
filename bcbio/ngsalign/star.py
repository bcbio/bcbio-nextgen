from os import path, listdir

from bcbio.pipeline import config_utils
from bcbio.utils import safe_makedir, file_exists, get_in
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction

CLEANUP_FILES = ["Aligned.out.sam", "Log.out", "Log.progress.out"]



def align(fastq_file, pair_file, ref_file, out_base, align_dir, data,
          names=None):
    config = data["config"]
    out_prefix = path.join(align_dir, out_base)
    out_file = out_prefix + "Aligned.out.sam"
    if file_exists(out_file):
        return out_file
    star_path = config_utils.get_program("STAR", config)
    fastq = " ".join([fastq_file, pair_file])
    num_cores = config["algorithm"].get("num_cores", 1)

    safe_makedir(align_dir)
    cmd = ("{star_path} --genomeDir {ref_file} --readFilesIn {fastq} "
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

def index(item):
    """
    create a genome file with splicing information
    """
    genome_dir = item["align_ref"]
    if path.exists(genome_dir) and listdir(genome_dir):
        return genome_dir

    config = item["config"]
    star_path = config_utils.get_program("star", config, default="STAR")
    overhang = 99 # default to 100 bp reads
    fasta_file = item["sam_ref"]
    gtf_file = get_in(item, ("genome_resources", "rnaseq", "transcripts"), None)
    num_cores = config["algorithm"].get("num_cores", 1)

    with file_transaction(genome_dir) as tx_genome_dir:
        safe_makedir(tx_genome_dir)
        cmd = ("{star_path} --genomeDir {tx_genome_dir} --genomeFastaFiles {fasta_file} "
               " --runThreadN {num_cores} --runMode genomeGenerate "
               "--sjdbOverhang {overhang} ")
        if gtf_file:
            cmd + " --sjdbGTFfile {gtf_file}"
        run_message = ("Generating a STAR index for %s using %s as transcripts with an "
                       "overhang of %d") % (fasta_file, gtf_file, overhang)
        do.run(cmd.format(**locals()), run_message, None)
    return genome_dir

def job_requirements(cores, memory):
    MIN_STAR_MEMORY = 30.0
    if not memory or cores * memory < MIN_STAR_MEMORY:
        memory = MIN_STAR_MEMORY / cores
    return cores, memory

align.job_requirements = job_requirements
