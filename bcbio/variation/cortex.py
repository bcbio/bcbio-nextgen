"""Perform regional de-novo assembly calling with cortex_var.

Using a pre-mapped set of reads and BED file of regions, performs de-novo
assembly and variant calling against the reference sequence in each region.
This avoids whole genome costs while gaining the advantage of de-novo
prediction.

http://cortexassembler.sourceforge.net/index_cortex_var.html
"""
import os
import glob
import subprocess
from contextlib import closing

import pysam
from Bio import Seq

from bcbio import broad
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.utils import file_exists

def run_cortex(align_bam, ref_file, config, dbsnp=None, region=None,
               out_file=None):
    """Top level entry to regional de-novo based variant calling with cortex_var.
    """
    if out_file is None:
        out_file = "%s-cortex.vcf" % os.path.splitext(align_bam)[0]
    if not file_exists(out_file):
        variant_regions = config["algorithm"].get("variant_regions", None)
        if not variant_regions:
            raise ValueError("Only regional variant calling with cortex_var is supported. Set variant_regions")
        target_regions = subset_variant_regions(variant_regions, region, out_file)
        with open(target_regions) as in_handle:
            regional_vcfs = [_run_cortex_on_region(x.strip().split("\t")[:3], align_bam,
                                                   ref_file, out_file, config)
                             for x in in_handle]
    return out_file

def _run_cortex_on_region(region, align_bam, ref_file, out_file_base, config):
    """Run cortex on a specified chromosome start/end region.
    """
    kmers = [31]
    cortex_dir = config["program"].get("cortex")
    stampy_dir = config["program"].get("stampy")
    vcftools_dir = config["program"].get("vcftools")
    if cortex_dir is None or stampy_dir is None:
        raise ValueError("cortex_var requires path to pre-built cortex and stampy")
    out_vcf_base = apply("{0}-{1}-{2}-{3}".format, [os.path.splitext(out_file_base)[0]] + region)
    fastq = _get_fastq_in_region(region, align_bam, out_vcf_base)
    local_ref, genome_size = _get_local_ref(region, ref_file, out_vcf_base)
    indexes = _index_local_ref(local_ref, cortex_dir, stampy_dir, kmers)
    _run_cortex(fastq, indexes, {"kmers": kmers, "genome_size": genome_size},
                out_vcf_base, {"cortex": cortex_dir, "stampy": stampy_dir,
                               "vcftools": vcftools_dir},
                config)
    
    print region, align_bam, fastq, indexes
    raise NotImplementedError

def _run_cortex(fastq, indexes, params, out_base, dirs, config):
    """Run cortex_var run_calls.pl, producing a VCF variant file.
    """
    assert len(params["kmers"]) == 1, "Currently only support single kmer workflow"
    fastaq_index = "{0}.fastaq_index".format(out_base)
    se_fastq_index = "{0}.se_fastq".format(out_base)
    pe_fastq_index = "{0}.pe_fastq".format(out_base)
    reffasta_index = "{0}.list_ref_fasta".format(out_base)
    with open(se_fastq_index, "w") as out_handle:
        out_handle.write(fastq + "\n")
    with open(pe_fastq_index, "w") as out_handle:
        out_handle.write("")
    with open(fastaq_index, "w") as out_handle:
        out_handle.write("{0}\t{1}\t{2}\t{2}\n".format(os.path.basename(out_base), se_fastq_index, pe_fastq_index))
    with open(reffasta_index, "w") as out_handle:
        for x in indexes["cortex"]:
            out_handle.write(x + "\n")
    subprocess.check_call(["perl", os.path.join(dirs["cortex"], "scripts", "calling", "run_calls.pl"),
                           "--first_kmer", str(params["kmers"][0]), "--fastaq_index", fastaq_index,
                           "--auto_cleaning", "yes", "--bc", "yes", "--pd", "yes",
                           "--outdir", os.path.dirname(out_base), "--outvcf", os.path.basename(out_base),
                           "--ploidy", str(config["algorithm"].get("ploidy", 2)),
                           "--stampy_hash", indexes["stampy"],
                           "--stampy_bin", os.path.join(dirs["stampy"], "stampy.py"),
                           "--refbindir", os.path.dirname(indexes["cortex"][0]),
                           "--list_ref_fasta",  reffasta_index,
                           "--genome_size", str(params["genome_size"]),
                           "--max_read_len", "10000",
                           "--format", "FASTQ", "--qthresh", "5", "--do_union", "yes",
                           "--mem_height", "17", "--mem_width", "100",
                           "--ref", "CoordinatesAndInCalling", "--workflow", "independent",
                           "--vcftools_dir", dirs["vcftools"],
                           "--logfile", "{0}.logfile,f".format(out_base)])

def _get_cortex_binary(kmer, cortex_dir):
    cortex_bin = None
    for check_bin in sorted(glob.glob(os.path.join(cortex_dir, "bin", "cortex_var_*"))):
        kmer_check = int(os.path.basename(check_bin).split("_")[2])
        if kmer_check >= kmer:
            cortex_bin = check_bin
            break
    assert cortex_bin is not None, \
        "Could not find cortex_var executable in %s for kmer %s" % (cortex_dir, kmer)
    return cortex_bin

def _index_local_ref(fasta_file, cortex_dir, stampy_dir, kmers):
    """Pre-index a generated local reference sequence with cortex_var and stampy.
    """
    base_out = os.path.splitext(fasta_file)[0]
    cindexes = []
    for kmer in kmers:
        out_file = "{0}.k{1}.ctx".format(base_out, kmer)
        if not file_exists(out_file):
            file_list = "{0}.se_list".format(base_out)
            with open(file_list, "w") as out_handle:
                out_handle.write(fasta_file + "\n")
            subprocess.check_call([_get_cortex_binary(kmer, cortex_dir),
                                   "--kmer_size", str(kmer), "--mem_height", "17",
                                   "--se_list", file_list, "--format", "FASTA",
                                   "--max_read_len", "10000", "--sample_id", base_out,
                                   "--dump_binary", out_file])
        cindexes.append(out_file)
    if not file_exists("{0}.stidx".format(base_out)):
        subprocess.check_call([os.path.join(stampy_dir, "stampy.py"), "-G",
                               base_out, fasta_file])
        subprocess.check_call([os.path.join(stampy_dir, "stampy.py"), "-g",
                               base_out, "-H", base_out])
    return {"stampy": base_out,
            "cortex": cindexes}

def _get_local_ref(region, ref_file, out_vcf_base):
    """Retrieve a local FASTA file corresponding to the specified region.
    """
    out_file = "{0}.fa".format(out_vcf_base)
    if not file_exists(out_file):
        with closing(pysam.Fastafile(ref_file)) as in_pysam:
            contig, start, end = region
            seq = in_pysam.fetch(contig, int(start) - 1, int(end))
            with open(out_file, "w") as out_handle:
                out_handle.write(">{0}-{1}-{2}\n{3}".format(contig, start, end,
                                                              str(seq)))
    with open(out_file) as in_handle:
        in_handle.readline()
        size = len(in_handle.readline().strip())
    return out_file, size

def _get_fastq_in_region(region, align_bam, out_base):
    """Retrieve fastq files in region as single end.
    Paired end is more complicated since pairs can map off the region, so focus
    on local only assembly since we've previously used paired information for mapping.
    """
    out_file = "{0}.fastq".format(out_base)
    if not file_exists(out_file):
        with closing(pysam.Samfile(align_bam, "rb")) as in_pysam:
            with file_transaction(out_file) as tx_out_file:
                with open(out_file, "w") as out_handle:
                    contig, start, end = region
                    for read in in_pysam.fetch(contig, int(start) - 1, int(end)):
                        seq = Seq.Seq(read.seq)
                        qual = list(read.qual)
                        if read.is_reverse:
                            seq = seq.reverse_complement()
                            qual.reverse()
                        out_handle.write("@{name}\n{seq}\n+\n{qual}\n".format(
                                name=read.qname, seq=str(seq), qual="".join(qual)))
    return out_file
