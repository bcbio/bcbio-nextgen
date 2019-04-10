"""Perform regional de-novo assembly calling with cortex_var.

Using a pre-mapped set of reads and BED file of regions, performs de-novo
assembly and variant calling against the reference sequence in each region.
This avoids whole genome costs while gaining the advantage of de-novo
prediction.

http://cortexassembler.sourceforge.net/index_cortex_var.html
"""
from __future__ import print_function
import os
import glob
import subprocess
import itertools
import shutil

import pysam
from Bio import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from bcbio import bam
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import config_utils
from bcbio.pipeline.shared import subset_variant_regions
from bcbio.utils import file_exists, safe_makedir
from bcbio.variation import vcfutils

def run_cortex(align_bams, items, ref_file, assoc_files, region=None,
               out_file=None):
    """Top level entry to regional de-novo based variant calling with cortex_var.
    """
    raise NotImplementedError("Cortex currently out of date and needs reworking.")
    if len(align_bams) == 1:
        align_bam = align_bams[0]
        config = items[0]["config"]
    else:
        raise NotImplementedError("Need to add multisample calling for cortex_var")
    if out_file is None:
        out_file = "%s-cortex.vcf" % os.path.splitext(align_bam)[0]
    if region is not None:
        work_dir = safe_makedir(os.path.join(os.path.dirname(out_file),
                                             region.replace(".", "_")))
    else:
        work_dir = os.path.dirname(out_file)
    if not file_exists(out_file):
        bam.index(align_bam, config)
        variant_regions = config["algorithm"].get("variant_regions", None)
        if not variant_regions:
            raise ValueError("Only support regional variant calling with cortex_var: set variant_regions")
        target_regions = subset_variant_regions(variant_regions, region, out_file)
        if os.path.isfile(target_regions):
            with open(target_regions) as in_handle:
                regional_vcfs = [_run_cortex_on_region(x.strip().split("\t")[:3], align_bam,
                                                       ref_file, work_dir, out_file, config)
                                 for x in in_handle]

            combine_file = "{0}-raw{1}".format(*os.path.splitext(out_file))
            _combine_variants(regional_vcfs, combine_file, ref_file, config)
            _select_final_variants(combine_file, out_file, config)
        else:
            vcfutils.write_empty_vcf(out_file)
    return out_file

def _passes_cortex_depth(line, min_depth):
    """Do any genotypes in the cortex_var VCF line passes the minimum depth requirement?
    """
    parts = line.split("\t")
    cov_index = parts[8].split(":").index("COV")
    passes_depth = False
    for gt in parts[9:]:
        cur_cov = gt.split(":")[cov_index]
        cur_depth = sum(int(x) for x in cur_cov.split(","))
        if cur_depth >= min_depth:
            passes_depth = True
    return passes_depth

def _select_final_variants(base_vcf, out_vcf, config):
    """Filter input file, removing items with low depth of support.

    cortex_var calls are tricky to filter by depth. Count information is in
    the COV FORMAT field grouped by alleles, so we need to sum up values and
    compare.
    """
    min_depth = int(config["algorithm"].get("min_depth", 4))
    with file_transaction(out_vcf) as tx_out_file:
        with open(base_vcf) as in_handle:
            with open(tx_out_file, "w") as out_handle:
                for line in in_handle:
                    if line.startswith("#"):
                        passes = True
                    else:
                        passes = _passes_cortex_depth(line, min_depth)
                    if passes:
                        out_handle.write(line)
    return out_vcf

def _combine_variants(in_vcfs, out_file, ref_file, config):
    """Combine variant files, writing the header from the first non-empty input.

    in_vcfs is a list with each item starting with the chromosome regions,
    and ending with the input file.
    We sort by these regions to ensure the output file is in the expected order.
    """
    in_vcfs.sort()
    wrote_header = False
    with open(out_file, "w") as out_handle:
        for in_vcf in (x[-1] for x in in_vcfs):
            with open(in_vcf) as in_handle:
                header = list(itertools.takewhile(lambda x: x.startswith("#"),
                                                  in_handle))
                if not header[0].startswith("##fileformat=VCFv4"):
                    raise ValueError("Unexpected VCF file: %s" % in_vcf)
                for line in in_handle:
                    if not wrote_header:
                        wrote_header = True
                        out_handle.write("".join(header))
                    out_handle.write(line)
        if not wrote_header:
            out_handle.write("".join(header))
    return out_file

def _run_cortex_on_region(region, align_bam, ref_file, work_dir, out_file_base, config):
    """Run cortex on a specified chromosome start/end region.
    """
    kmers = [31, 51, 71]
    min_reads = 1750
    cortex_dir = config_utils.get_program("cortex", config, "dir")
    stampy_dir = config_utils.get_program("stampy", config, "dir")
    vcftools_dir = config_utils.get_program("vcftools", config, "dir")
    if cortex_dir is None or stampy_dir is None:
        raise ValueError("cortex_var requires path to pre-built cortex and stampy")
    region_str = "{0}-{1}-{2}".format(*region)
    base_dir = safe_makedir(os.path.join(work_dir, region_str))
    try:
        out_vcf_base = os.path.join(base_dir, "{0}-{1}".format(
                    os.path.splitext(os.path.basename(out_file_base))[0], region_str))
        out_file = os.path.join(work_dir, os.path.basename("{0}.vcf".format(out_vcf_base)))
        if not file_exists(out_file):
            fastq = _get_fastq_in_region(region, align_bam, out_vcf_base)
            if _count_fastq_reads(fastq, min_reads) < min_reads:
                vcfutils.write_empty_vcf(out_file)
            else:
                local_ref, genome_size = _get_local_ref(region, ref_file, out_vcf_base)
                indexes = _index_local_ref(local_ref, cortex_dir, stampy_dir, kmers)
                cortex_out = _run_cortex(fastq, indexes, {"kmers": kmers, "genome_size": genome_size,
                                                          "sample": get_sample_name(align_bam)},
                                         out_vcf_base, {"cortex": cortex_dir, "stampy": stampy_dir,
                                                        "vcftools": vcftools_dir},
                                         config)
                if cortex_out:
                    _remap_cortex_out(cortex_out, region, out_file)
                else:
                    vcfutils.write_empty_vcf(out_file)
    finally:
        if os.path.exists(base_dir):
            shutil.rmtree(base_dir)
    return [region[0], int(region[1]), int(region[2]), out_file]

def _remap_cortex_out(cortex_out, region, out_file):
    """Remap coordinates in local cortex variant calls to the original global region.
    """
    def _remap_vcf_line(line, contig, start):
        parts = line.split("\t")
        if parts[0] == "" or parts[1] == "":
            return None
        parts[0] = contig
        try:
            parts[1] = str(int(parts[1]) + start)
        except ValueError:
            raise ValueError("Problem in {0} with \n{1}".format(
                    cortex_out, parts))
        return "\t".join(parts)
    def _not_filtered(line):
        parts = line.split("\t")
        return parts[6] == "PASS"
    contig, start, _ = region
    start = int(start)
    with open(cortex_out) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                if line.startswith("##fileDate"):
                    pass
                elif line.startswith("#"):
                    out_handle.write(line)
                elif _not_filtered(line):
                    update_line = _remap_vcf_line(line, contig, start)
                    if update_line:
                        out_handle.write(update_line)

def _run_cortex(fastq, indexes, params, out_base, dirs, config):
    """Run cortex_var run_calls.pl, producing a VCF variant file.
    """
    print(out_base)
    fastaq_index = "{0}.fastaq_index".format(out_base)
    se_fastq_index = "{0}.se_fastq".format(out_base)
    pe_fastq_index = "{0}.pe_fastq".format(out_base)
    reffasta_index = "{0}.list_ref_fasta".format(out_base)
    with open(se_fastq_index, "w") as out_handle:
        out_handle.write(fastq + "\n")
    with open(pe_fastq_index, "w") as out_handle:
        out_handle.write("")
    with open(fastaq_index, "w") as out_handle:
        out_handle.write("{0}\t{1}\t{2}\t{2}\n".format(params["sample"], se_fastq_index,
                                                       pe_fastq_index))
    with open(reffasta_index, "w") as out_handle:
        for x in indexes["fasta"]:
            out_handle.write(x + "\n")
    os.environ["PERL5LIB"] = "{0}:{1}:{2}".format(
        os.path.join(dirs["cortex"], "scripts/calling"),
        os.path.join(dirs["cortex"], "scripts/analyse_variants/bioinf-perl/lib"),
        os.environ.get("PERL5LIB", ""))
    kmers = sorted(params["kmers"])
    kmer_info = ["--first_kmer", str(kmers[0])]
    if len(kmers) > 1:
        kmer_info += ["--last_kmer", str(kmers[-1]),
                      "--kmer_step", str(kmers[1] - kmers[0])]
    subprocess.check_call(["perl", os.path.join(dirs["cortex"], "scripts", "calling", "run_calls.pl"),
                           "--fastaq_index", fastaq_index,
                           "--auto_cleaning", "yes", "--bc", "yes", "--pd", "yes",
                           "--outdir", os.path.dirname(out_base), "--outvcf", os.path.basename(out_base),
                           "--ploidy", str(config["algorithm"].get("ploidy", 2)),
                           "--stampy_hash", indexes["stampy"],
                           "--stampy_bin", os.path.join(dirs["stampy"], "stampy.py"),
                           "--refbindir", os.path.dirname(indexes["cortex"][0]),
                           "--list_ref_fasta",  reffasta_index,
                           "--genome_size", str(params["genome_size"]),
                           "--max_read_len", "30000",
                           #"--max_var_len", "4000",
                           "--format", "FASTQ", "--qthresh", "5", "--do_union", "yes",
                           "--mem_height", "17", "--mem_width", "100",
                           "--ref", "CoordinatesAndInCalling", "--workflow", "independent",
                           "--vcftools_dir", dirs["vcftools"],
                           "--logfile", "{0}.logfile,f".format(out_base)]
                          + kmer_info)
    final = glob.glob(os.path.join(os.path.dirname(out_base), "vcfs",
                                   "{0}*FINALcombined_BC*decomp.vcf".format(os.path.basename(out_base))))
    # No calls, need to setup an empty file
    if len(final) != 1:
        print("Did not find output VCF file for {0}".format(out_base))
        return None
    else:
        return final[0]

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
                                   "--max_read_len", "30000",
			           "--sample_id", base_out,
                                   "--dump_binary", out_file])
        cindexes.append(out_file)
    if not file_exists("{0}.stidx".format(base_out)):
        subprocess.check_call([os.path.join(stampy_dir, "stampy.py"), "-G",
                               base_out, fasta_file])
        subprocess.check_call([os.path.join(stampy_dir, "stampy.py"), "-g",
                               base_out, "-H", base_out])
    return {"stampy": base_out,
            "cortex": cindexes,
            "fasta": [fasta_file]}

def _get_local_ref(region, ref_file, out_vcf_base):
    """Retrieve a local FASTA file corresponding to the specified region.
    """
    out_file = "{0}.fa".format(out_vcf_base)
    if not file_exists(out_file):
        with pysam.Fastafile(ref_file) as in_pysam:
            contig, start, end = region
            seq = in_pysam.fetch(contig, int(start), int(end))
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
        with pysam.Samfile(align_bam, "rb") as in_pysam:
            with file_transaction(out_file) as tx_out_file:
                with open(tx_out_file, "w") as out_handle:
                    contig, start, end = region
                    for read in in_pysam.fetch(contig, int(start), int(end)):
                        seq = Seq.Seq(read.seq)
                        qual = list(read.qual)
                        if read.is_reverse:
                            seq = seq.reverse_complement()
                            qual.reverse()
                        out_handle.write("@{name}\n{seq}\n+\n{qual}\n".format(
                                name=read.qname, seq=str(seq), qual="".join(qual)))
    return out_file

## Utility functions

def _count_fastq_reads(in_fastq, min_reads):
    """Count the number of fastq reads in a file, stopping after reaching min_reads.
    """
    with open(in_fastq) as in_handle:
        items = list(itertools.takewhile(lambda i : i <= min_reads,
                                         (i for i, _ in enumerate(FastqGeneralIterator(in_handle)))))
    return len(items)

def get_sample_name(align_bam):
    with pysam.Samfile(align_bam, "rb") as in_pysam:
        if "RG" in in_pysam.header:
            return in_pysam.header["RG"][0]["SM"]
