"""Detect viral infections via bwa alignment of unaligned reads.

This is primarily useful for cancer samples where viral infection can
inform treatment.
"""
import glob
import os
import subprocess

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils
from bcbio.heterogeneity import chromhacks

def run(bam_file, data, out_dir):
    """Run viral QC analysis:
       1. Extract the unmapped reads
       2. BWA-MEM to the viral sequences from GDC database https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
       3. Report viruses that are in more than 50% covered by at least 5x
    """
    source_link = 'https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files'
    viral_target = "gdc-viral"
    out = {}
    viral_refs = [x for x in dd.get_viral_files(data) if os.path.basename(x) == "%s.fa" % viral_target]
    if viral_refs and utils.file_exists(viral_refs[0]):
        viral_ref = viral_refs[0]
        viral_bam = os.path.join(utils.safe_makedir(out_dir),
                                    "%s-%s.bam" % (dd.get_sample_name(data),
                                                utils.splitext_plus(os.path.basename(viral_ref))[0]))
        out_file = "%s-completeness.txt" % utils.splitext_plus(viral_bam)[0]
        cores = dd.get_num_cores(data)
        if not utils.file_uptodate(out_file, bam_file):
            if not utils.file_uptodate(viral_bam, bam_file):
                with file_transaction(data, viral_bam) as tx_out_file:
                    tmpfile = "%s-tmp" % utils.splitext_plus(tx_out_file)[0]
                    cmd = ("samtools view -u -f 4 {bam_file} | "
                            "bamtofastq collate=0 | "
                            "bwa mem -t {cores} {viral_ref} - | "
                            "bamsort tmpfile={tmpfile} inputthreads={cores} outputthreads={cores} "
                            "inputformat=sam index=1 indexfilename={tx_out_file}.bai O={tx_out_file}")
                    do.run(cmd.format(**locals()), "Align unmapped reads to viral genome")

            total_reads = _count_reads(bam_file)
            assert total_reads > 0, 'Reads count is {total_reads}, is there a bug in counting the read count? {bam_file}'.format(**locals())
            with file_transaction(data, out_file) as tx_out_file:
                sample_name = dd.get_sample_name(data)
                mosdepth_prefix = os.path.splitext(viral_bam)[0]
                cmd = ("mosdepth -t {cores} {mosdepth_prefix} {viral_bam} -n --thresholds 1,5,25 --by "
                       "<(awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}' {viral_ref}.fai) && "
                       "echo '## Viral sequences (from {source_link}) found in unmapped reads' > {tx_out_file} &&"
                       "echo '## Sample: {sample_name}' >> {tx_out_file} && "
                       "echo '#virus\tsize\tdepth\t1x\t5x\t25x\treads\treads_pct' >> {tx_out_file} && "
                       "paste "
                       "<(zcat {mosdepth_prefix}.regions.bed.gz) "
                       "<(zgrep -v ^# {mosdepth_prefix}.thresholds.bed.gz) "
                       "<(samtools idxstats {viral_bam} | grep -v '*') | "
                       "awk 'BEGIN {{FS=\"\\t\"}} {{ print $1 FS $3 FS $4 FS $10/$3 FS $11/$3 FS $12/$3 FS $15 FS $15/{total_reads}}}' | "
                       "sort -n -r -k 5,5 >> {tx_out_file}")
                do.run(cmd.format(**locals()), "Analyse coverage of viral genomes")
                if chromhacks.get_EBV(data):
                    ref_file = dd.get_ref_file(data)
                    work_bam = dd.get_work_bam(data)
                    ebv = chromhacks.get_EBV(data)
                    mosdepth_prefix = os.path.splitext(work_bam)[0] + "-EBV"
                    cmd = ("mosdepth -t {cores} {mosdepth_prefix} {work_bam} -n --thresholds 1,5,25 --by "
                            "<(grep {ebv} {ref_file}.fai | awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}') && "
                            "paste "
                            "<(zcat {mosdepth_prefix}.regions.bed.gz) "
                            "<(zgrep -v ^# {mosdepth_prefix}.thresholds.bed.gz) "
                            "<(samtools idxstats {work_bam} | grep {ebv}) | "
                            "awk 'BEGIN {{FS=\"\\t\"}} {{ print $1 FS $3 FS $4 FS $10/$3 FS $11/$3 FS $12/$3 FS $15 FS $15/{total_reads}}}' | "
                            "sort -n -r -k 5,5 >> {tx_out_file}")
                    do.run(cmd.format(**locals()), "Analyse coverage of EBV")

        out["base"] = out_file
        out["secondary"] = []
    return out

def get_files(data):
    """Retrieve pre-installed viral reference files.
    """
    all_files = glob.glob(os.path.normpath(os.path.join(os.path.dirname(dd.get_ref_file(data)),
                                                        os.pardir, "viral", "*")))
    return sorted(all_files)

def _count_reads(bam_file):
    cmd = "samtools idxstats %s | awk '{sum += $3 + $4} END {print sum}'"
    count = subprocess.check_output(cmd % bam_file, shell=True)
    return int(count.strip())
