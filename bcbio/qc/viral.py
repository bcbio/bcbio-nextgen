"""Detect viral infections via bwa alignment of unaligned reads.

This is primarily useful for cancer samples where viral infection can
inform treatment.
"""
import glob
import os

from bcbio import bam, utils
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.variation import vcfutils

def run(bam_file, data, out_dir):
    """Run viral QC analysis:
       1. Extract the unmapped reads
       2. BWA-MEM to the viral sequences from GDC database https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
       3. Report viruses that are in more than 50% covered by at least 5x
    """
    source_link = 'https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files'
    viral_target = "gdc-viral"
    out = {}
    if vcfutils.get_paired_phenotype(data):
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
                with file_transaction(data, out_file) as tx_out_file:
                    sample_name = dd.get_sample_name(data)
                    mosdepth_prefix = os.path.splitext(viral_bam)[0]
                    cmd = ("mosdepth -t {cores} {mosdepth_prefix} {viral_bam} -n --thresholds 1,5,25 --by "
                           "<(awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}' {viral_ref}.fai) && "
                           "echo '## Viral sequences (from {source_link}) found in unmapped reads' > {tx_out_file} &&"
                           "echo '## Sample: {sample_name}' >> {tx_out_file} && "
                           "echo '#virus\tsize\tdepth\t1x\t5x\t25x' >> {tx_out_file} && "
                           "paste <(zcat {mosdepth_prefix}.regions.bed.gz) <(zgrep -v ^# {mosdepth_prefix}.thresholds.bed.gz) | "
                           "awk 'BEGIN {{FS=\"\\t\"}} {{ print $1 FS $3 FS $4 FS $10/$3 FS $11/$3 FS $12/$3}}' | "
                           "sort -n -r -k 5,5 >> {tx_out_file}")
                    do.run(cmd.format(**locals()), "Analyse coverage of viral genomes")
            out["base"] = out_file
            out["secondary"] = []
    return out

def get_files(data):
    """Retrieve pre-installed viral reference files.
    """
    all_files = glob.glob(os.path.normpath(os.path.join(os.path.dirname(dd.get_ref_file(data)),
                                                        os.pardir, "viral", "*")))
    return sorted(all_files)
