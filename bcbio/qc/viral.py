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
    """Run viral QC analysis.
    """
    viral_target = "gdc-viral"
    out = {}
    if vcfutils.get_paired_phenotype(data):
        viral_refs = [x for x in dd.get_viral_files(data) if os.path.basename(x) == "%s.fa" % viral_target]
        if viral_refs and utils.file_exists(viral_refs[0]):
            viral_ref = viral_refs[0]
            viral_bam = os.path.join(utils.safe_makedir(out_dir),
                                     "%s-%s.bam" % (dd.get_sample_name(data),
                                                    utils.splitext_plus(os.path.basename(viral_ref))[0]))
            out_file = "%s-counts.txt" % utils.splitext_plus(viral_bam)[0]
            if not utils.file_uptodate(out_file, bam_file):
                if not utils.file_uptodate(viral_bam, bam_file):
                    with file_transaction(data, viral_bam) as tx_out_file:
                        cores = dd.get_num_cores(data)
                        tmpfile = "%s-tmp" % utils.splitext_plus(tx_out_file)[0]
                        cmd = ("samtools view -u -f 4 {bam_file} | "
                               "bamtofastq collate=0 | "
                               "bwa mem -t {cores} {viral_ref} - | "
                               "bamsort tmpfile={tmpfile} inputthreads={cores} outputthreads={cores} "
                               "inputformat=sam index=1 indexfilename={tx_out_file}.bai O={tx_out_file}")
                        do.run(cmd.format(**locals()), "Compare unmapped reads to viral genome")
                with file_transaction(data, out_file) as tx_out_file:
                    with open(tx_out_file, "w") as out_handle:
                        out_handle.write("# sample\t%s\n" % dd.get_sample_name(data))
                        for info in bam.idxstats(viral_bam, data):
                            if info.aligned > 0:
                                out_handle.write("%s\t%s\n" % (info.contig, info.aligned))
            out["base"] = out_file
    return out

def get_files(data):
    """Retrieve pre-installed viral reference files.
    """
    all_files = glob.glob(os.path.normpath(os.path.join(os.path.dirname(dd.get_ref_file(data)),
                                                        os.pardir, "viral", "*")))
    return sorted(all_files)
