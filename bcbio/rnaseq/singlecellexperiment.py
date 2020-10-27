from __future__ import print_function

import os
import subprocess

from bcbio.rnaseq import gtf
from bcbio.utils import file_exists, Rscript_cmd, R_sitelib
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.pipeline import datadict as dd
from bcbio.log import logger

def make_scrnaseq_object(samples):
    """
    load the initial se.rda object using sinclecell-experiment
    """
    local_sitelib = R_sitelib()
    counts_dir = os.path.dirname(dd.get_in_samples(samples, dd.get_combined_counts))
    gtf_file = dd.get_in_samples(samples, dd.get_transcriptome_gtf)
    if not gtf_file:
        gtf_file = dd.get_in_samples(samples, dd.get_gtf_file)
    rda_file = os.path.join(counts_dir, "se.rda")
    if not file_exists(rda_file):
        with file_transaction(rda_file) as tx_out_file:
            rcode = "%s-run.R" % os.path.splitext(rda_file)[0]
            rrna_file = "%s-rrna.txt" % os.path.splitext(rda_file)[0]
            rrna_file = _find_rRNA_genes(gtf_file, rrna_file)
            with open(rcode, "w") as out_handle:
                out_handle.write(_script.format(**locals()))
            rscript = Rscript_cmd()
            try:
                # do.run([rscript, "--vanilla", rcode],
                #        "SingleCellExperiment",
                #        log_error=False)
                rda_file = rcode
            except subprocess.CalledProcessError as msg:
                logger.exception()


def _find_rRNA_genes(gtf_file, rrna_file):
    rrna_features = gtf.get_rRNA(gtf_file)
    transcripts = set([x[0] for x in rrna_features if x])
    with open(rrna_file, 'w') as outh:
        outh.write("\n".join(transcripts))
    return rrna_file


_script = """
library(SingleCellExperiment)
library(Matrix)

counts = readMM("{counts_dir}/tagcounts.mtx")
rownames = read.csv("{counts_dir}/tagcounts.mtx.rownames", header = F)[["V1"]]
rownames = as.character(rownames)
colnames = read.csv("{counts_dir}/tagcounts.mtx.colnames", header = F)[["V1"]]
colnames = make.names(as.character(colnames))
reads = read.csv("{counts_dir}/cb-histogram.txt", header = F, sep="\t", row.names = 1)
rownames(reads) = make.names(rownames(reads))

counts =  as(counts, "dgCMatrix")
rownames(counts) = rownames
colnames(counts) = colnames
metadata = read.csv("{counts_dir}/tagcounts.mtx.metadata")
rownames(metadata) = colnames
metadata[["nUMI"]] = colSums(counts)
metadata[["nGenes"]] = colSums(counts>0)
metadata[["log10GenesPerUMI"]] = log10(metadata$nGene) / log10(metadata$nUMI)
metadata[["nReads"]] = reads[colnames,]

rrna = read.csv("{rrna_file}", header=F, stringsAsFactors = F)
metadata[["mtUMI"]] = colSums(counts[rrna[["V1"]],], na.rm = T)
metadata[["mtUMI"]][is.na(metadata[["mtUMI"]])] = 0
metadata[["mitoRatio"]] = metadata$mtUMI/metadata$nUMI

se = SingleCellExperiment(assays=list(raw=counts), colData = metadata)
save(se, file = "{counts_dir}/se.rda")
"""
