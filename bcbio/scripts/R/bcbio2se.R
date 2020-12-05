# create SummarizedExperiment object (se) from bcbio output
# usage:
# Rscript bcbio2se.R work_dir out_file
# output: out_file a RDS file with the SE object in it

library(tidyverse)
library(tximport)
library(SummarizedExperiment)
library(janitor)
library(DESeq2)

args <- commandArgs(trailingOnly = T)
work_dir <- args[1]
out_file <- args[2]

metadata <- read_csv(file.path(work_dir, "metadata.csv"))
colnames(metadata)[1] <- "sample"
metadata$sample = make_clean_names(metadata$sample)

metrics <- read_tsv(file.path(work_dir, "qc",
                  "multiqc",
                  "multiqc_data",
                  "multiqc_bcbio_metrics.txt")) %>%
  clean_names(case = "snake")
salmon_dir = file.path(work_dir, "salmon")
salmon_files = file.path(salmon_dir, metrics$sample, "quant.sf")
metrics$sample = make_clean_names(metrics$sample)
names(salmon_files) = metrics$sample

transcripts2genes_file <- file.path(work_dir, "inputs", "transcriptome", "tx2gene.csv")
transcripts2genes <- read_csv(transcripts2genes_file,
                              col_names = c("transcript_id", "gene_id"))

txi_salmon <- tximport(salmon_files, type = "salmon",
                       tx2gene = transcripts2genes,
                       countsFromAbundance = "lengthScaledTPM",
                       dropInfReps=TRUE)

col_data <- metadata %>%
  column_to_rownames(var = "sample")
col_data$sample <- rownames(col_data)
col_data = col_data[names(salmon_files),]

raw_counts <- round(data.frame(txi_salmon$counts, check.names = FALSE), 0) %>% as.matrix()

se_metadata <- list(metrics = metrics,
                 countsFromAbundance = txi_salmon$countsFromAbundance)

vst <- vst(raw_counts)

se <- SummarizedExperiment(assays = list(raw = raw_counts,
                                         tpm = txi_salmon$abundance,
                                         length = txi_salmon$length,
                                         vst = vst),
                           colData = col_data,
                           metadata = se_metadata)
saveRDS(se, out_file)
