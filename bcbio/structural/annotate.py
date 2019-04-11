"""Annotate structural variant calls with associated genes.
"""
import os
import sys

from bcbio import utils
from bcbio.bam import ref
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline import datadict as dd
from bcbio.provenance import do
from bcbio.structural import regions
from bcbio.variation import bedutils
from bcbio.pipeline import config_utils

import pybedtools

def add_genes(in_file, data, max_distance=10000, work_dir=None):
    """Add gene annotations to a BED file from pre-prepared RNA-seq data.

    max_distance -- only keep annotations within this distance of event
    """
    gene_file = regions.get_sv_bed(data, "exons", out_dir=os.path.dirname(in_file))
    if gene_file and utils.file_exists(in_file):
        out_file = "%s-annotated.bed" % utils.splitext_plus(in_file)[0]
        if work_dir:
            out_file = os.path.join(work_dir, os.path.basename(out_file))
        if not utils.file_uptodate(out_file, in_file):
            fai_file = ref.fasta_idx(dd.get_ref_file(data))
            with file_transaction(data, out_file) as tx_out_file:
                _add_genes_to_bed(in_file, gene_file, fai_file, tx_out_file, data, max_distance)
        return out_file
    else:
        return in_file

def _add_genes_to_bed(in_file, gene_file, fai_file, out_file, data, max_distance=10000):
    """Re-usable subcomponent that annotates BED file genes from another BED
    """
    try:
        input_rec = next(iter(pybedtools.BedTool(in_file)))
    except StopIteration:  # empty file
        utils.copy_plus(in_file, out_file)
        return
    # keep everything after standard chrom/start/end, 1-based
    extra_fields = list(range(4, len(input_rec.fields) + 1))
    # keep the new gene annotation
    gene_index = len(input_rec.fields) + 4
    extra_fields.append(gene_index)
    columns = ",".join([str(x) for x in extra_fields])
    max_column = max(extra_fields) + 1
    ops = ",".join(["distinct"] * len(extra_fields))
    # swap over gene name to '.' if beyond maximum distance
    # cut removes the last distance column which can cause issues
    # with bedtools merge: 'ERROR: illegal character '.' found in integer conversion of string'
    distance_filter = (r"""awk -F$'\t' -v OFS='\t' '{if ($NF > %s || $NF < -%s) $%s = "."} {print}'""" %
                       (max_distance, max_distance, gene_index))
    sort_cmd = bedutils.get_sort_cmd(os.path.dirname(out_file))
    cat_cmd = "zcat" if in_file.endswith(".gz") else "cat"
    # Ensure gene transcripts match reference genome
    ready_gene_file = os.path.join(os.path.dirname(out_file), "%s-genomeonly.bed" %
                                   (utils.splitext_plus(os.path.basename(gene_file))[0]))
    ready_gene_file = bedutils.subset_to_genome(gene_file, ready_gene_file, data)
    exports = "export TMPDIR=%s && %s" % (os.path.dirname(out_file), utils.local_path_export())
    bcbio_py = sys.executable
    gsort = config_utils.get_program("gsort", data)
    cmd = ("{exports}{cat_cmd} {in_file} | grep -v ^track | grep -v ^browser | grep -v ^# | "
           "{bcbio_py} -c 'from bcbio.variation import bedutils; bedutils.remove_bad()' | "
           "{gsort} - {fai_file} | "
            "bedtools closest -g {fai_file} "
            "-D ref -t first -a - -b <({gsort} {ready_gene_file} {fai_file}) | "
            "{distance_filter} | cut -f 1-{max_column} | "
            "bedtools merge -i - -c {columns} -o {ops} -delim ',' -d -10 > {out_file}")
    do.run(cmd.format(**locals()), "Annotate BED file with gene info")

def gene_one_per_line(in_file, data):
    """Split comma-separated gene annotations (after add_genes). Leads to duplicated records.
       Input:
          chr1 100 200 F1,F2
       Output:
          chr1 100 200 F1
          chr1 100 200 F2
    """
    if in_file:
        # Report all duplicated annotations one-per-line
        one_per_line_file = "%s-opl.bed" % utils.splitext_plus(in_file)[0]
        if not utils.file_uptodate(one_per_line_file, in_file):
            with file_transaction(data, one_per_line_file) as tx_out_file:
                with open(tx_out_file, 'w') as out:
                    for r in pybedtools.BedTool(in_file):
                        for g in r.name.split(','):
                            out.write('\t'.join(map(str, [r.chrom, r.start, r.end, g])) + '\n')
        return one_per_line_file

def count_genes(in_file, data):
    if pybedtools.BedTool(in_file).field_count() <= 3:
        ann_bed = add_genes(in_file, data)
        ann_bed = gene_one_per_line(ann_bed, data)
    else:
        ann_bed = in_file
    if ann_bed:
        return len(list(set(r.name for r in pybedtools.BedTool(ann_bed)
            if r.name and r.name != ".")))
