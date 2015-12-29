import os

from bcbio.rnaseq import gtf, cpat
from bcbio.bam import fasta
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger

def annotate_novel_coding(assembled_gtf, ref_gtf, ref_fasta, data, out_file=None):
    if not out_file:
        out_file = os.path.splitext(assembled_gtf)[0] + ".annotated.gtf"
    if file_exists(out_file):
        return out_file
    classification = cpat.classify_with_cpat(assembled_gtf, ref_gtf,
                                             ref_fasta, data)
    if not classification:
        logger.info("Protein coding classification of %s was skipped because "
                    "CPAT was not found." % assembled_gtf)
        return assembled_gtf
    ref_db = gtf.get_gtf_db(ref_gtf)
    known_transcript = {feature['transcript_id'][0]: feature.source for feature in
                        gtf.complete_features(ref_db)}
    assembled_db = gtf.get_gtf_db(assembled_gtf)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for feature in gtf.complete_features(assembled_db):
                transcript_id = feature['transcript_id'][0]
                if transcript_id not in known_transcript:
                    feature.source = classification[transcript_id]
                else:
                    feature.source = known_transcript[transcript_id]
                out_handle.write(str(feature) + "\n")
    return out_file

def cleanup_transcripts(assembled_gtf, ref_gtf, ref_fasta, out_file=None):
    """
    Clean up a GTF file of assembled transcripts
    1) if a known gene is known to code for a protein, remove any *novel*
    isoforms of the that do not also code for a protein.
    2) if a new gene has been annotated and none of its isoforms are protein
    coding and it is > 200 bp, mark it as a lincRNA. < 200 bp mark it as ncRNA
    """

    if not out_file:
        out_file = os.path.splitext(assembled_gtf)[0] + ".cleaned.gtf"
    if file_exists(out_file):
        return out_file
    ref_db = gtf.get_gtf_db(ref_gtf)
    known_transcript = {feature['transcript_id'][0]: feature.source for feature
                        in gtf.complete_features(ref_db)}
    ref_gene_to_source = gtf.get_gene_source_set(ref_gtf)
    assembled_db = gtf.get_gtf_db(assembled_gtf)
    assembled_fasta = gtf.gtf_to_fasta(assembled_gtf, ref_fasta)
    lengths = fasta.sequence_length(assembled_fasta)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as out_handle:
            for feature in gtf.complete_features(assembled_db):
                transcript_id = feature['transcript_id'][0]
                gene_id = feature['gene_id'][0]
                if transcript_id in known_transcript:
                    out_handle.write(str(feature) + "\n")
                    continue
                known_coding = "protein_coding" in ref_gene_to_source.get(gene_id, [None])
                if known_coding and feature.source != "protein_coding":
                    continue
                if feature.source != "protein_coding":
                    if lengths[transcript_id] > 200:
                        feature.source = "lincRNA"
                    else:
                        feature.source = "ncRNA"
                out_handle.write(str(feature) + "\n")
    return out_file
