"""
run the Coding Potential Assessment Tool (CPAT)
http://nar.oxfordjournals.org/content/early/2013/01/17/nar.gkt006.full
"""
import numpy
import shutil
import tempfile
import os
import sys

from bcbio import utils
from bcbio.rnaseq import gtf
from bcbio.utils import file_exists, safe_makedir
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.bam import fasta

def classify_with_cpat(assembled_gtf, ref_gtf, ref_fasta):
    cpat_cmd = _find_executable("cpat.py")
    if not cpat_cmd:
        return {}
    cutoff, hexamer, logit = get_coding_potential_cutoff(ref_gtf, ref_fasta)
    assembled_fasta = gtf.gtf_to_fasta(assembled_gtf, ref_fasta)
    cpat_fn = cpat(assembled_fasta, hexamer, logit)
    coding_probabilities = load_cpat_coding_prob(cpat_fn)
    lengths = fasta.sequence_length(assembled_fasta)
    classification = {}
    for transcript, prob in coding_probabilities.items():
        if prob > cutoff:
            classification[transcript] = "protein_coding"
        if lengths[transcript] > 200:
            classification[transcript] = "lncRNA"
        else:
            classification[transcript] = "ncRNA"
    return classification

def _find_executable(name):
    in_path = utils.which(name)
    if in_path:
        return in_path
    else:
        in_conda = os.path.join(os.path.dirname(sys.executable), name)
        if os.path.exists(in_conda):
            return in_conda
        else:
            return None

def cpat(assembled_fasta, hexamer, logit, out_file=None):
    if out_file and file_exists(out_file):
        return out_file
    if not out_file:
        out_file = tempfile.NamedTemporaryFile(delete=False, suffix=".cpat").name

    cpat_cmd = _find_executable("cpat.py")
    cmd = ("{cpat_cmd} --gene={assembled_fasta} --hex={hexamer} "
           "--logitModel={logit} --outfile={tx_out_file}")
    message = "Predicing coding potential of %s." % (assembled_fasta)
    with file_transaction(out_file) as tx_out_file:
        do.run(cmd.format(**locals()), message)
    return out_file

def load_cpat_coding_prob(cpat_file):
    with open(cpat_file) as in_handle:
        header = in_handle.next()
        return {line.split()[0]: float(line.split()[5]) for line in in_handle}

def load_cpat_orf_size(cpat_file):
    with open(cpat_file) as in_handle:
        header = in_handle.next()
        return {line.split()[0]: float(line.split()[2]) for line in in_handle}

def grade_cpat(coding_transcripts, noncoding_transcripts, cpat, cutoff):
    coding_tp = 0
    coding_fp = 0
    noncoding_tp = 0
    noncoding_fp = 0
    for transcript in coding_transcripts:
        if cpat[transcript] < cutoff:
            noncoding_fp += 1
        else:
            coding_tp += 1
    for transcript in noncoding_transcripts:
        if cpat[transcript] >= cutoff:
            coding_fp += 1
        else:
            noncoding_tp += 1
    tp = float(coding_tp)
    fp = float(coding_fp)
    tn = float(noncoding_tp)
    fn = float(noncoding_fp)
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision = tp / (tp + fp) if (tp + fp > 0) else -1
    return {"sensitivity": sensitivity, "specificity": specificity,
            "accuracy": accuracy, "precision": precision}

def make_logit_model(coding_fasta, noncoding_fasta, hexamers, out_dir=None):
    safe_makedir(out_dir)
    out_prefix = os.path.join(out_dir, "logit")
    out_file = out_prefix + ".logit.RData"
    if file_exists(out_file):
        return out_file

    tx_prefix = tempfile.NamedTemporaryFile(delete=False).name
    tx_out_file = tx_prefix +  ".logit.RData"

    logit_cmd = _find_executable("make_logitModel.py")
    cmd = ("{logit_cmd} --cgene={coding_fasta} --ngene={noncoding_fasta} "
           "--hex={hexamers} --outfile={tx_prefix}")
    message = "Building coding/noncoding logistical model."
    do.run(cmd.format(**locals()), message)

    shutil.move(tx_out_file, out_file)

    return out_file

def get_coding_potential_cutoff(ref_gtf, ref_fasta):
    """
    estimate the coding potential cutoff that best classifies
    coding/noncoding transcripts by splitting the reference
    annotation into a test and training set and determining
    the cutoff where the sensitivity and specificity meet
    """
    train_gtf, test_gtf = gtf.split_gtf(ref_gtf, sample_size=2000)
    coding_gtf = gtf.partition_gtf(train_gtf, coding=True)
    noncoding_gtf = gtf.partition_gtf(train_gtf)
    noncoding_fasta = gtf.gtf_to_fasta(noncoding_gtf, ref_fasta)
    cds_fasta = gtf.gtf_to_fasta(coding_gtf, ref_fasta, cds=True)
    hexamer_content = hexamer_table(cds_fasta, noncoding_fasta)
    coding_fasta = gtf.gtf_to_fasta(coding_gtf, ref_fasta)
    logit_model = make_logit_model(coding_fasta, noncoding_fasta,
                                       hexamer_content, "test_gtf")
    test_fasta = gtf.gtf_to_fasta(test_gtf, ref_fasta)
    cpat_fn = cpat(test_fasta, hexamer_content, logit_model)
    cpat_prob = load_cpat_coding_prob(cpat_fn)
    coding, noncoding = gtf.get_coding_noncoding_transcript_ids(test_gtf)
    best_score = 1
    best_cutoff = 0
    best_sensitivity = 0
    best_specificity = 0
    for cutoff in list(numpy.arange(0.1, 1, 0.01)):
        grade = grade_cpat(coding, noncoding, cpat_prob, cutoff)
        score = abs(grade["sensitivity"] - grade["specificity"])
        if score < best_score:
            best_score = score
            best_cutoff = cutoff
            best_sensitivity = grade["sensitivity"]
            best_specificity = grade["specificity"]
    return best_cutoff, hexamer_content, logit_model

def hexamer_table(cds_fasta, noncoding_fasta, out_file=None):
    if out_file and file_exists(out_file):
        return out_file
    if not out_file:
        out_file = tempfile.NamedTemporaryFile(delete=False, suffix=".hexamers").name
    hex_cmd = _find_executable("make_hexamer_tab.py")
    cmd = ("{hex_cmd} --cod={cds_fasta} --noncod={noncoding_fasta} "
           "> {tx_out_file}")
    with file_transaction(out_file) as tx_out_file:
        message = ("Calculating hexamer content in %s and %s."
                   % (cds_fasta, noncoding_fasta))
        do.run(cmd.format(**locals()), message)
    return out_file
