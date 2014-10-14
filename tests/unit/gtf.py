import os
import shutil
import numpy
from bcbio.rnaseq import gtf, cpat, annotate_gtf
from bcbio.utils import safe_makedir, file_exists
from pprint import pprint

# TEST_GTF = "../data/genomes/hg19/rnaseq/ref-transcripts.gtf"
# TEST_FASTA = "../data/genomes/hg19/seq/hg19.fa"
# TEST_GTF = "big_test/assembled.gtf"
#TEST_GTF = "/v-data/incubator/transcriptome-simulation/genomes/Hsapiens/chr22/rnaseq/ref-transcripts.gtf"
#TEST_FASTA = "/v-data/incubator/transcriptome-simulation/genomes/Hsapiens/chr22/seq/chr22.fa"
#TEST_ASSEMBLED = "big_test/assembled.gtf"
#TEST_GTF = "/usr/local/share/bcbio-nextgen/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf"
#TEST_FASTA = "/usr/local/share/bcbio-nextgen/genomes/Hsapiens/GRCh37/seq/GRCH37.fa"
TEST_FASTA = "/v-data/incubator/transcriptome-simulation/genomes/Hsapiens/chr22/seq/chr22.fa"
TEST_GTF = "/v-data/incubator/transcriptome-simulation/genomes/Hsapiens/chr22/rnaseq/ref-transcripts.gtf"
TEST_ASSEMBLED = "../../../incubator/transcriptome-simulation/flux-simulator/half/chr22_half_clean/2014-06-14_half/assembled.gtf"

# def test_get_gtf_db():
#     db = gtf.get_gtf_db(TEST_GTF)
#     assert list(db.all_features()) is not []

# def test_gtf_to_fasta():
#     safe_makedir("test_gtf")
#     out_file = os.path.join("test_gtf", "ref-transcripts.fa")

#     fa = gtf.gtf_to_fasta(TEST_GTF, TEST_FASTA, out_file)
#     assert(file_exists(fa))
#     #shutil.rmtree("test_gtf")

# def test_noncoding_gtf():
#     safe_makedir("test_gtf")
#     noncoding = gtf.partition_gtf(TEST_GTF)
#     assert(file_exists(noncoding))

# def test_coding_gtf():
#     safe_makedir("test_gtf")
#     coding = gtf.partition_gtf(TEST_GTF, coding=True)
#     assert(file_exists(coding))

def test_annotation():
    safe_makedir("test_gtf")
    fixed = annotate_gtf.annotate_novel_coding(TEST_ASSEMBLED, TEST_GTF, TEST_FASTA)
    cleaned = annotate_gtf.remove_noncongruent_transcripts(fixed, TEST_GTF, TEST_FASTA)
    print cleaned


# def test_hexamer_table():
#     safe_makedir("test_gtf")
#     cutoff = cpat.get_coding_potential_cutoff(TEST_GTF, TEST_FASTA)
#     print cutoff

    # train_gtf, test_gtf = gtf.split_gtf(TEST_GTF, out_dir="test_gtf")
    # coding_gtf = gtf.partition_gtf(train_gtf, coding=True)
    # noncoding_gtf = gtf.partition_gtf(train_gtf)
    # noncoding_fasta = gtf.gtf_to_fasta(noncoding_gtf, TEST_FASTA)
    # cds_fasta = gtf.gtf_to_fasta(coding_gtf, TEST_FASTA, cds=True)
    # hexamer_content = gtf.hexamer_table(cds_fasta, noncoding_fasta)
    # coding_fasta = gtf.gtf_to_fasta(coding_gtf, TEST_FASTA)
    # logit_model = gtf.make_logit_model(coding_fasta, noncoding_fasta,
    #                                    hexamer_content, "test_gtf")
    # test_fasta = gtf.gtf_to_fasta(test_gtf, TEST_FASTA)
    # cpat_fn = gtf.cpat(test_fasta, hexamer_content, logit_model)
    # cpat = gtf.load_cpat_coding_prob(cpat_fn)
    # coding, noncoding = gtf.get_coding_noncoding_transcript_ids(test_gtf)
    # best_score = 1
    # best_cutoff = 0
    # best_sensitivity = 0
    # best_specificity = 0
    # for cutoff in list(numpy.arange(0.1, 1, 0.01)):
    #     grade = gtf.grade_cpat(coding, noncoding, cpat, cutoff)
    #     score = abs(grade["sensitivity"] - grade["specificity"])
    #     if score < best_score:
    #         best_score = score
    #         best_cutoff = cutoff
    #         best_sensitivity = grade["sensitivity"]
    #         best_specificity = grade["specificity"]
    # print "cutoff: %f, score: %f" % (best_cutoff, best_score)
    # print "sensitivity: %f, specificity: %f " % (best_sensitivity, best_specificity)
