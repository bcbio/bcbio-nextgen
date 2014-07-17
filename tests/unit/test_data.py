import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "bcbio-nextgen-test-data", "data")
COUNT_DIR = os.path.join(DATA_DIR, "count")
MOUSE_DIR = os.path.join(DATA_DIR, "organisms", "mouse")

GTF_FILE = os.path.join(MOUSE_DIR, "mouse.gtf")
DEXSEQ_GFF = os.path.join(MOUSE_DIR, "mouse.dexseq.gff3")
DEXSEQ_COUNT_FILES = [os.path.join(COUNT_DIR, "dexseq", "dexseq-counts.txt"),
                      os.path.join(COUNT_DIR, "dexseq", "dexseq-counts2.txt")]
BAM_FILE = os.path.join(MOUSE_DIR, "mouse.bam")
BED_FILE = os.path.join(MOUSE_DIR, "mouse.bed")
