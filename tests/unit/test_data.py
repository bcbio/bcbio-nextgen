import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "bcbio-nextgen-test-data", "data")
MOUSE_DIR = os.path.join(DATA_DIR, "organisms", "mouse")

GTF_FILE = os.path.join(MOUSE_DIR, "mouse.gtf")
DEXSEQ_GFF = os.path.join(MOUSE_DIR, "mouse.dexseq.gff3")
BAM_FILE = os.path.join(MOUSE_DIR, "mouse.bam")
BED_FILE = os.path.join(MOUSE_DIR, "mouse.bed")
