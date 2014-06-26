import gffutils
from bcbio.utils import file_exists

def get_gtf_db(gtf, in_memory=False):
    db_file = ":memory:" if in_memory else gtf + ".db"
    if in_memory or not file_exists(db_file):
        db = gffutils.create_db(gtf, dbfn=db_file)
    if in_memory:
        return db
    else:
        return gffutils.FeatureDB(db_file)


