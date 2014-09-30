import gffutils
from argparse import ArgumentParser
from gffutils.iterators import DataIterator

dialect = {'field separator': '; ',
           'fmt': 'gtf',
           'keyval separator': ' ',
           'leading semicolon': False,
           'multival separator': ',',
           'quoted GFF2 values': True,
           'order': ['gene_id', 'transcript_id'],
           'repeated keys': False,
           'trailing semicolon': True}

if __name__ == "__main__":
    parser = ArgumentParser(description="GFF3->GTF converter for use with bcbio_setup_genome")
    parser.add_argument("GFF3", help="GFF3 file to convert")

    args = parser.parse_args()
    db = gffutils.create_db(args.GFF3, ":memory:")

    for feature in DataIterator(db.features_of_type("exon"), dialect=dialect):
        transcript_id = feature["Parent"][0]
        gene_id = db[transcript_id]["Parent"][0]
        attr = {"transcript_id": transcript_id, "gene_id": gene_id}
        attributes = gffutils.attributes.Attributes(attr)
        feature.attributes = attributes
        print feature
