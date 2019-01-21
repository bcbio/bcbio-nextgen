#!/usr/bin/env python
"""
Convert genbank to gtf.

USAGE:
gb2genome.py --gbk file.gbk --prefix ecoli_k12

NOTE:
It's designed to work with gb files coming from GenBank. gene is used as gene_id and transcript_id (locus_tag if gene not present).
Only entries having types in allowedTypes = ['gene','CDS','tRNA','tmRNA','rRNA','ncRNA'] are stored in GTF. Need to include exon processing.
No frame info is processed. Need to be included in order to process genes having introns!

AUTHOR:
Leszek Pryszcz
lpryszcz@crg.eu

MAINTAINERS
Andreas Sjodin
Lorena Pantano

Version 0.2

"""

import sys
from datetime import datetime
from Bio import SeqIO
from argparse import ArgumentParser


def gb2gtf( gbk_fn, gtf_fn, source='gb2gtf', allowedTypes=set(['gene','CDS','tRNA','tmRNA','rRNA','ncRNA']) ):
     with open(gtf_fn, 'w') as out_handle:
         for gb in SeqIO.parse( gbk_fn,'gb' ):
            acc     = gb.id #gb.name #gb.description # # 
            skipped = 0
            skippedTypes = set()
            for f in gb.features:

              #process only gene and CDS entries
              if f.type not in allowedTypes:
                skipped += 1
                skippedTypes.add( f.type )
                continue

              #Extract gene id and transcript id to generate comments field 
              if 'locus_tag' in f.qualifiers:
                #use locul tag as gene_id/transcript_id
                gene_id = f.qualifiers['locus_tag'][0].replace(" ", ".").replace("'","")
                transcript_id = 'T' + gene_id
              elif 'protein_id' in f.qualifiers:
                gene_id = f.qualifiers['protein_id'][0].replace(" ", ".").replace("'","")
                transcript_id = 'T' + gene_id
              elif 'gene' in f.qualifiers:
                gene_id = f.qualifiers['gene'][0].replace(" ", ".").replace("'","")
                transcript_id = 'T' + gene_id
              elif 'label' in f.qualifiers:
                gene_id = f.qualifiers['label'][0].replace(" ", ".").replace("'","")
                transcript_id = 'T' + gene_id

              #Extract gene name and transcript name to generate comments field 
              comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )

              if 'gene' in f.qualifiers:
                clean_name = f.qualifiers['gene'][0].replace(" ", ".").replace("'","")
                comments += '; gene_name "%s"' % clean_name
                comments += '; transcript_name "%s-1"' % clean_name
              elif 'label' in f.qualifiers:
                clean_name = f.qualifiers['label'][0].replace(" ", ".").replace("'","")
                comments += '; gene_name "%s"' % clean_name
                comments += '; transcript_name "%s-1"' % clean_name
              if 'protein_id' in f.qualifiers:
                clean_name = f.qualifiers['protein_id'][0].replace(" ", ".").replace("'","")
                comments += '; protein_id "%s"' % clean_name

              #code strand as +/- (in genbank 1 or -1)
              if int(f.strand)>0: strand = '+'
              else:               strand = '-'

              #Define source
              if f.type == 'CDS':
                source = 'protein_coding'
                f.type = 'exon'
              if f.type == 'tRNA':
                source = 'tRNA'
              if f.type == 'tmRNA':
                source = 'tmRNA'
              if f.type == 'rRNA':
                source = 'rRNA'
              if f.type == 'ncRNA':
                source = 'ncRNA'
              #define gb
              """
              seqname - The name of the sequence. Must be a chromosome or scaffold.
              source - The program that generated this feature.
              feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
              start - The starting position of the feature in the sequence. The first base is numbered 1.
              end - The ending position of the feature (inclusive).
              score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
              strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
              frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
              comments - gene_id "Em:U62317.C22.6.mRNA"; transcript_id "Em:U62317.C22.6.mRNA"; exon_number 1
              """

              if f.type == 'gene':
                continue
              elif f.type == 'CDS':
                continue
              else:
                comments += '; gene_biotype "%s"' % source
                gtf_gene = '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s;' % ( acc, source, 'gene', f.location.start.position+1, f.location.end.position, strand, comments )
                print >>out_handle, gtf_gene
                gtf_tx = '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s;' % ( acc, source, 'transcript', f.location.start.position+1, f.location.end.position, strand, comments )
                print >>out_handle, gtf_tx
                comments += '; exon_number "1"'
                gtf_exon = '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s;' % ( acc, source, 'exon', f.location.start.position+1, f.location.end.position, strand, comments )
                print >>out_handle, gtf_exon

            sys.stderr.write( "%s\tSkipped %s entries having types: %s.\n" % ( gb.id,skipped, ', '.join(skippedTypes) ) )

if __name__=='__main__':
    description = ("Convert GeneBank files to FASTA and GTF bcbio ready files.")

    parser = ArgumentParser(description=description)
    parser.add_argument("--gbk", required=True, help="GBK files")
    parser.add_argument("--prefix", required=True, help="prefix")
    args = parser.parse_args()

    t0=datetime.now()

    count = SeqIO.convert(args.gbk, "genbank", args.prefix + "_tmp.fa", "fasta")
    with open(args.prefix + ".fa", "w") as out_handle:
        with open(args.prefix + "_tmp.fa") as in_handle:
            header = next(in_handle)
            print >>out_handle, header.split()[0]
            for line in in_handle:
                print >>out_handle, line.strip()

    gb2gtf(args.gbk, args.prefix + ".gtf")
    dt=datetime.now() - t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
