#!/bin/bash


# Reverse engineered what is in the genomes/Organism/build/rnaseq to make it easier to add new genomes 
# and/or get and build updated versions to bcbio-nextgen

genomegrab=$1
ensmble_name=$2

#these passed at update --genome, the ensemble_name may need work but this is a quick hack
#genomegrab=rn5
#ensmble_name=rattus_norvegicus

# Get latest GTF from ensembl

NEWGTF=$(curl -s "ftp://ftp.ensembl.org/pub/current_gtf/$ensmble_name/" --list-only | grep gtf.gz)

echo "Found $NEWGTF file! downloading to ref-transcripts.gtf"
curl -s  "ftp://ftp.ensembl.org/pub/current_gtf/$ensmble_name/$NEWGTF" | gunzip - > ref-transcripts.gtf

# get ensembl bed from UCSC - @jpeden1 loves my pipes and such, just trying to match what I see in the bcbio files
#   no...
#   mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D -N $genomegrab -e 'select * from refFlat'
#   This doesnt match whatever the refFlat in bcbio-next is, so I cut-n-awk

echo "Downloading ensembl genes for $genomegrab in to ref-transcripts.refFlat"
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N  -A -D $genomegrab -e 'select * from ensGene'  | cut -f2- | awk '{print$1"\t"$0}' > ref-transcripts.refFlat

# make the ref-transcripts.genePred, right?
# needs gtfToGenePred binary from kent lab - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

echo "make the ref-transcripts.genePred with Kent's gtfToGenePred"
gtfToGenePred  -genePredExt $NEWGTF  ref-transcripts.genePred

## then for the mask I assume we do something like:

echo "Make the ref-transcripts-mask.gtf and rRNA.gtf files with a little awk love"
awk '$1 ~ /MT/ || $2 ~ /tRNA/ || $2 ~ /rRNA/ || $2 ~ /misc_RNA/' ref-transcripts.gtf > ref-transcripts-mask.gtf
awk '$2 ~ /rRNA/' ref-transcripts.gtf > rRNA.gtf

# now get crazy and make the interval_list file!

echo "Making the rRNA.interval_list with a wicked rad oneliner"
paste <(cut -f1,4,5,7 ref-transcripts-mask.gtf ) <(grep -oP 'transcript_id\s\"*+.*?\d?\"' rRNA.gtf | sed -e 's/\"//g' -e 's/transcript_id //g') |  cat ../seq/$genomegrab.dict - > rRNA.interval_list

#need tophat 2.10.0 to do this without reads... So I'll make one

echo "bcbio has us at topat 2.0.9 so we have to provide a read, so we make just one"

echo "@BCFLAMETRHOWER:183563:AWESOMEFLOWCELL:1:1102:10813:163728 1:N:0:" > read1.fq
echo -e "GTTGGTGGTGTGAGTGTGTTTGTGGTGTGTGTGTGTGTGTGTGAGTGTGTTTGTGGTGTGTGTGTGCGTGTTTGTGGTGTGTGTGTGTTTGTGG\n+" >> read1.fq
echo -e "FEF.F:@A8A>DDDAGCGDFFFCFF<FDFDE>EDECFCFCFFG8GFFBEECFFBDD@EBE<?@B>EEDE8?:@@AFFECEEEDBADDD6DDBB?" >> read1.fq

echo "here is the read 'read1.fq'"
cat read1.fq

echo -e "\n And now we'll build the tophat transcriptome with $NEWGTF for $genomegrab !"
tophat -o /dev/null \
-G $NEWGTF \
--transcriptome-index=./tophat \
../bowtie2/$genomegrab \
read1.fq

echo "did we do it? exit says $?"

rm read1.fq

exit 0
