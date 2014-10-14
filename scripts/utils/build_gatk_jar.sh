#!/usr/bin/env bash
# build a GATK jar without embedded dependencies from current git
mkdir -p ~/tmp/gatkbuild
cd ~/tmp/gatkbuild
git clone https://github.com/chapmanb/gatk.git
cd gatk
ant
cd dist
mkdir -p nodeps
for x in GenomeAnalysisTK.jar Aligner.jar Queue.jar StingUtils.jar picard-private-parts*jar
do
    cp $x nodeps
done
cd nodeps
for x in *.jar
do
    jar xvf $x
    rm -f $x
done
jar cvmf META-INF/MANIFEST.MF ../gatk.jar *
