#!/bin/bash

if [[ -n $1 ]]
then
    PATH=$1
else
    PATH=$(pwd)/genomes/NA24385/
fi

mkdir -p $PATH

SOURCE_PATH=ftp://anonymous:anonymous@ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/
SOURCE_FILES="README HG002.hs37d5.2x250.bam HG002.hs37d5.2x250.bam.bai"

for file in ${SOURCE_FILES}
do
  wget -P ${PATH} ${SOURCE_PATH}/${file}
done
