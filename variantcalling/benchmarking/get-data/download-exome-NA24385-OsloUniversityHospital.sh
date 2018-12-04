#!/bin/bash

if [[ -n $1 ]]
then
    PATH=$1
else
    PATH=$(pwd)/genomes/NA24385/
fi

mkdir -p $PATH

SOURCE_PATH=ftp://anonymous:anonymous@ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/
SOURCE_FILES="151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai"

for file in ${SOURCE_FILES}
do
  wget -P ${PATH} ${SOURCE_PATH}/${file}
done
