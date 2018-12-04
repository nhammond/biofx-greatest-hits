#!/bin/bash

if [[ -n $1 ]]
then
    PATH=$1
else
    PATH=$(pwd)/truthsets/NA24385/
fi

mkdir -p $PATH

SOURCE_PATH=ftp://anonymous:anonymous@ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/
SOURCE_FILES="md5sum HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi"

for file in ${SOURCE_FILES}
do
  wget -P ${PATH} ${SOURCE_PATH}/${file}
done
