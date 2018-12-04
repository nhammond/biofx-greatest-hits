#!/bin/bash

if [[ -n $1 ]]
then
    PATH=$1
else
    PATH=$(pwd)/reference_genomes/human_g1k_v37_decoy/
fi

wget --timestamping 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta' -O ${PATH}/human_g1k_v37_decoy.fasta
wget --timestamping 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta' -O ${PATH}/human_g1k_v37_decoy.dict

samtools faidx ${PATH}/human_g1k_v37_decoy.fasta
