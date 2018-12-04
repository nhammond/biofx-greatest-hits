#!/bin/bash

if [[ -n $1 ]]
then
    PATH=$1
else
    PATH=$(pwd)/reference_genomes/human_g1k_v37_decoy/
fi

wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.g
z' -O ${PATH}/chromFa.tar.gz
sudo tar -xf ${PATH}/chromFa.tar.gz

(cd ${PATH} && cat chrM.fa chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa chr1_gl000191_random.fa chr1_gl000192_random.fa chr4_ctg9_hap1.fa chr4_gl000193_random.fa chr4_gl000194_random.fa chr6_apd_hap1.fa chr6_cox_hap2.fa chr6_dbb_hap3.fa chr6_mann_hap4.fa chr6_mcf_hap5.fa chr6_qbl_hap6.fa chr6_ssto_hap7.fa chr7_gl000195_random.fa chr8_gl000196_random.fa chr8_gl000197_random.fa chr9_gl000198_random.fa chr9_gl000199_random.fa chr9_gl000200_random.fa chr9_gl000201_random.fa chr11_gl000202_random.fa chr17_ctg5_hap1.fa chr17_gl000203_random.fa chr17_gl000204_random.fa chr17_gl000205_random.fa chr17_gl000206_random.fa chr18_gl000207_random.fa chr19_gl000208_random.fa chr19_gl000209_random.fa chr21_gl000210_random.fa chrUn_gl000211.fa chrUn_gl000212.fa chrUn_gl000213.fa chrUn_gl000214.fa chrUn_gl000215.fa chrUn_gl000216.fa chrUn_gl000217.fa chrUn_gl000218.fa chrUn_gl000219.fa chrUn_gl000220.fa chrUn_gl000221.fa chrUn_gl000222.fa chrUn_gl000223.fa chrUn_gl000224.fa chrUn_gl000225.fa chrUn_gl000226.fa chrUn_gl000227.fa chrUn_gl000228.fa chrUn_gl000229.fa chrUn_gl000230.fa chrUn_gl000231.fa chrUn_gl000232.fa chrUn_gl000233.fa chrUn_gl000234.fa chrUn_gl000235.fa chrUn_gl000236.fa chrUn_gl000237.fa chrUn_gl000238.fa chrUn_gl000239.fa chrUn_gl000240.fa chrUn_gl000241.fa chrUn_gl000242.fa chrUn_gl000243.fa chrUn_gl000244.fa chrUn_gl000245.fa chrUn_gl000246.fa chrUn_gl000247.fa chrUn_gl000248.fa chrUn_gl000249.fa > hg19.fa) 

samtools faidx ${PATH}/hg19.fa
