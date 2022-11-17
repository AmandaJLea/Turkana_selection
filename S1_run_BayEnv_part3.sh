#!/bin/bash

# NOTE: CHROMNUM was replaced with the chromosome number

# WGS order = El Molo Ik LWK Pokot_Ngitepes_Karamojong Rendille Samburu Turkana YRI

SNPFILE=/Genomics/ayroleslab2/alea2/28May20_BayEnv_chrCHROMNUM_wLWKYRI.txt
ENVFILE=/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/bayenv_files/wgs_environment.txt
MATFILE=/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/bayenv_files/array_cov_matrix_8pops_final.txt
POPNUM=8
ITNUM=100000
ENVNUM=1
rnd1=$(perl -e 'printf("%05d",rand(99999))')
rnd2=$(perl -e 'printf("%05d",rand(99999))')
rnd3=$(perl -e 'printf("%05d",rand(99999))')

rm 28May20_BayEnv_chrCHROMNUM_output*
f=$(< "$SNPFILE" wc -l)
seq 2 2 $f > temp_chrCHROMNUM_loci.txt

for f in `cat temp_chrCHROMNUM_loci.txt`; 
do head -$f $SNPFILE | tail -2 > 28May20_temp_chrCHROMNUM_$f.txt; 
~/programs/tguenther-bayenv2_public-2b2b7f20bb62/bayenv2 -i 28May20_temp_chrCHROMNUM_$f.txt -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -X -r $rnd1 -o 28May20_BayEnv_chrCHROMNUM_output1; 
~/programs/tguenther-bayenv2_public-2b2b7f20bb62/bayenv2 -i 28May20_temp_chrCHROMNUM_$f.txt -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -X -r $rnd2 -o 28May20_BayEnv_chrCHROMNUM_output2; 
~/programs/tguenther-bayenv2_public-2b2b7f20bb62/bayenv2 -i 28May20_temp_chrCHROMNUM_$f.txt -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -X -r $rnd3 -o 28May20_BayEnv_chrCHROMNUM_output3; 
rm 28May20_temp_chrCHROMNUM_$f.txt; rm 28May20_temp_chrCHROMNUM_$f.txt.freqs; done

