#!/bin/bash

# NOTE: CHROMNUM was replaced with the chromosome number and STARTNUM/ENDNUM were replaced with the start and end locations of 5Mb intervals

chr=CHROMNUM
start=STARTNUM	
end=ENDNUM

gen=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/chr${chr}.hg19.low_cov.GQ10.gen.gz

haps=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz
legend=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/1000GP_Phase3_chr${chr}.v2.legend
map=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt

phased_haps2=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased_v2.haps
phased_leg=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased_v2.leg

# impute using Turkana + 1000 Genomes reference panels
/Genomics/grid/users/alea/programs/impute_v2.3.2_x86_64_static/impute2 \
 -m $map \
 -h $phased_haps2 $haps \
 -l $phased_leg $legend  \
 -g $gen \
 -int $start $end \
 -Ne 20000 \
 -o /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute2panel.${chr}.${start}.${end}.hg19_GQ10

# impute using Turkana reference panel only
 /Genomics/grid/users/alea/programs/impute_v2.3.2_x86_64_static/impute2 \
 -m $map \
 -h $phased_haps2 \
 -l $phased_leg  \
 -g $gen \
 -int $start $end \
 -Ne 20000 \
 -o /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.${start}.${end}.hg19_GQ10