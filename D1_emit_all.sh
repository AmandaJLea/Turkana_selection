#########
# RUN IN BASH
# NOTE, INTERVAL was replaced with a genomic interval such as chr10:1-5000000
#########

#!/bin/bash

# for f in `cat /Genomics/ayroleslab2/alea2/chr_pieces.txt` ; do cat high_cov2.sh | sed -e s/INTERVAL/$f/g > high_cov2.$f.sh; done
# rm commands1.sh; touch commands1.sh; for f in `cat /Genomics/ayroleslab2/alea2/chr_pieces.txt` ; do echo "sh high_cov2.$f.sh" >> commands1.sh; done

module load java
module load samtools

chrom=INTERVAL
in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_database=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/$chrom
in_samples=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/all_high_cov.sample_map
out_vcf_gz=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/$chrom.all_high_cov.emit_all.vcf.gz

echo $chrom > temp.${chrom}.list

$in_gatk/gatk --java-options "-Xmx80g" GenotypeGVCFs -R $in_genome -L temp.${chrom}.list -V gendb://$in_database -O $out_vcf_gz --tmp-dir=/scratch/tmp/ayroles/alea_files --max-alternate-alleles 2 --include-non-variant-sites 

#########
# RUN IN BASH
# NOTE, CHROMNUM was replaced with an autosomal chromsome number
#########

module load bcftools

f=CHROMNUM
bcftools concat -a -Oz -o /scratch/tmp/alea/chr${f}.all_high_cov.emit_all.vcf.gz /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/chr${f}:*.all_high_cov.emit_all.vcf.gz

module load vcftools

vcftools --gzvcf chrCHROMNAME.all_high_cov.emit_all.vcf.gz --counts2 --out chrCHROMNAME

sed '1d' chrCHROMNAME.frq.count | cut -f 6 | sort -g | uniq -c > chrCHROMNAME.sfs.txt

