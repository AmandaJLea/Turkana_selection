#!/bin/bash
#SBATCH --get-user-env

# NOTE: SAMPLEINFO is replaced with the appropriate sample name

f=SAMPLEINFO
R1=/Genomics/ayroleslab2/alea/archive_raw_fastq/WilsonLab_Turkana_WGS/SAMPLEINFO_R1_001.fastq.gz
R2=/Genomics/ayroleslab2/alea/archive_raw_fastq/WilsonLab_Turkana_WGS/SAMPLEINFO_R2_001.fastq.gz
R1_trim=/scratch/tmp/ayroles/alea_files/turkana_ASU/$f.R1.trim.fastq.gz
R2_trim=/scratch/tmp/ayroles/alea_files/turkana_ASU/$f.R2.trim.fastq.gz

cutadapt --nextseq-trim 20 -e 0.05 --overlap 2 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -a "G{101}" --minimum-length=20 --trim-n -o $R1_trim -p $R2_trim $R1 $R2

echo 'trimming done'

bwa mem -t 10 /Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_all_chr.fa $R1_trim $R2_trim > /scratch/tmp/ayroles/alea_files/turkana_ASU/$f.hg38.sam

echo 'mapping done'

module load samtools
samtools view -Sbq 1 /scratch/tmp/ayroles/alea_files/turkana_ASU/$f.hg38.sam > /scratch/tmp/ayroles/alea_files/turkana_ASU/$f.hg38.bam
echo $f

rm /scratch/tmp/ayroles/alea_files/turkana_ASU/$f.hg38.sam
rm $R1_trim
rm $R2_trim

