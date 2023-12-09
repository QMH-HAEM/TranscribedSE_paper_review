#!/bin/bash
################################################################################
# This pipeline processes ChIP-seq data available in GEO and outputs analysis-
# ready de-duplicated BAM file for MACS2.
#
# Need to modify the 2colour/quality and length options for Trim-Galore.
#
# Created on: 06/07/2021
# Modified on: 03/01/2023
################################################################################
set -e
set -u
set -o pipefail

echo "$(date +%Y-%m-%d\ %H:%M:%S): Performing Trim-Galore for $1..."
trim_galore --gzip --2colour 20 --cores 4 --length 25 --paired ./${1}_1.fastq ./${1}_2.fastq --output_dir ./
#trim_galore --gzip --quality 20 --cores 6 --length 25 ./${1}.fastq.gz --output_dir ./
echo "$(date +%Y-%m-%d\ %H:%M:%S): Performed Trim-Galore for $1."

echo "$(date +%Y-%m-%d\ %H:%M:%S): Performing BWA for $1..."
bwa mem -t 4 -M ~/Ref/GRCh38.primary_assembly.genome.fa ${1}_1_val_1.fq.gz ${1}_2_val_2.fq.gz > ${1}.sam
echo "$(date +%Y-%m-%d\ %H:%M:%S): Performed BWA for $1."

echo "$(date +%Y-%m-%d\ %H:%M:%S): Performing SortSam for $1..."
gatk --java-options "-Xmx8g -Xms8g" SortSam -I ${1}.sam -O ${1}_sort.bam --SORT_ORDER coordinate
echo "$(date +%Y-%m-%d\ %H:%M:%S): Performed SortSam for $1."

echo "$(date +%Y-%m-%d\ %H:%M:%S): Performing AddOrReplaceReadGroups for $1..."
gatk --java-options "-Xmx8g -Xms8g" AddOrReplaceReadGroups -I ${1}_sort.bam -O ${1}_sort_RG.bam --RGID ${1} --RGLB CHIPSEQ --RGPL ILLUMINA --RGPU unit1 --RGSM ${1}
echo "$(date +%Y-%m-%d\ %H:%M:%S): Performed AddOrReplaceReadGroups for $1."

echo "$(date +%Y-%m-%d\ %H:%M:%S): Performing MarkDuplicates for $1..."
gatk --java-options "-Xmx8g -Xms8g" MarkDuplicates -I ${1}_sort_RG.bam -O ${1}_sort_RG_MD.bam -M ${1}_MD_metrics.txt --CREATE_INDEX true --REMOVE_DUPLICATES true
echo "$(date +%Y-%m-%d\ %H:%M:%S): Performed MarkDuplicates for $1."

#rm -f ${1}.fastq.gz
rm -f ${1}.*_trimming_report.txt
rm -f ${1}_sort.bam
rm -f ${1}_sort_RG.bam
rm -f ${1}_*.fq.gz
rm -f ${1}.sam
rm -f ${1}_*.fastq
