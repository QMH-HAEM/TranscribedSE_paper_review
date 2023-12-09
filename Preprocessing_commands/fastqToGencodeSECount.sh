#!/bin/bash
set -e
set -u
set -o pipefail

CPU=16
GTF=~/Ref/gencode.v30.primary_assembly.annotation.gtf

sample_list=$1
sample_names=($(cut -f 2 "$sample_list"))

for SN in ${sample_names[@]}
do
   echo $(date +%Y/%m/%d\ %H:%M): "Processing ${SN}..."
   start=$(date +%Y%m%d_%H%M)
   mkdir -p ./${SN}_$start
   DATAOUT=./${SN}_$start

   trim_galore --phred33 --quality 20 --cores $CPU --stringency 5 --length 20 --paired ${SN}_1.fastq.gz ${SN}_2.fastq.gz
   mv ${SN}_1_val_1.fq.gz ./$DATAOUT/
   mv ${SN}_2_val_2.fq.gz ./$DATAOUT/
   mv ${SN}_1.fastq.gz_trimming_report.txt ./$DATAOUT/
   mv ${SN}_2.fastq.gz_trimming_report.txt ./$DATAOUT/

   STAR --runThreadN $CPU --genomeDir ~/Ref/GRCh38_gencode30_STAR_50/ --readFilesIn ./$DATAOUT/${SN}_1_val_1.fq.gz ./$DATAOUT/${SN}_2_val_2.fq.gz --readFilesCommand zcat --sjdbGTFfile $GTF BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 5 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic --runMode alignReads --outFileNamePrefix ./$DATAOUT/${SN}. --outBAMcompression 1 --outBAMsortingThreadN 10 --outWigType bedGraph --outWigNorm RPM

   featureCounts -p -B -T $CPU -F SAF -a ~/Projects/38.superenhancer/sEnh_gencode30.annotation.saf -o ./$DATAOUT/${SN}.count -s 0 $DATAOUT/"${SN}.Aligned.sortedByCoord.out.bam"
   featureCounts -p -B -T $CPU -a $GTF -t exon -g gene_id -o ./$DATAOUT/${SN}_gencode30.count -s 2 $DATAOUT/"${SN}.Aligned.sortedByCoord.out.bam"

   echo $(date +%Y/%m/%d\ %H:%M): "Finished processing ${SN}."
done
