#!/usr/bin/env bash

index_file=$1
R1=$2
R2=$3
label=$4


COL3=$label

module load bwa/0.7.16a
module load samtools/1.7
module load bedtools/2.29.2

bwa mem -t 8 $index_file $R1 $R2 | samtools view -bS - > ${COL3}.bam

samtools sort -o ${COL3}.st.bam ${COL3}.bam

samtools index ${COL3}.st.bam

rm ${COL3}.bam

samtools view -b -F 2048 -F 256 -f 2 ${COL3}.st.bam > ${COL3}.filter.bam

samtools index ${COL3}.filter.bam

bedtools bamtobed -i ${COL3}.filter.bam | grep "/1" > ${COL3}.R1.bed



