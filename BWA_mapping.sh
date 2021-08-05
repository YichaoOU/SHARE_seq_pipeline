#!/usr/bin/env bash

index_file=$1
R1=$2
R2=$3
label=$4
rseqc_bed=$5
genome=$6

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


# Read distribution

module load conda3/202011

source activate /home/yli11/.conda/envs/rseqc

read_distribution.py -i ${COL3}.filter.bam -r $rseqc_bed > ${COL3}.read_distribution.tsv

source activate /home/yli11/.conda/envs/r_env

run_ATACseqQC.R ${COL3}.filter.bam $genome

