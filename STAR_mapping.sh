#!/usr/bin/env bash

index_file=$1
R1=$2
R2=$3
label=$4
gtf=$5

COL3=$label

module load star/2.5.3a
module load samtools/1.7
module load bedtools/2.29.2

STAR --chimOutType WithinBAM --runThreadN 8 --readFilesCommand gunzip -c --genomeDir $index_file --readFilesIn $R1 --outFileNamePrefix ${COL3}_ --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.06 --limitOutSJcollapsed 2000000 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx

featureCounts -a $gtf -o gene_assigned -R BAM ${COL3}_Aligned.sortedByCoord.out.bam -T 4

mv ${COL3}_Aligned.sortedByCoord.out.bam.featureCounts.bam ${COL3}.bam

samtools sort -o ${COL3}.st.bam ${COL3}.bam

samtools index ${COL3}.st.bam

rm ${COL3}.bam

umi_tools dedup --stdin=${COL3}.st.bam --log=${COL3}.dedup.log --output-stats=${COL3}.stats.tsv --per-cell --per-gene --gene-tag=XT  > ${COL3}.dedup.bam

samtools index ${COL3}.dedup.bam

umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --wide-format-cell-counts --per-cell -I ${COL3}.dedup.bam -S ${COL3}.UMI_counts.tsv.gz

bedtools bamtobed -i ${COL3}.st.bam > ${COL3}.R1.bed





