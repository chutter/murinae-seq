#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/seq-run-2/53930_HTVLWBCX2_TGCGTCC_S18_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/seq-run-2/53930_HTVLWBCX2_TGCGTCC_S18_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/seq-run-1/53930_HTVLKBCX2_TGCGTCC_S18_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/seq-run-1/53930_HTVLKBCX2_TGCGTCC_S18_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Conilurus-albipes/Conilurus-albipes-flagstats.txt
