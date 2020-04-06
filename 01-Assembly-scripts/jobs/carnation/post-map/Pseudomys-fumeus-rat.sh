#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/seq-run-2/52929_AH7T7NBCX2_TAATCAT_S12_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/seq-run-2/52929_AH7T7NBCX2_TAATCAT_S12_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/seq-run-1/52929_AH7T2JBCX2_TAATCAT_S12_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-fumeus/Pseudomys-fumeus-flagstats.txt