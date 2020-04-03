#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/seq-run-2/53922_HTVLWBCX2_CCGGTAC_S10_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/seq-run-2/53922_HTVLWBCX2_CCGGTAC_S10_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/seq-run-1/53922_HTVLKBCX2_CCGGTAC_S10_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/seq-run-1/53922_HTVLKBCX2_CCGGTAC_S10_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-bolami/Pseudomys-bolami-flagstats.txt