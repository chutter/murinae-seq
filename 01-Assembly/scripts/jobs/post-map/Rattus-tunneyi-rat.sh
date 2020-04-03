#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/seq-run-2/52952_BH7T77BCX2_GACTTCT_S15_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/seq-run-2/52952_BH7T77BCX2_GACTTCT_S15_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/seq-run-1/52952_BH7T33BCX2_GACTTCT_S15_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-tunneyi/Rattus-tunneyi-flagstats.txt
