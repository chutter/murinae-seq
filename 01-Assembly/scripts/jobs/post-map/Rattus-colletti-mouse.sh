#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/seq-run-2/52951_BH7T77BCX2_GCTAATC_S14_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/seq-run-2/52951_BH7T77BCX2_GCTAATC_S14_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/seq-run-1/52951_BH7T33BCX2_GCTAATC_S14_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-colletti/Rattus-colletti-flagstats.txt