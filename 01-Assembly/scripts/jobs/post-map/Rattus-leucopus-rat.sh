#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/hiseq-paired/leucopusWGA_ATGCCG_L008_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/hiseq-paired/leucopus-WGA_ATGCCG_L007_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/hiseq2500-paired/lane2_ATGCCG_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Rattus-leucopus/Rattus-leucopus-flagstats.txt
