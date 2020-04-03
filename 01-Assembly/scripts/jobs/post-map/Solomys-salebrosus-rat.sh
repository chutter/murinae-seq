#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/hiseq4000-paired-2/SM076_CKDL190143345-1a-D705-AK1780_H7275BBXX_L4_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Solomys-salebrosus/Solomys-salebrosus-flagstats.txt
