#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/hiseq4000-paired-1/SM085_CKDL190143343-1a-DY0088-AK1543_H7275BBXX_L2_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Carpomys-melanurus/Carpomys-melanurus-flagstats.txt