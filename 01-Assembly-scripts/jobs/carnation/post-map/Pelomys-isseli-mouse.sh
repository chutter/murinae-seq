#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/hiseq4000-paired-1/SM068_CKDL190143343-1a-D708-AK1681_H7275BBXX_L2_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pelomys-isseli/Pelomys-isseli-flagstats.txt