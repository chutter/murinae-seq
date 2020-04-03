#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/hiseq4000-paired-1/SM056_CKDL190143343-1a-D711-AK1680_H7275BBXX_L2_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Niviventer-rapit/Niviventer-rapit-flagstats.txt
