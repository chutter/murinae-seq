#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/hiseq4000-paired-1/REC10_CKDL190143343-1a-D708-AK1544_H7275BBXX_L2_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Hylomyscus-stella/Hylomyscus-stella-flagstats.txt