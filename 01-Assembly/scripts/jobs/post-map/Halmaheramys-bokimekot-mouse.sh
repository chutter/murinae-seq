#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/hiseq4000-paired-2/REC22_CKDL190143345-1a-D710-AK1545_H7275BBXX_L4_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Halmaheramys-bokimekot/Halmaheramys-bokimekot-flagstats.txt
