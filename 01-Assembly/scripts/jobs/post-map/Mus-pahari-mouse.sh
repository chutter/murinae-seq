#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/hiseq2000-paired-2/lane1_TCGCAG_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/hiseq2000-paired-1/10460X1_131030_SN141_0734_AC2RLAACXX_7_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Mus-pahari/Mus-pahari-flagstats.txt
