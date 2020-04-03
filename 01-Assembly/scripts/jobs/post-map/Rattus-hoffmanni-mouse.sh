#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/hiseq-paired-no-WGA/hoffmani_GATCTC_L007_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/hiseq-paired-no-WGA/hoffmani_GATCTC_L008_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/hiseq4000-paired-2/SM074_CKDL190143345-1a-D709-AK1682_H7275BBXX_L4_.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/hiseq-paired/hoffmaniWGA_CCGATT_L008_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/hiseq-paired/hoffmani-WGA_CCGATT_L007_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/hiseq2500-paired-no-WGA/lane2_GATCTC_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/hiseq2500-paired/lane2_CCGATT_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Rattus-hoffmanni/Rattus-hoffmanni-flagstats.txt
