#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/hiseq4000-paired-1/SM059_CKDL190143344-1a-D705-AK1682_H7275BBXX_L3_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Genus-sp-nov/Genus-sp-nov-flagstats.txt