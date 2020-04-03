#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-1/JAE4870_S66_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-1/JAE4870_S66_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-1/JAE4870_S66_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-1/JAE4870_S66_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-paired-2/JAE4870_S35_L004_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-paired-2/JAE4870_S35_L003_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-paired-2/JAE4870_S35_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-paired-2/JAE4870_S35_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/miseq-paired/JAE4870_S40_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-2/JAE4870_S66_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-2/JAE4870_S66_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-2/JAE4870_S66_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-single-2/JAE4870_S66_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/nextseq-paired-1/JAE4870_S40_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Paruromys-dominator/Paruromys-dominator-flagstats.txt
