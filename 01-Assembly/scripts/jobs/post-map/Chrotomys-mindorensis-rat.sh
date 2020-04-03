#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-1/JAE520_S27_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-1/JAE520_S27_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-1/JAE520_S27_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-1/JAE520_S27_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-paired-2/JAE520_S34_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-paired-2/JAE520_S34_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-paired-2/JAE520_S34_L004_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-paired-2/JAE520_S34_L003_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/miseq-paired/JAE520_S4_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-2/JAE520_S27_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-2/JAE520_S27_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-2/JAE520_S27_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-single-2/JAE520_S27_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/nextseq-paired-1/JAE520_S4_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Chrotomys-mindorensis/Chrotomys-mindorensis-flagstats.txt
