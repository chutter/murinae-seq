#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-1/Z5248_S59_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-1/Z5248_S59_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-1/Z5248_S59_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-1/Z5248_S59_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-paired-2/Z5248_S27_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-paired-2/Z5248_S27_L003_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-paired-2/Z5248_S27_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-paired-2/Z5248_S27_L004_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/miseq-paired/Z5248_S33_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-2/Z5248_S59_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-2/Z5248_S59_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-2/Z5248_S59_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-single-2/Z5248_S59_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/nextseq-paired-1/Z5248_S33_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Hydromys-chrysogaster/Hydromys-chrysogaster-flagstats.txt
