#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/seq-run-2/52926_AH7T7NBCX2_TCAGCTT_S9_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/seq-run-2/52926_AH7T7NBCX2_TCAGCTT_S9_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/seq-run-1/52926_AH7T2JBCX2_TCAGCTT_S9_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-chapmani/Pseudomys-chapmani-flagstats.txt
