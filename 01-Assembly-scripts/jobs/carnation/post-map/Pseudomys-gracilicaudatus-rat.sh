#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/seq-run-2/52934_AH7T7NBCX2_TTCAACC_S17_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/seq-run-2/52934_AH7T7NBCX2_TTCAACC_S17_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/seq-run-1/52934_AH7T2JBCX2_TTCAACC_S17_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-gracilicaudatus/Pseudomys-gracilicaudatus-flagstats.txt