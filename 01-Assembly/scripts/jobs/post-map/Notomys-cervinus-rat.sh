#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/seq-run-2/52922_AH7T7NBCX2_AAGCTAA_S5_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/seq-run-2/52922_AH7T7NBCX2_AAGCTAA_S5_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/seq-run-1/52922_AH7T2JBCX2_AAGCTAA_S5_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Notomys-cervinus/Notomys-cervinus-flagstats.txt
