#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/seq-run-2/53927_HTVLWBCX2_CGACCTG_S15_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/seq-run-2/53927_HTVLWBCX2_CGACCTG_S15_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/seq-run-1/53927_HTVLKBCX2_CGACCTG_S15_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/seq-run-1/53927_HTVLKBCX2_CGACCTG_S15_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-novaehollandiae/Pseudomys-novaehollandiae-flagstats.txt
