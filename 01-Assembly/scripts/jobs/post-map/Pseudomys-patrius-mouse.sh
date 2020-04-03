#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/seq-run-2/52940_BH7T77BCX2_CCAACCT_S3_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/seq-run-2/52940_BH7T77BCX2_CCAACCT_S3_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/seq-run-1/52940_BH7T33BCX2_CCAACCT_S3_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-patrius/Pseudomys-patrius-flagstats.txt
