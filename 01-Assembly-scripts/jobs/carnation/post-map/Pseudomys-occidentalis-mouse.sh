#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/seq-run-2/52939_BH7T77BCX2_AGCCTTG_S2_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/seq-run-2/52939_BH7T77BCX2_AGCCTTG_S2_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/seq-run-1/52939_BH7T33BCX2_AGCCTTG_S2_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Pseudomys-occidentalis/Pseudomys-occidentalis-flagstats.txt