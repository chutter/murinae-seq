#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/seq-run-2/52935_AH7T7NBCX2_TTAACTC_S18_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/seq-run-2/52935_AH7T7NBCX2_TTAACTC_S18_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/seq-run-1/52935_AH7T2JBCX2_TTAACTC_S18_L001_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Pseudomys-johnsoni/Pseudomys-johnsoni-flagstats.txt