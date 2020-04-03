#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/seq-run-2/53920_HTVLWBCX2_GCTCGAA_S8_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/seq-run-2/53920_HTVLWBCX2_GCTCGAA_S8_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/seq-run-1/53920_HTVLKBCX2_GCTCGAA_S8_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/seq-run-1/53920_HTVLKBCX2_GCTCGAA_S8_L002_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-mm10-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedMouse/Conilurus-penicillatus/Conilurus-penicillatus-flagstats.txt
