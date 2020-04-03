#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-1/SMG3823_S60_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-1/SMG3823_S60_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-1/SMG3823_S60_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-1/SMG3823_S60_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-paired-2/SMG3823_S28_L002_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-paired-2/SMG3823_S28_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-paired-2/SMG3823_S28_L004_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-paired-2/SMG3823_S28_L003_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/miseq-paired/SMG3823_S34_L001_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-2/SMG3823_S60_L003_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-2/SMG3823_S60_L002_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-2/SMG3823_S60_L004_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-single-2/SMG3823_S60_L001_R1_001.fastp.decon.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/nextseq-paired-1/SMG3823_S34_001.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Apodemus-sylvaticus/Apodemus-sylvaticus-flagstats.txt
