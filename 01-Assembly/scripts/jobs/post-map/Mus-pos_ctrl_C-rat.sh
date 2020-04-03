#!/bin/bash
samtools merge -r /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C.sorted.bam /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/hiseq4000-paired-1/M44_C11_CKDL190143344-1a-D711-AK1682_H7275BBXX_L3_.fastp.decon.bam
java -jar ~/bin/picard.jar MarkDuplicates I=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C.sorted.bam O=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C.sorted.mkdup.bam M=/scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-mkdup-metrics.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/targets-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-target-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-target-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-avg-target-depth.txt
samtools depth -b /scratch/gregg_thomas/Murinae-seq/Targets/tiles-rnor6-coords.bed /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-tile-depth.tab
awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-tile-depth.tab > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-avg-tile-depth.txt
samtools flagstat /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C.sorted.mkdup.bam > /scratch/gregg_thomas/Murinae-seq/03B-MappedRat/Mus-pos_ctrl_C/Mus-pos_ctrl_C-flagstats.txt
