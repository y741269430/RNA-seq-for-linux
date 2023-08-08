#!/bin/bash
## mapping (hisat2) ##

mm39="/home/yangjiajun/downloads/genome/mm39_GRCm39/hisat2_idx/GRCm39"

cat filenames | while read i; 
do
nohup hisat2 -p 4 \
-x ${mm39} \
-1 ${i}_1.clean.fq.gz \
-2 ${i}_2.clean.fq.gz \
-S ./bam/${i}.sam 2> ./bam/${i}_map.txt & 

done