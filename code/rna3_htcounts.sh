#!/bin/bash
## calculate rawcounts (htseq-count) ##

gtf="/home/jjyang/downloads/genome/mm39_GRCm39/gencode.vM27.annotation.gtf"   # gtf所在位置

cat filenames | while read i; 
do
nohup htseq-count -n 10 \
-f bam \
-r name \
-s no ./bam/${i}-sorted-name.bam ${gtf} > ./rawcounts/${i}.count &
done
