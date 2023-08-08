#!/bin/bash
## calculate rawcounts (htseq-count) ##

gtf="/home/yangjiajun/downloads/genome/mm39_GRCm39/gencode.vM27.annotation.gtf"

cat filenames | while read i; 
do
nohup htseq-count -n 4 \
-f bam \
-r name \
-s no ./bam/${i}-sorted-name.bam ${gtf} > ./rawcounts/${i}.count &
done