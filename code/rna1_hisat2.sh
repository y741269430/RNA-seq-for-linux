#!/bin/bash
## mapping (hisat2) ##

mm39="/home/jjyang/downloads/genome/mm39_GRCm39/hisat2_idx/GRCm39"  # 基因组所在位置

cat filenames | while read i; 
do
nohup hisat2 -p 4 \
-x ${mm39} \
-1 XK-250716T/${i}_1.fq.gz \
-2 XK-250716T/${i}_2.fq.gz \
-S ./bam/${i}.sam 2> ./mapinfo/${i}_map.txt & 

# 以下是单端比对的代码
# single end
# nohup hisat2 -p 4 \
# -x ${mm39} \
# -U ${i}_1.fq.gz \
# -S ./bam/${i}.sam 2> ./mapinfo/${i}_map.txt &
    
done
