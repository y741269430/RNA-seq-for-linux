#!/bin/bash
## sam to bam (samtools) ##
## sorted by read name (samtools) for rowcounts ##
## sorted by position (samtools) for rmats ##

cat filenames | while read i; 
do
nohup samtools view -@ 4 -S ./bam/${i}.sam -b | samtools sort -@ 4 -n -o ./bam/${i}-sorted-name.bam &  

done