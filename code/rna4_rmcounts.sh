#!/bin/bash
## remove redundant rows ##

tail -n 5 ./rawcounts/* > ./rawcounts/total.info

cat filenames | while read i; 
do
sed -i '/process/d;/__/d;/retrieve/d' ./rawcounts/${i}.count & 
done
