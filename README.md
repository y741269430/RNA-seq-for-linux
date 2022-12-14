# RNA-seq

## 0. Build source used for RNA-seq  

    conda create -n rnaseq python=3.7
    conda activate rnaseq
    conda install -c conda-forge aria2
    pip install HTSeq==2.0.0 -i https://pypi.tuna.tsinghua.edu.cn/simple
    conda install -c bioconda hisat2=2.2.1
    conda install -c bioconda samtools
    conda install -c bioconda seqtk
    
## 0. Build the hisat2 reference genome index (mm39)  

其实hisat2-buld在运行的时候也会自己寻找exons和splice_sites，但是先做的目的是为了提高运行效率  
    
    mkdir hisat2_idx
    
    cd /home/yangjiajun/downloads/genome/mm39_GRCm39/
    
    nohup /home/yangjiajun/miniconda3/envs/rnaseq/bin/hisat2_extract_exons.py gencode.vM27.annotation.gtf > vM27.exons.gtf &
    nohup /home/yangjiajun/miniconda3/envs/rnaseq/bin/hisat2_extract_splice_sites.py gencode.vM27.annotation.gtf > vM27.splice_sites.gtf &

建立index， 必选项是基因组所在文件路径和输出的前缀  

    nohup hisat2-build -p 60 --ss vM27.splice_sites.gtf --exon vM27.exons.gtf ./ucsc_fa/GRCm39.genome.fa ./hisat2_idx/GRCm39 &

## 1. Activate the source and create the folder  
    
    conda activate rnaseq  
    
    mkdir bam rawcounts 
    
## 2. Write the filenames  

    ls *1.clean* |cut -d "_" -f 1 > filenames

## 3. Alignment to mm39  

    vim rna1_hisat2.sh

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

    # single end
    # nohup hisat2 -p 4 \
    -x ${mm39} \
    -U ${i}_1.clean.fq.gz \
    -S ./bam/${i}.sam 2> ./bam/${i}_map.txt &
    
    done

## 4. sam to bam    

    vim rna2_sam2bam.sh

    #!/bin/bash
    ## sam to bam (samtools) ##
    ## sorted by read name (samtools) for rowcounts ##
    ## sorted by position (samtools) for rmats ##

    cat filenames | while read i; 
    do
    nohup samtools view -@ 4 -S ./bam/${i}.sam -b > ./bam/${i}.bam &&
    samtools sort -@ 4 -n ./bam/${i}.bam -o ./bam/${i}-sorted-name.bam &

    done

## 5. htseq-count    

    vim rna3_htcounts.sh

    #!/bin/bash
    ## calculate rawcounts (htseq-count) ##

    gtf="/home/yangjiajun/downloads/genome/mm39_GRCm39/gencode.vM27.annotation.gtf"

    cat filenames | while read i; 
    do

    nohup htseq-count -n 10 \
    -f bam \
    -r name \
    -s no ./bam/${i}-sorted-name.bam ${gtf} > ./rawcounts/${i}.count &

    done

## 6. Remove row      

    vim rna4_rmcounts.sh

    #!/bin/bash
    ## remove redundant rows ##

    tail -n 5 ./rawcounts/* > ./rawcounts/total.info

    cat filenames | while read i; 
    do
    sed -i '/process/d;/__/d;/retrieve/d' ./rawcounts/${i}.count & 
    done




