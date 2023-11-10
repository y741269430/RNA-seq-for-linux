# RNA-seq （Linux上游分析，获取rawcounts）

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
先到网上下载：https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/  
    
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
    # -x ${mm39} \
    # -U ${i}_1.clean.fq.gz \
    # -S ./bam/${i}.sam 2> ./bam/${i}_map.txt &
    
    done

## 4. sam to bam    

    vim rna2_sam2bam.sh

    #!/bin/bash
    ## sam to bam (samtools) ##
    ## sorted by read name (samtools) for rowcounts ##
    ## sorted by position (samtools) for rmats ##

    cat filenames | while read i; 
    do
    nohup samtools view -@ 4 -S ./bam/${i}.sam -b | samtools sort -@ 4 -n -o ./bam/${i}-sorted-name.bam &  

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

## 6. Remove redundant rows     

    vim rna4_rmcounts.sh

    #!/bin/bash
    ## remove redundant rows ##

    tail -n 5 ./rawcounts/* > ./rawcounts/total.info

    cat filenames | while read i; 
    do
    sed -i '/process/d;/__/d;/retrieve/d' ./rawcounts/${i}.count & 
    done

## hisat2输出解读  

    21800552 reads; of these:  #一共读取的reads数目，共分三部分
      21800552 (100.00%) were paired; of these:  #第一部分是paired-end模式下比对结果一致的结果
        692176 (3.18%) aligned concordantly 0 times  #在比对一致的结果中，不太合理的比对。即read1和read2都能比对上但是不合理（方向不对或者片段长度不对）
        20143802 (92.40%) aligned concordantly exactly 1 time   #在这些比对一致结果中，合理的比对中有92.4是reads只匹配到一处
        964574 (4.42%) aligned concordantly >1 times  #在这些比对的结果中，合理的比对中，有4.42reads匹配到了多处
        ----
        692176 pairs aligned concordantly 0 times; of these:  #paired-end 模式下比对不合理的比对
          245626 (35.49%) aligned discordantly 1 time
        ----
        446550 pairs aligned 0 times concordantly or discordantly; of these: # 剩余的reads
          893100 mates make up the pairs; of these:
            620889 (69.52%) aligned 0 times
            223389 (25.01%) aligned exactly 1 time
            48822 (5.47%) aligned >1 times
    98.58% overall alignment rate  #有98.58的匹配率

## rawcounts输出解读  

https://htseq.readthedocs.io/en/master/count.html#usage  

    __no_feature: reads which could not be assigned to any feature.
    __ambiguous: reads which could have been assigned to more than one feature and hence were not counted for any of these.
    __too_low_aQual: reads which were skipped due to the -a option.
    __not_aligned: reads in the SAM file without alignment.
    __alignment_not_unique: reads with more than one reported alignment. These reads are recognized from the NH optional SAM field tag.

    __no_feature              #不能对应到任何单位类型的reads数
    __ambiguous               #不能判断落在那个单位类型的reads数
    __too_low_aQual           #低于-a设定的reads mapping质量的reads数
    __not_aligned             #存在于SAM文件，但没有比对上的reads数
    __alignment_not_unique    #比对到多个位置的reads数
    
## seqtk提取随机10000条reads  
    
    nohup seqtk sample -s100 BL6_6001_1.clean.fq.gz 10000 | gzip > exp_6001_1.clean.fq.gz &
    nohup seqtk sample -s100 BL6_6001_2.clean.fq.gz 10000 | gzip > exp_6001_2.clean.fq.gz &
    nohup seqtk sample -s100 BL6_6002_1.clean.fq.gz 10000 | gzip > exp_6002_1.clean.fq.gz &
    nohup seqtk sample -s100 BL6_6002_2.clean.fq.gz 10000 | gzip > exp_6002_2.clean.fq.gz &

## rnabash.sh 报错 ##
把所有代码合并，一起跑

    vim rnabash.sh
    #!/bin/bash
    threads=2
    
    mm39="~/downloads/genome/mm39_GRCm39/hisat2_idx/GRCm39"
    gtf="/home/jjyang/downloads/genome/mm39_GRCm39/gencode.vM27.annotation.gtf"
    
    cat filenames | while read i; 
    do
        nohup hisat2 -p ${threads} \
        -x ${mm39} \
        -1 ./nuohe_raw/${i}_1.fq.gz \
        -2 ./nuohe_raw/${i}_2.fq.gz \
        -S ./bam/${i}.sam 2> ./bam/${i}_map.txt & 
    done
    
    wait  
    
    cat filenames | while read i; 
    do
        nohup samtools view -@ ${threads} -S ./bam/${i}.sam -b | samtools sort -@ 8 -n -o ./bam/${i}-sorted-name.bam &  
    done
    
    wait  
    
    cat filenames | while read i; 
    do
        nohup htseq-count -n ${threads} \
        -f bam \
        -r name \
        -s no ./bam/${i}-sorted-name.bam ${gtf} > ./rawcounts/${i}.count &
    done
    
    wait  
    
    tail -n 5 ./rawcounts/* > ./rawcounts/total.info
    
    cat filenames | while read i; 
    do
        sed -i '/process/d;/__/d;/retrieve/d' ./rawcounts/${i}.count & 
    done
