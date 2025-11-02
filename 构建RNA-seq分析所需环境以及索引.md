# 构建RNA-seq分析所需环境以及索引      
## 1. 使用本人部署的环境进行安装    
- [environment.yml文件下载](https://github.com/y741269430/RNA-seq-for-linux/blob/main/environment.yml)
```bash
# conda env export > environment.yml
conda env create -f environment.yml
```
该环境是由以下命令构建的    
- [RNA-seq(5):序列比对：Hisat2](https://www.jianshu.com/p/479c7b576e6f)    
```bash
conda create -n rnaseq python=3.7
conda activate rnaseq
conda install -c conda-forge aria2
pip install HTSeq==2.0.0 -i https://pypi.tuna.tsinghua.edu.cn/simple
conda install -c bioconda hisat2=2.2.1
conda install -c bioconda samtools
conda install -c bioconda seqtk
```
## 2.构建小鼠基因组mm39的index 
其实hisat2-buld在运行的时候也会自己寻找exons和splice_sites，但是先做的目的是为了提高运行效率  
先到网上下载小鼠mm39的基因组：ftp://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/
```bash
cd downloads/genome/mm_v112
mkdir hisat2_idx

# 这一步是提取gtf文件中的外显子和可变剪切的位点(这两个py脚本藏在了conda的目录里面)    
nohup ~/.conda/envs/rnaseq/bin/hisat2_extract_exons.py Mus_musculus.GRCm39.112.chr.gtf > v112.exons.gtf &
nohup ~/.conda/envs/rnaseq/bin/hisat2_extract_splice_sites.py Mus_musculus.GRCm39.112.chr.gtf > v112.splice_sites.gtf &

# 构建index， 必选项是基因组所在文件路径和输出的前缀  
nohup hisat2-build -p 60 --ss v112.splice_sites.gtf --exon v112.exons.gtf Mus_musculus.GRCm39.dna_sm.chromosome.chr.fa ./hisat2_idx/GRCm39 &
```
