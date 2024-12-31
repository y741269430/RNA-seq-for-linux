# 名称暂定

#### 1. 从人类基因组注释文件（GTF格式）中提取外显子和剪接位点信息，用于后续建立索引。    
这一步会生成两个文件：一个包含所有外显子的信息，另一个包含所有剪接位点的信息。    
    nohup /home/jjyang/.conda/envs/rnaseq/bin/hisat2_extract_exons.py gencode.v47.annotation.gtf > v47.exons.gtf &
    nohup /home/jjyang/.conda/envs/rnaseq/bin/hisat2_extract_splice_sites.py gencode.v47.annotation.gtf > v47.splice_sites.gtf &

#### 2. 使用上述提取的外显子和剪接位点信息，以及人类基因组序列，构建hisat2比对工具所需的索引。    
索引文件是一次性创建的，之后可以直接使用，无需重复构建。    
    nohup hisat2-build -p 40 --ss v47.splice_sites.gtf --exon v47.exons.gtf ./ucsc_fa/GRCh38.p14.genome.fa ./hisat2_idx/GRCh38 &

#### 3. 构建BLAST数据库，以便将RNA测序数据与自定义的序列进行比对。    
我们使用的是核苷酸（DNA/RNA）类型的数据库，并指定了输出路径。    
    makeblastdb -in S_genome.fa -dbtype nucl -parse_seqids -out mydna/mydna &

#### 4. 检查BLAST数据库是否成功构建。    
    blastdbcmd -db mydna/mydna -dbtype nucl -info

#### 5. 使用hisat2比对工具将RNA测序数据与人类基因组进行比对。    
这里我们调用了一个shell脚本`rna1_hisat2.sh`来执行比对操作。    
    bash rna1_hisat2.sh

#### 6. 将比对结果从SAM格式转换为FASTA格式，方便后续处理。    
SAM是比对结果的标准格式，而FASTA是序列的标准格式。    
    nohup samtools fasta -@ 20 -n MP.sam > MP.fa &
    nohup samtools fasta -@ 20 -n P.sam > P.fa &

#### 7. 使用BLAST工具将转换后的FASTA序列与自定义的序列库进行比对。    
输出格式设置为6，这是一种简洁的表格格式，便于解析。    
    blastn -db mydna/mydna -query bam/MP.fa -outfmt 6 > test_MP_1.txt &
    blastn -db mydna/mydna -query bam/P.fa -outfmt 6 > test_P_1.txt &

#### （ex）从自定义的序列中，提取比对后的序列    
    blastdbcmd -db mydna/mydna -entry_batch <(cut -f2 test_MP_1.txt) -out output_MP.fa &
    blastdbcmd -db mydna/mydna -entry_batch <(cut -f2 test_P_1.txt) -out output_P.fa &

#### 8. 从BLAST比对结果中提取唯一的查询ID，并根据这些ID从FASTA文件中筛选出对应的序列。    
    cut -f1 test_P_1.txt | sort -u > query_ids_P.txt
    seqkit grep -f query_ids_P.txt bam/P.fa > extracted_P.fa &
    
    cut -f1 test_MP_1.txt | sort -u > query_ids_MP.txt
    seqkit grep -f query_ids_MP.txt bam/MP.fa > extracted_MP.fa &

#### 9. 将筛选出来的序列再次与人类基因组进行比对，生成新的SAM文件。    
这一步是为了确认这些序列在人类基因组中的具体位置。    
    hisat2 -p 8 -x /home/jjyang/downloads/genome/h38_GRCh38/hisat2_idx/GRCh38 -f extracted_P.fa -S extracted_P.sam 2> ./extracted_P_map.txt &
    hisat2 -p 8 -x /home/jjyang/downloads/genome/h38_GRCh38/hisat2_idx/GRCh38 -f extracted_MP.fa -S extracted_MP.sam 2> ./extracted_MP_map.txt &

#### 10. 将新生成的SAM文件转换为BAM格式，并按照名称排序。    
BAM是压缩后的SAM格式，更加节省空间，也更适合进一步分析。    
    nohup samtools view -@ 8 -S extracted_MP.sam -b | samtools sort -@ 8 -n -o extracted_MP-sorted-name.bam &
    nohup samtools view -@ 8 -S extracted_P.sam -b | samtools sort -@ 8 -n -o extracted_P-sorted-name.bam &

#### 11. 使用HTSeq工具对排序后的BAM文件进行基因表达定量。    
结果会保存到指定的文件夹中，供后续分析使用。    
    nohup htseq-count -n 8 -f bam -r name -s no ./bam/extracted_MP-sorted-name.bam "/home/jjyang/downloads/genome/h38_GRCh38/gencode.v47.annotation.gtf" > ./rawcounts/extracted_MP.count &
    nohup htseq-count -n 8 -f bam -r name -s no ./bam/extracted_P-sorted-name.bam "/home/jjyang/downloads/genome/h38_GRCh38/gencode.v47.annotation.gtf" > ./rawcounts/extracted_P.count &

#### 12. 删除无用的行 得到count文件。    
    sed -i '/process/d;/__/d;/retrieve/d' ./rawcounts/extracted_MP.count &
    sed -i '/process/d;/__/d;/retrieve/d' ./rawcounts/extracted_P.count &

#### 13. 过滤表达为0的行 得到count文件。    
    awk '$2 != 0' extracted_P.count > filtered_P.txt
    awk '$2 != 0' extracted_MP.count > filtered_MP.txt
