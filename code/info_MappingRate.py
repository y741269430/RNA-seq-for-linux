import os
import re
import sys
import pandas as pd

def parse_hisat2_output(file_path):
    """解析单个Hisat2输出文件，提取关键统计指标及其百分比"""
    with open(file_path, 'r') as f:
        content = f.read()
    
    results = {'filename': os.path.basename(file_path)}
    
    # 总reads数
    total_reads = re.search(r'(\d+) reads; of these:', content)
    if total_reads:
        results['total_reads'] = int(total_reads.group(1))
    
    # 完全不对齐的reads
    concordantly_0 = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly 0 times', content)
    if concordantly_0:
        results['concordantly_0_times'] = int(concordantly_0.group(1))
        results['concordantly_0_times_%'] = float(concordantly_0.group(2))
    
    # 完全对齐一次的reads
    concordantly_1 = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time', content)
    if concordantly_1:
        results['concordantly_exactly_1_time'] = int(concordantly_1.group(1))
        results['concordantly_exactly_1_time_%'] = float(concordantly_1.group(2))
    
    # 对齐多次的reads
    concordantly_gt1 = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly >1 times', content)
    if concordantly_gt1:
        results['concordantly_>1_times'] = int(concordantly_gt1.group(1))
        results['concordantly_>1_times_%'] = float(concordantly_gt1.group(2))
    
    # 不一致对齐的reads
    discordantly_1 = re.search(r'(\d+) \(([\d.]+)%\) aligned discordantly 1 time', content)
    if discordantly_1:
        results['discordantly_1_time'] = int(discordantly_1.group(1))
        results['discordantly_1_time_%'] = float(discordantly_1.group(2))
    
    # 总对齐率
    overall_alignment = re.search(r'(\d+\.\d+)% overall alignment rate', content)
    if overall_alignment:
        results['overall_alignment_rate_%'] = float(overall_alignment.group(1))
    
    return results

def process_mapinfo_directory(directory, output_file):
    """处理mapinfo目录下的所有Hisat2输出文件并保存结果"""
    data = []
    
    for filename in os.listdir(directory):
        if filename.endswith('_map.txt'):
            file_path = os.path.join(directory, filename)
            data.append(parse_hisat2_output(file_path))
    
    if not data:
        print(f"错误: 在目录 {directory} 中没有找到任何以'_map.txt'结尾的文件")
        sys.exit(1)
    
    # 创建DataFrame并保存为CSV或Excel
    df = pd.DataFrame(data)
    
    # 重新排列列顺序
    columns_order = [
        'filename', 'total_reads',
        'concordantly_0_times', 'concordantly_0_times_%',
        'concordantly_exactly_1_time', 'concordantly_exactly_1_time_%',
        'concordantly_>1_times', 'concordantly_>1_times_%',
        'discordantly_1_time', 'discordantly_1_time_%',
        'overall_alignment_rate_%'
    ]
    
    # 只保留存在的列
    columns_order = [col for col in columns_order if col in df.columns]
    df = df[columns_order]
    
    # 根据输出文件扩展名决定保存格式
    if output_file.lower().endswith('.csv'):
        df.to_csv(output_file, index=False)
    elif output_file.lower().endswith(('.xlsx', '.xls')):
        df.to_excel(output_file, index=False)
    else:
        print("错误: 输出文件格式不支持，请使用.csv或.xlsx")
        sys.exit(1)
    
    print(f"结果已保存到 {output_file}")
    return df

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("用法: python MappingRate.py 输入目录 输出文件名.csv 或 .xlsx")
        print("示例: python MappingRate.py mapinfo alignment_results.csv")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    # 检查输入目录是否存在
    if not os.path.isdir(input_dir):
        print(f"错误: 输入目录 {input_dir} 不存在")
        sys.exit(1)
    
    # 处理文件
    try:
        df = process_mapinfo_directory(input_dir, output_file)
        print("\n整合后的数据预览:")
        print(df.head())
    except Exception as e:
        print(f"处理过程中发生错误: {str(e)}")
        sys.exit(1)