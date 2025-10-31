# 检查gtf文件结构

```bash
vim gtf_type_count.sh
```
写入以下内容
```bash
#!/bin/bash

# 默认统计第3列
COL=3

# 判断用户是否提供了列号参数（即总共传了2个参数：文件路径 + 列号）
if [ $# -ge 2 ]; then
    COL="$2"
    # 简单检查列号是否是数字
    if ! [[ "$COL" =~ ^[0-9]+$ ]]; then
        echo "错误：列号必须是数字，比如 3"
        exit 1
    fi
fi

# 检查是否提供了GTF文件路径（至少传1个参数）
if [ $# -lt 1 ]; then
    echo "用法: $0 <GTF文件路径> [列号，默认为3]"
    echo "示例: $0 file.gtf       # 统计第3列"
    echo "示例: $0 file.gtf 5     # 统计第5列"
    exit 1
fi

GTF_FILE="$1"

# 检查文件是否存在
if [ ! -f "$GTF_FILE" ]; then
    echo "错误：文件 '$GTF_FILE' 不存在！"
    exit 1
fi

# 提取指定列，统计并彩色输出
awk -F'\t' -v col="$COL" '{print $col}' "$GTF_FILE" | sort | uniq -c | awk '{printf "\033[1;32m%s\033[0m\t%d\n", $2, $1}'
```
运行，默认为第三列，可选其他列，统计结果（类似R里面的table功能）
```bash
bash gtf_type_count.sh xxx.gtf 3
```
