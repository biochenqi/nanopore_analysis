# nanopore_analysis
## pipline

**nanopore_v3.0.sh** 是一个用于分析三代测序nanopore下机数据用于微生物检测的流程.

## Installation
Install **hurst** module with      
`pip install -e git+https://github.com/Mottl/hurst#egg=hurst`

## nanopore_v3.0.sh Usage
```
Usage: nanopore_v3.0.sh [-a 参数A] [-b 参数B] [-h]
  -o outdir  : 输出目录
  -i inputdir  : 输入目录（注:该目录下数据后缀为fastq）
  -e excute  : 选择是否直接执行[y:执行|n：不执行;默认：n]
  -n thread  : 设置线程数量[默认：1]
  -c covthreahold  : 设置阈值用于过滤掉覆盖度过低的结果[默认：0.1]
  -x excutenonecov  : 选择是否执行不过滤覆盖度的top10重新比对分析[y:执行|n:不执行；默认：y]
  -h        : 显示帮助信息 
```

## nanopore_tools.py
```
usage: nanopore_tools.py [-h] [-v]
                         {remap_comp_map,fastv_comp_remap,merge_each_sample_speceis,fastv_json,extract_species,merge_result,comp_cov_none}
                         ...

Author : chenqi
Email  : 2040463170@qq.com
Date   : 2023-06-16
Description:
    专门用于处理nanopore比对数据下的各种分析

positional arguments:
  {remap_comp_map,fastv_comp_remap,merge_each_sample_speceis,fastv_json,extract_species,merge_result,comp_cov_none}
    remap_comp_map      compare remap and map infomation
    fastv_comp_remap    compare fastv and remap infomation
    merge_each_sample_speceis
                        merge same speceis reads
    fastv_json          统计fastv 相关结果
    extract_species     从minimap2的比对结果中提取最终比对结果
    merge_result        合并所有样本的minimap2最终比对结果
    comp_cov_none       输出比较covfilt与不过滤的结果

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```
