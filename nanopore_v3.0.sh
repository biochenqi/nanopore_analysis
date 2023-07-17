#!/bin/bash
# **********************************************************
# * Author        : chenqi
# * Email         : chenqi@gooalgene.com
# * Create time   : 2023-05-17 13:44
# * Last modified : 2023-05-17 13:44
# * Filename      : nanopore.sh
# * Description   : 
# **********************************************************
usage() {
  echo "Usage: $0 [-a 参数A] [-b 参数B] [-h]"
  echo "  -o outdir  : 输出目录"
  echo "  -i inputdir  : 输入目录（注:该目录下数据后缀为fastq）"
  echo "  -e excute  : 选择是否直接执行[y:执行|n：不执行;默认：n]"
  echo "  -n thread  : 设置线程数量[默认：1]"
  echo "  -c covthreahold  : 设置阈值用于过滤掉覆盖度过低的结果[默认：0.1]"
  echo "  -x excutenonecov  : 选择是否执行不过滤覆盖度的top10重新比对分析[y:执行|n:不执行；默认：y]"
  echo "  -h        : 显示帮助信息"
}

arg_c=0.1
arg_n=1
arg_e='n'
arg_x='y'
while getopts ":o:i:c:n:e:x:h" opt; do
  case $opt in
    o)
      arg_o="$OPTARG"
      ;;
    i)
      arg_i="$OPTARG"
      ;;
    c)
      arg_c="$OPTARG"
      ;;
    n)
      arg_n="$OPTARG"
      ;;
    e)
      arg_e="$OPTARG"
      ;;
    x)
      arg_e="$OPTARG"
      ;;
    h)
      usage
      exit 0
      ;;
    \?)
      echo "无效选项: -$OPTARG" >&2
      usage
      exit 1
      ;;
    :)
      echo "选项 -$OPTARG 需要参数." >&2
      usage
      exit 1
      ;;
  esac
done

if [ ! -n "$arg_o" ] || [ ! -n "$arg_i" ];then
    usage
    exit 0
fi

echo start at time `date +%F'  '%H:%M:%S`
#software
nanoplot='/home/liyingqiang/anaconda3/envs/nanoplot/bin/NanoPlot'
minimap2='/haplox/haprs/liyq/soft/minimap2-2.17_x64-linux/minimap2'
samtools='/haplox/haprs/chenqi/software/anaconda3/bin/samtools'
seqkit='/haplox/haprs/chenqi/software/anaconda3/bin/seqkit'
fastv='/haplox/haprs/chenqi/software/fastv/fastv'
nanofilt="/home/liyingqiang/anaconda3/envs/nanoplot/bin/NanoFilt"
taxonkit='/haplox/haprs/chenqi/software/anaconda3/bin/taxonkit'

#bin-script
bin_dir='/haplox/haprs/chenqi/project/Nanopore/bin'
total_script="$bin_dir/nanopore_tools.py"

#database
# bwahuman='/haplox/haprs/liyq/bin/VIcall/references/references/hg19_wesplus.fa'
micro_db='/haplox/haprs/chenqi/test/test_mNGS/nanopore/ref/com_rcs2_s16_s23.fa'
# gb_to_txid='/haplox/haprs/chenqi/database/nt/nucl_gb.accession2taxid'
# nt='/haplox/haprs/chenqi/database/nt/nt/nt'
nc_gb='/haplox/haprs/chenqi/database/nt/nucl_gb.accession2taxid'

##生成目录结构
date_now=`date +%Y%m%d`
mkdir -p $arg_o
mkdir -p $arg_o/01.nanoplot
mkdir -p $arg_o/02.mapp
mkdir -p $arg_o/03.depth
mkdir -p $arg_o/04.stats
mkdir -p $arg_o/Result_${date_now}
mkdir -p $arg_o/Result_${date_now}/allsample
mkdir -p $arg_o/Shell
shell_dir="$arg_o/Shell"

#generate-script
nanoplot_sc="${shell_dir}/01.nanoplot.sh"
mapp_sc="${shell_dir}/02.mapp.sh"
depth_sc="${shell_dir}/03.depth.sh"
stats_sc="${shell_dir}/04.stats.sh"
fastv_remap="${shell_dir}/05.fastv_remap.sh"
stats_remap_sc="${shell_dir}/06.stats_remap.sh"
result_sc="$arg_o/Shell/result.sh"

echo "#01.nanoplot" >$nanoplot_sc
echo "#02.mapp" >$mapp_sc
echo "#03.depth" >$depth_sc
echo "#04.stats" >$stats_sc
echo "#05.fastv_remap" >$fastv_remap
echo "#06.stats_remap" >$stats_remap_sc
echo "#result" >$result_sc

pwd_dir=`pwd`

#执行各个单个的shell脚本时需要先进入工作目录
echo "cd $pwd_dir" >>$nanoplot_sc
echo "cd $pwd_dir" >>$mapp_sc
echo "cd $pwd_dir" >>$depth_sc
echo "cd $pwd_dir" >>$stats_sc
echo "cd $pwd_dir" >>$fastv_remap
echo "cd $pwd_dir" >>$stats_remap_sc
echo "cd $pwd_dir" >>$result_sc

#对nanoplot结果进行一个汇总
echo "echo -e \"sample\tMean_len\tMean_qua\tMedian_len\tMedian_qua\tNum_reads\tRead_N50\tTotal_bases\t>Q5\t>Q7\t>Q10\t>Q12\t>Q15\th1\th2\th3\th4\th5\tl1\tl2\tl3\tl4\tl5\t\" >$arg_o/Result_${date_now}/01.nanoplot.stat.xls " >> $result_sc
echo "echo -e \"sample\tMean_len\tMean_qua\tMedian_len\tMedian_qua\tNum_reads\tRead_N50\tTotal_bases\t>Q5\t>Q7\t>Q10\t>Q12\t>Q15\th1\th2\th3\th4\th5\tl1\tl2\tl3\tl4\tl5\t\" >$arg_o/Result_${date_now}/01.nanoplot.filt.stat.xls " >> $result_sc
echo "echo -e \"sample\tspecies\tmapp_reads\tmapp_reads_ratio_mapping\tmapping_reads\tmapp_reads_ratio_all\tall_reads\tcov\tmean_depth\" >$arg_o/Result_${date_now}/02.mapp.stat.xls " >>$result_sc
echo "echo -e \"sample\tspecies\tmapp_reads\tmapp_reads_ratio_mapping\tmapping_reads\tmapp_reads_ratio_all\tall_reads\tcov\tmean_depth\" >$arg_o/Result_${date_now}/02.mapp.stat.covfilt.xls " >>$result_sc
echo "echo -e \"sample\tspecies\tmapp_reads\tref_len\ttotal_depth\tcov mean_depth\" >$arg_o/Result_${date_now}/03.remap.stat.xls " >>$result_sc
echo "echo -e \"sample\tspecies\tmapp_reads\tref_len\ttotal_depth\tcov mean_depth\" >$arg_o/Result_${date_now}/03.remap.none.stat.xls " >>$result_sc
echo "echo -e \"sample\tspeceis\tread_num\tcov\tmean_depth\"  >$arg_o/Result_${date_now}/03.fastv.stat.xls " >>$result_sc

#各步骤脚本生成
ls $arg_i/*fastq|while read line;
do
    prefix=`python -c "import sys;print(sys.argv[1].split('/')[-1][:-6])" ${line}`
    echo "##################${prefix}####################" >> $nanoplot_sc
    echo "##################${prefix}####################" >> $mapp_sc
    echo "##################${prefix}####################" >> $depth_sc
    echo "##################${prefix}####################" >> $stats_sc
    echo "##################${prefix}####################" >> $fastv_remap
    echo "##################${prefix}####################" >> $stats_remap_sc
    echo "##################${prefix}####################" >> $result_sc

    #0.1nanoplot
    mkdir -p $arg_o/01.nanoplot/${prefix}
    nanoplot_dir="$arg_o/01.nanoplot/${prefix}"
    filt_nanoplot_dir="$arg_o/01.nanoplot/${prefix}/filt_nanoplot"
    mkdir -p $filt_nanoplot_dir
    echo $nanoplot -t $arg_n --fastq ${line} --plots hex dot -o $nanoplot_dir >> $nanoplot_sc
    #质控去除掉碱基序列小于5，序列长度低于200的reads
    echo "$nanofilt -q 5 -l 200 ${line} > ${nanoplot_dir}/${prefix}.filt.fq" >> $nanoplot_sc
    #统计质控后的序列情况
    echo $nanoplot -t $arg_n --fastq ${nanoplot_dir}/${prefix}.filt.fq --plots hex dot -o $filt_nanoplot_dir >> $nanoplot_sc

    #0.2mapp
    mkdir -p $arg_o/02.mapp/${prefix}
    mapp_dir="$arg_o/02.mapp/${prefix}"
    #无过滤结果
    # echo "$minimap2 -ax map-ont -N 0 -t $arg_n $micro_db ${line} -o ${mapp_dir}/${prefix}.sam " >>$mapp_sc
    #原始序列过滤比对结果
    echo "$minimap2 -ax map-ont -N 0 -t $arg_n $micro_db ${nanoplot_dir}/${prefix}.filt.fq -o ${mapp_dir}/${prefix}_s23s16.filt.sam " >>$mapp_sc
    #直接提取出所有reads的唯一比对结果，有多重比对结果的情况首选flag==2,如果没有则选择遇到的第一个结果
    echo "grep -v '^@' ${mapp_dir}/${prefix}_s23s16.filt.sam |cut -f 1 |sort |uniq -d  >${mapp_dir}/list_id.txt" >>$mapp_sc
    echo "awk -F '\t' 'ARGIND==1{info[\$1]=0}ARGIND==2{if(\$1 in info){if(\$2==0){info[\$1]=\$0}else{if(info[\$1]==0){info[\$1]=\$0}}}else{print \$0}}END{for(key in info)print info[key]}' ${mapp_dir}/list_id.txt ${mapp_dir}/${prefix}_s23s16.filt.sam |$samtools view -bS -@ $arg_n - |$samtools sort -@ $arg_n -o ${mapp_dir}/${prefix}.sort.bam" >>$mapp_sc
    #比对nt库结果
    # echo "$minimap2 -ax map-ont -N 0 -t $arg_n $nt ${nanoplot_dir}/${prefix}.filt.fq -o ${mapp_dir}/${prefix}_nt.filt.sam " >>$mapp_sc
    # echo "grep -v '^@' ${mapp_dir}/${prefix}_nt.filt.sam |cut -f 1 |sort |uniq -d  >${mapp_dir}/list_nt_id.txt" >>$mapp_sc
    # echo "awk -F '\t' 'ARGIND==1{info[\$1]=0}ARGIND==2{if(\$1 in info){if(\$2==0){info[\$1]=\$0}else{if(info[\$1]==0){info[\$1]=\$0}}}else{print \$0}}END{for(key in info)print info[key]}' ${mapp_dir}/list_nt_id.txt ${mapp_dir}/${prefix}_nt.filt.sam >${mapp_dir}/${prefix}_nt.redup.sam" >>$mapp_sc
    # echo "awk -F '\t' '{if(\$0!~/^@/ && \$2!=4){if(!(\$3 in info)){info[\$3]=0};info[\$3]++}}END{for(key in info){print key\"\t\"info[key]}}' ${mapp_dir}/${prefix}_nt.redup.sam |sort -nrk 2 >${mapp_dir}/test.result" >>$mapp_sc
    # echo "cut -f 1 ${mapp_dir}/test.result  > ${mapp_dir}/gb_name.txt" >>$mapp_sc
    # echo "grep -w -f ${mapp_dir}/gb_name.txt $nc_gb >${mapp_dir}/txid.txt" >>$mapp_sc
    # echo "cut -f 3 ${mapp_dir}/txid.txt|$taxonkit lineage -L -n - |awk -F '\t' 'ARGIND==1{gb_num[\$1]=\$2}ARGIND==2{gb_txid[\$2]=\$3}ARGIND==3{txid_sc[\$1]=\$2}END{for(key in gb_num){print txid_sc[gb_txid[key]]\"\t\"gb_num[key]}}' ${mapp_dir}/test.result ${mapp_dir}/txid.txt - |sort -t \$'\t' -nrk 2 > ${mapp_dir}/${prefix}_nt.species.result" >>$mapp_sc

    #0.3depth
    mkdir -p $arg_o/03.depth/${prefix}
    depth_dir="$arg_o/03.depth/${prefix}"
    echo "$samtools depth -a ${mapp_dir}/${prefix}.sort.bam >${depth_dir}/${prefix}.depth" >> $depth_sc

    #0.4stats
    mkdir -p $arg_o/04.stats/${prefix}
    stats_dir="$arg_o/04.stats/${prefix}"
    #计算每个物种的参考序列长度，总深度, 覆盖度，平均深度 pos:储存每个物种的碱基序列长度；count:储存每个物种的覆盖到的位点总数；info：储存每个物种深度之和
    echo "awk -F '\t' '{pos[\$1]=\$2;if(!(\$1 in count)){count[\$1]=0};if(\$3!=0)count[\$1]++;if(!(\$1 in info)){info[\$1]=0};info[\$1]=info[\$1] + \$3}END{for(key in pos){print key\"\t\"pos[key]\"\t\"info[key]\"\t\"count[key]/pos[key]\"\t\"info[key]/pos[key]}}' ${depth_dir}/${prefix}.depth > ${stats_dir}/result.txt" >>$stats_sc
    #统计比对reads条数
    echo "$samtools view ${mapp_dir}/${prefix}.sort.bam|awk -F '\t' '{if(\$2!=4){if(!(\$3 in info)){info[\$3]=0};info[\$3]++}}END{for(key in info){print key\"\t\"info[key]}}' >${stats_dir}/species_result.txt" >>$stats_sc
    echo "awk -F '\t' 'ARGIND==1{info[\$1]=\$2}ARGIND==2{split(\$1,name,\";\");if(\$1 in info)print \$1\"\t\"info[\$1]\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5}'  ${stats_dir}/species_result.txt ${stats_dir}/result.txt >${stats_dir}/all_stats.txt" >>$stats_sc
    echo "python $total_script extract_species -i ${stats_dir}/all_stats.txt >${stats_dir}/${prefix}.result.csv" >> $stats_sc

    #0.5fastv and remap
    mkdir -p $arg_o/05.fastv_remap/${prefix}
    fastv_remap_dir="$arg_o/05.fastv_remap/${prefix}"
    if [ $arg_x == "y" ];then
    ##不进行过滤的结果
    #对结果文件result.csv先将S23和s16内同一物种的序列合并，然后按照reads数量进行排序，并且最终只选取S23S16各排名前十的序列,保留这所有结果，不进行过滤
    echo "python $total_script merge_each_sample_speceis -i ${stats_dir}/${prefix}.result.csv -c 0 >${fastv_remap_dir}/extract.txt" >> $fastv_remap
    #通过seqkit提取出16s23s中对应的序列
    echo "cut -f 1 ${fastv_remap_dir}/extract.txt|while read line;do $seqkit grep -i -n $micro_db -p \${line};done |awk '{if(\$0~/^>/){print \$0;next};gsub(\"[Uu]\",\"T\");print \$0}' >$fastv_remap_dir/${prefix}_ref.fa" >> $fastv_remap
    #fastv进行
    echo "awk 'NR%4==2{gsub(\"[Uu]\",\"T\");print \$0}NR%4!=2' ${nanoplot_dir}/${prefix}.filt.fq  >${fastv_remap_dir}/${prefix}_transe.filt.fq" >> $fastv_remap
    echo "$fastv -i ${fastv_remap_dir}/${prefix}_transe.filt.fq -g $fastv_remap_dir/${prefix}_ref.fa -h ${fastv_remap_dir}/${prefix}.html -j ${fastv_remap_dir}/${prefix}.json" >> $fastv_remap
    #remap进行
    echo "$minimap2 -ax map-ont -N 0 -t $arg_n $fastv_remap_dir/${prefix}_ref.fa ${fastv_remap_dir}/${prefix}_transe.filt.fq -o ${fastv_remap_dir}/${prefix}.filt.sam " >>$fastv_remap
    #直接提取出所有reads的唯一比对结果，有多重比对结果的情况首选flag==2,如果没有则选择遇到的第一个结果
    echo "grep -v '^@' ${fastv_remap_dir}/${prefix}.filt.sam |cut -f 1 |sort |uniq -d  >${fastv_remap_dir}/list_id.txt" >>$fastv_remap
    echo "awk -F '\t' 'ARGIND == 1{info[\$1]=0}ARGIND == 2{if(\$1 in info){if(\$2==0){info[\$1]=\$0}else{if(info[\$1]==0){info[\$1]=\$0}}}else{print \$0}}END{for(key in info)print info[key]}' ${fastv_remap_dir}/list_id.txt ${fastv_remap_dir}/${prefix}.filt.sam |$samtools view -bS -@ $2 - |$samtools sort -@ $2 -o ${fastv_remap_dir}/${prefix}.sort.bam" >>$fastv_remap
    fi

    ##进行过滤的结果
    #对结果文件result.csv先将S23和s16内同一物种的序列合并，然后按照reads数量进行排序，并且最终只选取S23S16各排名前十的序列,设置过滤标准，覆盖度没有达到0.1的都去除掉
    echo "python $total_script merge_each_sample_speceis -i ${stats_dir}/${prefix}.result.csv -c $arg_c >${fastv_remap_dir}/extract.covfilt.txt" >> $fastv_remap
    #通过seqkit提取出16s23s中对应的序列
    echo "cut -f 1 ${fastv_remap_dir}/extract.covfilt.txt|while read line;do $seqkit grep -i -n $micro_db -p \${line};done |awk '{if(\$0~/^>/){print \$0;next};gsub(\"[Uu]\",\"T\");print \$0}' >$fastv_remap_dir/${prefix}_ref.covfilt.fa" >> $fastv_remap
    #fastv进行
    echo "awk 'NR%4==2{gsub(\"[Uu]\",\"T\");print \$0}NR%4!=2' ${nanoplot_dir}/${prefix}.filt.fq  >${fastv_remap_dir}/${prefix}_transe.covfilt.filt.fq" >> $fastv_remap
    echo "$fastv -i ${fastv_remap_dir}/${prefix}_transe.filt.fq -g $fastv_remap_dir/${prefix}_ref.fa -h ${fastv_remap_dir}/${prefix}.html -j ${fastv_remap_dir}/${prefix}.json" >> $fastv_remap
    #remap进行
    echo "$minimap2 -ax map-ont -N 0 -t $arg_n $fastv_remap_dir/${prefix}_ref.covfilt.fa ${fastv_remap_dir}/${prefix}_transe.covfilt.filt.fq -o ${fastv_remap_dir}/${prefix}.covfilt.filt.sam " >>$fastv_remap
    #直接提取出所有reads的唯一比对结果，有多重比对结果的情况首选flag==2,如果没有则选择遇到的第一个结果
    echo "grep -v '^@' ${fastv_remap_dir}/${prefix}.covfilt.filt.sam |cut -f 1 |sort |uniq -d  >${fastv_remap_dir}/list_id.covfilt.txt" >>$fastv_remap
    echo "awk -F '\t' 'ARGIND == 1{info[\$1]=0}ARGIND == 2{if(\$1 in info){if(\$2==0){info[\$1]=\$0}else{if(info[\$1]==0){info[\$1]=\$0}}}else{print \$0}}END{for(key in info)print info[key]}' ${fastv_remap_dir}/list_id.covfilt.txt ${fastv_remap_dir}/${prefix}.covfilt.filt.sam |$samtools view -bS -@ $2 - |$samtools sort -@ $2 -o ${fastv_remap_dir}/${prefix}.covfilt.sort.bam" >>$fastv_remap

    #0.6stats_remap
    mkdir -p $arg_o/06.stats_remap/${prefix}
    stats_remap_dir="$arg_o/06.stats_remap/${prefix}"
    if [ $arg_x == "y" ];then
    ##不进行过滤的结果
    echo "$samtools depth -a ${fastv_remap_dir}/${prefix}.sort.bam >${stats_remap_dir}/${prefix}.depth" >> $stats_remap_sc
    echo "awk -F '\t' '{pos[\$1]=\$2;if(!(\$1 in count)){count[\$1]=0};if(\$3!=0)count[\$1]++;if(!(\$1 in info)){info[\$1]=0};info[\$1]=info[\$1] + \$3}END{for(key in pos){print key\"\t\"pos[key]\"\t\"info[key]\"\t\"count[key]/pos[key]\"\t\"info[key]/pos[key]}}' ${stats_remap_dir}/${prefix}.depth > ${stats_remap_dir}/result.txt" >>$stats_remap_sc
    #统计比对reads条数
    echo "$samtools view ${fastv_remap_dir}/${prefix}.sort.bam|awk -F '\t' '{if(\$2!=4){if(!(\$3 in info)){info[\$3]=0};info[\$3]++}}END{for(key in info){print key\"\t\"info[key]}}' >${stats_remap_dir}/species_result.txt" >>$stats_remap_sc
    echo "awk -F '\t' 'ARGIND==1{info[\$1]=\$2}ARGIND==2{split(\$1,name,\";\");if(\$1 in info)print \$1\"\t\"info[\$1]\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5}'  ${stats_remap_dir}/species_result.txt ${stats_remap_dir}/result.txt >${stats_remap_dir}/all_stats.txt" >>$stats_remap_sc
    echo "python $total_script extract_species -i ${stats_remap_dir}/all_stats.txt >${stats_remap_dir}/${prefix}.result.csv" >> $stats_remap_sc
    fi

    ##进行过滤的结果
    echo "$samtools depth -a ${fastv_remap_dir}/${prefix}.covfilt.sort.bam >${stats_remap_dir}/${prefix}.covfilt.depth" >> $stats_remap_sc
    echo "awk -F '\t' '{pos[\$1]=\$2;if(!(\$1 in count)){count[\$1]=0};if(\$3!=0)count[\$1]++;if(!(\$1 in info)){info[\$1]=0};info[\$1]=info[\$1] + \$3}END{for(key in pos){print key\"\t\"pos[key]\"\t\"info[key]\"\t\"count[key]/pos[key]\"\t\"info[key]/pos[key]}}' ${stats_remap_dir}/${prefix}.covfilt.depth > ${stats_remap_dir}/result.covfilt.txt" >>$stats_remap_sc
    #统计比对reads条数
    echo "$samtools view ${fastv_remap_dir}/${prefix}.covfilt.sort.bam|awk -F '\t' '{if(\$2!=4){if(!(\$3 in info)){info[\$3]=0};info[\$3]++}}END{for(key in info){print key\"\t\"info[key]}}' >${stats_remap_dir}/species_result.covfilt.txt" >>$stats_remap_sc
    echo "awk -F '\t' 'ARGIND==1{info[\$1]=\$2}ARGIND==2{split(\$1,name,\";\");if(\$1 in info)print \$1\"\t\"info[\$1]\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5}'  ${stats_remap_dir}/species_result.covfilt.txt ${stats_remap_dir}/result.covfilt.txt >${stats_remap_dir}/all_stats.covfilt.txt" >>$stats_remap_sc
    echo "python $total_script extract_species -i ${stats_remap_dir}/all_stats.covfilt.txt >${stats_remap_dir}/${prefix}.covfilt.result.csv" >> $stats_remap_sc

    #result
    #nanoplot信息统计
    echo "awk -v name=$prefix -F ':' 'BEGIN{}{if(\$0~/^General/)next;if(\$0~/:/){gsub(/^[ \t]+|[ \t]+\$/, \"\", \$2);name=name\"\t\"\$2}}END{print name\"\t\"}' ${nanoplot_dir}/NanoStats.txt >>$arg_o/Result_${date_now}/01.nanoplot.stat.xls" >>$result_sc
    echo "awk -v name=$prefix -F ':' 'BEGIN{}{if(\$0~/^General/)next;if(\$0~/:/){gsub(/^[ \t]+|[ \t]+\$/, \"\", \$2);name=name\"\t\"\$2}}END{print name\"\t\"}' ${filt_nanoplot_dir}/NanoStats.txt >>$arg_o/Result_${date_now}/01.nanoplot.filt.stat.xls" >>$result_sc
    #初次比对16s23s信息统计
    # echo "awk -v name=$prefix '{print name\"\t\"\$0}' ${stats_dir}/${prefix}.result.csv >> $arg_o/Result_${date_now}/02.mapp.stat.xls" >>$result_sc
    #获取不过滤合并物种后的结果
    echo "total_reads=\`grep -w \"Number of reads\" $filt_nanoplot_dir/NanoStats.txt |cut -d \":\" -f 2|awk '{sub(\",\",\"\",\$1);print \$1}'\`" >>$result_sc

    echo "python $total_script merge_result -i ${stats_dir}/${prefix}.result.csv -r \$total_reads -c y -p $prefix >> $arg_o/Result_${date_now}/02.mapp.stat.xls" >>$result_sc
    echo "echo -e \"species\tmapp_reads\tmapp_reads_ratio_mapping\tmapping_reads\tmapp_reads_ratio_all\tall_reads\tcov\tmean_depth\" > $arg_o/Result_${date_now}/allsample/$prefix.mapp.xls" >>$result_sc
    echo "python $total_script merge_result -i ${stats_dir}/${prefix}.result.csv -r \$total_reads -c n -p $prefix >> $arg_o/Result_${date_now}/allsample/$prefix.mapp.xls" >>$result_sc
    #获取过滤覆盖度后物种合并的结果
    echo "python $total_script merge_result -i ${stats_dir}/${prefix}.result.csv -r \$total_reads -c y -p $prefix -t 0.1 >> $arg_o/Result_${date_now}/02.mapp.stat.covfilt.xls" >>$result_sc
    echo "echo -e \"species\tmapp_reads\tmapp_reads_ratio_mapping\tmapping_reads\tmapp_reads_ratio_all\tall_reads\tcov\tmean_depth\" > $arg_o/Result_${date_now}/allsample/$prefix.mapp.covfilt.xls" >>$result_sc
    echo "python $total_script merge_result -i ${stats_dir}/${prefix}.result.csv -r \$total_reads -c n -p $prefix -t 0.1 >> $arg_o/Result_${date_now}/allsample/$prefix.mapp.covfilt.xls" >>$result_sc
    ##提取出最优比对top10的16s23s参考序列后重新比对
    #fastv结果
    echo "python $total_script fastv_json -i ${fastv_remap_dir}/${prefix}.json -p ${prefix} -c y >> $arg_o/Result_${date_now}/03.fastv.stat.xls" >>$result_sc
    # echo "echo -e \"speceis\tread_num\tcov\tmean_depth\" > $arg_o/Result_${date_now}/allsample/${prefix}.fastv.xls" >>$result_sc
    # echo "python $total_script fastv_json -i ${fastv_remap_dir}/${prefix}.json -p ${prefix} -c n > $arg_o/Result_${date_now}/allsample/${prefix}.fastv.xls" >>$result_sc
    echo "cp ${fastv_remap_dir}/${prefix}.html $arg_o/Result_${date_now}/allsample/" >>$result_sc
    #remap结果
    echo "awk -v sap=$prefix '{print sap\"\t\"\$0}' ${stats_remap_dir}/${prefix}.result.csv >> $arg_o/Result_${date_now}/03.remap.none.stat.xls" >> $result_sc
    echo "awk -v sap=$prefix '{print sap\"\t\"\$0}' ${stats_remap_dir}/${prefix}.covfilt.result.csv >> $arg_o/Result_${date_now}/03.remap.stat.xls" >> $result_sc
    echo "echo -e \"species\tread_num\tref_len\tdepth_total\tcov\tmean_depth\" > $arg_o/Result_${date_now}/allsample/${prefix}.remap.xls" >> $result_sc
    echo "cat ${stats_remap_dir}/${prefix}.covfilt.result.csv >> $arg_o/Result_${date_now}/allsample/${prefix}.remap.xls" >> $result_sc
    

    echo -e "\n\n" >> $nanoplot_sc
    echo -e "\n\n" >> $mapp_sc
    echo -e "\n\n" >> $depth_sc
    echo -e "\n\n" >> $stats_sc
    echo -e "\n\n" >> $fastv_remap
    echo -e "\n\n" >> $stats_remap_sc
    echo -e "\n\n" >> $result_sc
done
#fastv比较remap结果、remap比较map结果、covfilt 比较nonefilt结果
echo "python $total_script remap_comp_map -r $arg_o/Result_${date_now}/03.remap.stat.xls -m $arg_o/Result_${date_now}/02.mapp.stat.xls > $arg_o/Result_${date_now}/remap_comp_mapp.xls">> $result_sc
echo "python $total_script fastv_comp_remap -f $arg_o/Result_${date_now}/03.fastv.stat.xls -r $arg_o/Result_${date_now}/03.remap.stat.xls > $arg_o/Result_${date_now}/fastv_comp_remap.xls">> $result_sc
echo "python $total_script comp_cov_none -f $arg_o/Result_${date_now}/03.remap.stat.xls -n $arg_o/Result_${date_now}/03.remap.none.stat.xls > $arg_o/Result_${date_now}/filt_comp_none.xls">> $result_sc
echo "rm $arg_o/Result_${date_now}/03.remap.none.stat.xls" >> $result_sc

if [ -n $arg_e ] && [ $arg_e == "y" ];then
    sh $nanoplot_sc
    sh $mapp_sc
    sh $depth_sc
    sh $stats_sc
    sh $fastv_remap
    sh $stats_remap_sc
    sh $result_sc
fi

#汇总所有脚本用于运行结果
pipline_sc="$arg_o/Shell/pipline.sh"
echo "cd $pwd_dir" >$pipline_sc
echo "sh $nanoplot_sc" >> $pipline_sc
echo "sh $mapp_sc" >> $pipline_sc
echo "sh $depth_sc" >> $pipline_sc
echo "sh $stats_sc" >> $pipline_sc
echo "sh $fastv_remap" >> $pipline_sc
echo "sh $stats_remap_sc" >> $pipline_sc
echo "sh $result_sc" >> $pipline_sc

echo finish at time `date +%F'  '%H:%M:%S`

