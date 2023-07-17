import sys,argparse,json
from collections import OrderedDict

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def usages(description):
    usage = '''
Author : chenqi
Email  : 2040463170@qq.com
Date   : 2023-06-16
Description:
    %s
''' % (description)
    return usage

class Nanopore_Tool():
    def __init__(self):
        self.description="专门用于处理nanopore比对数据下的各种分析"
        self.args = self.getopt()

    def getopt(self):
        parser = argparse.ArgumentParser(
            formatter_class=HelpFormatter, description=usages(self.description))
        ##设置分组
        parser.add_argument("-v", "--version", action="version", version="%(prog)s v1.0")
        subparsers = parser.add_subparsers(dest='command')

        remap_comp_map_parser = subparsers.add_parser("remap_comp_map", help="compare remap and  map infomation")
        fastv_comp_remap_parser = subparsers.add_parser("fastv_comp_remap", help="compare fastv and  remap infomation")
        merge_each_sample_speceis_parser = subparsers.add_parser("merge_each_sample_speceis", help="merge same speceis reads")
        fastv_json_parser = subparsers.add_parser("fastv_json", help="统计fastv 相关结果")
        extract_species_parser = subparsers.add_parser("extract_species", help="从minimap2的比对结果中提取最终比对结果")
        merge_result_parser = subparsers.add_parser("merge_result", help="合并所有样本的minimap2最终比对结果")
        comp_cov_none_parser = subparsers.add_parser("comp_cov_none", help="输出比较covfilt与不过滤的结果")
        # 其他子命令的解析器remap_comp_map
        remap_comp_map_parser.add_argument("-r", "--remapfile",dest='p1_remap', type=str,  help="input remap file",required=True)
        remap_comp_map_parser.add_argument("-m", "--mappfile", dest='p1_map',type=str,  help="input map file",required=True)
        #fastv_comp_remap
        fastv_comp_remap_parser.add_argument('-f',"--fastvfile",dest='p2_fastv', type=str,  help="input fastv file",required=True)
        fastv_comp_remap_parser.add_argument('-r',"--remapfile",dest='p2_reamp', type=str,  help="input remap file",required=True)
        #merge_each_sample_speceis
        merge_each_sample_speceis_parser.add_argument('-i',"--inputfile",dest='p3_file', type=str,  help="input result csv file",required=True)
        merge_each_sample_speceis_parser.add_argument('-c',"--covthreahold",dest='p3_cov', type=float,default=0.1,  help="set cov threahold")
        #fastv_json
        fastv_json_parser.add_argument('-i',"--inputfastvjson",dest='p4_json', type=str,  help="input fastv json file",required=True)
        fastv_json_parser.add_argument('-p',"--prefix",dest='p4_prefix', type=str,  help="set prefix of outfile",default='result')
        fastv_json_parser.add_argument('-c',"--outfmt",dest='p4_outfmt', choices=['y','n'] , help="chose y or n to out diff fmt",default='y')
        #extract_species
        extract_species_parser.add_argument('-i',"--inputfilt",dest='p5_file', type=str,  help="input minimap2 result",required=True)
        #merge_result
        merge_result_parser.add_argument('-i',"--inputfilt",dest='p6_file', type=str,  help="input minimap2 result",required=True)
        merge_result_parser.add_argument('-r',"--totalreads",dest='p6_read', type=float,  help="input the sample total reads num",required=True)
        merge_result_parser.add_argument('-c',"--outfmt",dest='p6_outfmt', choices=['y','n'] , help="chose y or n to out diff fmt",default='y')
        merge_result_parser.add_argument('-p',"--prefix",dest='p6_prefix', help="set prefix of outfile",default='result')
        merge_result_parser.add_argument('-t',"--covthreahold",dest='p6_covt', help="set threahold of cov",default=0,type=float)
        #comp_cov_none
        comp_cov_none_parser.add_argument('-f',"--inputfilt",dest='p7_covfilt', type=str,  help="input covfilt result",required=True)
        comp_cov_none_parser.add_argument('-n',"--inputnone",dest='p7_nonefilt', type=str,  help="input nonefilt result",required=True)
        args = parser.parse_args()
        return args

    #输出比较covfilt与covnone的结果
    def comp_cov_none(self):
        with open(self.args.p7_covfilt,'r') as f_cov,open(self.args.p7_nonefilt,'r') as f_none:
            dict_cov = {}
            for line in f_cov:
                if line.startswith('sample'):continue
                line = line.strip().split('\t')
                if line[0] not in dict_cov:dict_cov[line[0]] = {}
                dict_cov[line[0]][line[1]] = [line[2],line[5],line[6]]
            print('sample\tspeceis\treads_filt\treads_none\tcov_filt\tcov_none\tmean_depth_filt\tmean_depth_none\tratio_filt/none')
            for line in f_none:
                if line.startswith('sample'):continue
                line = line.strip().split('\t')
                if line[0] in dict_cov:
                    if line[1] in dict_cov[line[0]]:
                        fastv_info = dict_cov[line[0]][line[1]]
                        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(line[0],line[1],fastv_info[0],line[2],fastv_info[1],line[5],fastv_info[2],line[6],int(fastv_info[0])/int(line[2])))

    #合并最终的minimap2比对结果
    def merge_result(self):
        total_reads=int(self.args.p6_read)
        y_n = self.args.p6_outfmt
        with open(self.args.p6_file,'r') as f:
            mapp_reads = 0
            dict_total = OrderedDict()
            for line in f:
                line = line.strip().split('\t')
                cov = float(line[4])
                #设置覆盖度阈值
                if cov<self.args.p6_covt:continue
                #只按照reads总数进行排名，并且只提取S16S23排名前十的内容
                type_species = line[0].split('|')
                type_species = type_species[1].split('_')[0] if len(type_species) >=2 else type_species[0]
                subspecies = line[0].split(';')[-1]
                species = '_'.join(subspecies.split('|')[:2]).upper()
                if type_species not in dict_total:dict_total[type_species] = OrderedDict()
                if species not in  dict_total[type_species]:dict_total[type_species][species]=[0,0,0,0,0]
                dict_total[type_species][species][0] += int(line[1])
                dict_total[type_species][species][1] += float(line[4])*int(line[2])
                dict_total[type_species][species][2] += float(line[5])
                dict_total[type_species][species][3] += int(line[2])
                dict_total[type_species][species][4] += 1
                dict_total[type_species][species].append(line[0])
                mapp_reads += int(line[1])

        for key,value in dict_total.items():
            value = sorted(value.items(),key = lambda a:a[1][0],reverse=True)
            for k,v in value:
                if y_n == 'y':
                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(self.args.p6_prefix,v[5],v[0],v[0]/mapp_reads,mapp_reads,v[0]/total_reads,total_reads,v[1]/v[3],v[2]/v[4]))
                elif y_n == 'n':
                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(v[5],v[0],v[0]/mapp_reads,mapp_reads,v[0]/total_reads,total_reads,v[1]/v[3],v[2]/v[4]))

    #从minimap2的结果中提取出最终的比对结果
    def extract_species(self):
        with open(self.args.p5_file,'r') as f:
            dict_total = {}
            for line in f:
                line = line.strip().split('\t')
                type_species = line[0].split('|')
                type_species = type_species[1].split('_')[0] if len(type_species) >=2 else type_species[0]
                if type_species not in dict_total:dict_total[type_species] = {}
                subspecies = line[0].split(';')[-1]
                species = '_'.join(subspecies.split('|')[:2])
                if species not in dict_total[type_species]:dict_total[type_species][species] = []
                # dict_total[type_species][species].append('\t'.join(line))
                dict_total[type_species][species].append(line)

            for key,value in dict_total.items():
                for k,v in value.items():
                    dict_total[key][k] = sorted(v,key=lambda x:int(x[1]),reverse=True)
                dict_total[key] = sorted(value.items(),key = lambda item:int(item[1][0][1]),reverse=True)
            dict_total = sorted(dict_total.items(),key=lambda x:int(x[1][0][1][0][1]),reverse=True)

            for key,value in dict_total:
                for k,v in value:
                    v = ['\t'.join([str(n) for n in i]) for i in v]
                    print('\n'.join(v))

    #提取出fastv json文件中的信息
    def fastv_json(self):
        with open(self.args.p4_json,'r') as f:
            genome_cov = json.load(f)

        total_read=int(genome_cov['summary']['before_filtering']['total_reads'])
        mapp_read = 0
        dict_total = {}
        #genome_cov_info:type:list
        genome_cov_info = genome_cov['genome_mapping_result']['genome_coverage']
        for i in genome_cov_info:
            mapp_read += int(i['reads'])
            ref_len=int(i['size'])
            total_depth=int(i['bases'])
            # print("%s\t%s\t%s\t%s\t"%(i['name'],i['reads'],i['coverage_rate'],total_depth/ref_len))
            type_species = i['name'].split('|')
            type_species = type_species[1].split('_')[0] if len(type_species) >=2 else type_species[0]
            if type_species not in dict_total:dict_total[type_species] = {}
            subspecies = i['name'].split(';')[-1]
            species = '_'.join(subspecies.split('|')[:2])
            if species not in dict_total[type_species]:dict_total[type_species][species] = [i['name'],int(i['reads']),i['coverage_rate'],total_depth/ref_len]

        prefix=self.args.p4_prefix
        y_n = self.args.p4_outfmt
        for key,value in dict_total.items():
            value = sorted(value.items(),key=lambda a:a[1][1],reverse=True)
            for k,v in value:
                if y_n == 'y':
                    print("%s\t%s\t%s\t%s\t%s\t"%(prefix,v[0],v[1],v[2],v[3]))
                elif y_n =='n':
                    print("%s\t%s\t%s\t%s\t"%(v[0],v[1],v[2],v[3]))
                ##如果要输出比对率情况可以使用以下方式 speceis\tread_num\tcov\tmean_depth\tmapp_reads_ratio_mapping\tmapp_reads\tmapp_reads_ratio_all\tall_reads
                # if y_n == 'y':
                #     print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(prefix,v[0],v[1],v[2],v[3],v[1]/mapp_read,mapp_read,v[1]/total_read,total_read))
                # elif y_n =='n':
                #     print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(v[0],v[1],v[2],v[3],v[1]/mapp_read,mapp_read,v[1]/total_read,total_read))

    #合并相同物种的reads数量并提取出top10的物种
    def merge_each_sample_speceis(self):
        dict_total = OrderedDict({'S23':OrderedDict(),'S16':OrderedDict()})
        with open(self.args.p3_file,'r') as f:
            for line in f:
                line = line.strip().split('\t')
                cov = float(line[4])
                if cov<self.args.p3_cov:continue
                #只按照reads总数进行排名，并且只提取S16S23排名前十的内容
                type_species = line[0].split('|')
                type_species = type_species[1].split('_')[0] if len(type_species) >=2 else type_species[0]
                subspecies = line[0].split(';')[-1]
                species = '_'.join(subspecies.split('|')[:2]).upper()
                if type_species in dict_total:
                    if species not in  dict_total[type_species]:dict_total[type_species][species]=[0]
                    dict_total[type_species][species][0] += int(line[1])
                    dict_total[type_species][species].append(line[0])

        for key,value in dict_total.items():
            value = sorted(value.items(),key = lambda a:a[1][0],reverse=True)
            for k,v in value[:10]:
                print(v[1]+"\t%s"%v[0])
                pass

    def remap_comp_map(self):
        print('sample\tspeceis\treads_remap\treads_map\tcov_remap\tcov_map\tmean_depth_remap\tmean_depth_map\tread_ratio_remap/map')
        with open(self.args.p1_remap,'r') as f_remap,open(self.args.p1_map,'r') as f_map:
            dict_fastv = {}
            for line in f_remap:
                if line.startswith('sample'):continue
                line = line.strip().split('\t')
                if line[0] not in dict_fastv:dict_fastv[line[0]] = {}
                dict_fastv[line[0]][line[1]] = line
            for line in f_map:
                if line.startswith('sample'):continue
                line = line.strip().split('\t')
                if line[0] in dict_fastv:
                    if line[1] in dict_fastv[line[0]]:
                        fastv_info = dict_fastv[line[0]][line[1]]
                        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(line[0],line[1],fastv_info[2],line[2],fastv_info[5],line[7],fastv_info[6],line[8], int(fastv_info[2])/int(line[2])))

    def fastv_comp_remap(self):
        print('sample\tspeceis\treads_fastv\treads_remap\tcov_fastv\tcov_remap\tmean_depth_fast\tmean_depth_remap')
        with open(self.args.p2_fastv,'r') as f_fastv,open(self.args.p2_reamp,'r') as f_remap:
            dict_fastv = {}
            for line in f_fastv:
                if line.startswith('sample'):continue
                line = line.strip().split('\t')
                if line[0] not in dict_fastv:dict_fastv[line[0]] = {}
                dict_fastv[line[0]][line[1]] = line
            for line in f_remap:
                if line.startswith('sample'):continue
                line = line.strip().split('\t')
                if line[0] in dict_fastv:
                    if line[1] in dict_fastv[line[0]]:
                        fastv_info = dict_fastv[line[0]][line[1]]
                        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(line[0],line[1],fastv_info[2],line[2],fastv_info[3],line[5],fastv_info[4],line[6]))
    
    def main(self):
        if self.args.command == "remap_comp_map":
            # 调用 remap_comp_map 相关的函数
            self.remap_comp_map()
        elif self.args.command == "fastv_comp_remap":
            # 调用 fastv_comp_remap 相关的函数
            self.fastv_comp_remap()
        elif self.args.command == "merge_each_sample_speceis":
            # 调用 merge_each_sample_speceis 相关的函数
            self.merge_each_sample_speceis()
        elif self.args.command == "fastv_json":
            # 调用 fastv_json 相关的函数
            self.fastv_json()
        elif self.args.command == "extract_species":
            # 调用 extract_species 相关的函数
            self.extract_species()
        elif self.args.command == "merge_result":
            # 调用 extract_species 相关的函数
            self.merge_result()
        elif self.args.command == "comp_cov_none":
            # 调用 comp_cov_none 相关的函数
            self.comp_cov_none()
            
tool = Nanopore_Tool()
try:
    tool.main()
except BrokenPipeError:
    sys.exit(0)
