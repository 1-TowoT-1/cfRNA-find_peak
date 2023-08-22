#找到每个样本的目标基因的peak区间
import sys
import argparse
import re
import os
import numpy as np
from scipy.signal import find_peaks
import pandas as pd
from pybedtools import BedTool

doc="""
    ########
    chenzh: find_peak
    version: 1.0
    """
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='cfRNA find peak',epilog=doc)
parser.add_argument('-t','--target',help='long_trans.txt file',required=True)
parser.add_argument('-mp','--min_peak_h',help='Minimum peak value',required=True)
parser.add_argument('-n','--sample_min',help='Compliant samples',required=True)
parser.add_argument('-i','--input',required=True)
parser.add_argument('-d','--distance',help='Distance between peaks ; default=40',required=False,default=40)
parser.add_argument('-w','--wlen',help='Peak prominence ; default=80',required=False,default=80)
parser.add_argument('-m','--min_edge',help='Peak edge minimum ; default=5',required=False,default=5)

args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    exit(0)

target=args.target
min_peak_h=int(args.min_peak_h)
sample_min=int(args.sample_min)
outdir=args.input
distance=int(args.distance)
wlen=int(args.wlen)
min_edge=int(args.min_edge)



def merge_intervals(nums,chr):  #连接筛选区间
    merged_intervals = []
    start = nums[0]
    end = nums[0]
    for num in nums[1:]:
        if num == end + 1:
            end = num
        else:
            merged_intervals.append(f"{chr}\t{start}\t{end}")
            start = num
            end = num
    merged_intervals.append(f"{chr}\t{start}\t{end}")
    return merged_intervals

def reverse_complement(sequence): #负链反向互补
    complement = {'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c','N':"N","n":"n"}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement[base] for base in reverse_sequence)
    return reverse_complement_sequence


with open(target,'r') as t:
    for line in t:
        lst=line.strip().split('\t')
        gene_id=re.search(r"gene_id \"(.*?)\"",line).group(1)
        gene_name=re.search(r"gene_name \"(.*?)\"::",line).group(1)
        trans_id=re.search(r"transcript_id \"(.*?)\"",line).group(1)
        gene_direction=re.search(r"::(.)\t",line).group(1)
        indir=f'{outdir}/{gene_id}_{gene_name}/{trans_id}'

        dat=pd.read_table(f"{indir}/01.exon_depth.txt",header=0,sep='\t')

        chr=dat["chr"][0] #取一个染色体为后面bed做准备
        region_count_dic={} #记录每个位点的样本是否覆盖，键是位点pos，值为样本有覆盖计数器
        for sample_name in dat.columns[2:]:
            my_array=dat[sample_name]

            # 找到峰值的索引和峰值的高度
            peaks, _ = find_peaks(my_array,height=min_peak_h,distance=distance,wlen=wlen)

            count=0 #region计数器
            regionlist=[]
            # 打印峰值的索引和峰值的高度
            for peak_index in peaks:
                count+=1
                peak_height_max = my_array[peak_index]
                region_start=peak_index
                region_end=peak_index
                #判断区域左边，记录region_start
                for index in range(peak_index,0,-1):
                    if my_array[index]>=min_edge and my_array[index]>=peak_height_max/2:
                        region_start=index
                    else:
                        break
                
                #判断区域右边,记录region_end
                for index in range(peak_index,len(my_array),1):
                    if my_array[index]>=min_edge and my_array[index]>=peak_height_max/2:
                        region_end=index
                    else:
                        break
                chr_id=dat["chr"][peak_index]
                pos_start=dat["pos"][region_start]
                pos_end=dat["pos"][region_end]
                p=f"{chr_id}\t{pos_start}\t{pos_end}"
                regionlist.append(p)
            
            #pre.bed是单个样本可能会包含内含子的区间筛选
            with open(f'{indir}/{sample_name}_pre.bed',"w") as out: 
                regionstr='\n'.join(regionlist)
                regionbed=BedTool(regionstr,from_string=True)
                merge_bed=regionbed.sort().merge()
                for line in merge_bed:
                    start=int(line[1])
                    end=int(line[2])
                    for i in range(start,end+1):
                        if not region_count_dic.get(i):
                            region_count_dic[i]=1
                        else:
                            region_count_dic[i]+=1
                print(merge_bed,sep="\t",file=out)

        #输出符合筛选条件位点的样本个数信息的临时文件
        with open(f"{indir}/05.peak_{min_peak_h}_{sample_min}_detail.txt",'w') as p:
            filter_lst=[]
            for key,value in region_count_dic.items():    
                if value>=sample_min:
                    #筛选出符合条件的区间
                    filter_lst.append(key)
                    print(chr,key,value,sep='\t',file=p)
        
        # filter_lst=[key for key,value in region_count_dic.items() if value>=sample_min]
        sort_filter_lst=sorted(filter_lst)
        try:
            pre_bed = '\n'.join(merge_intervals(sort_filter_lst,chr))

            pre_bed=BedTool(pre_bed,from_string=True).sort()
            trans_exon=BedTool(f"{indir}/02.trans_exon_detail.txt").sort()

            with open(f"{indir}/03.peak_{min_peak_h}_{sample_min}.bed",'w') as out2:
                count=0
                for line in pre_bed:
                    count+=1
                    region=f'region_{count}'
                    line=str(line)
                    line=BedTool(line,from_string=True)
                    result=trans_exon.intersect(line)
                    for i in result:
                        i=str(i)
                        i=f'{i.strip()}\t{region}'
                        print(i,file=out2)
            
            #合并间隔外显子区间，并且筛选掉长度过短的序列，取"-"反向互补序列
            os.system("fastaFromBed -fi /project/personal/chenzh/genome/hg38_p13/hg38_p13.fa -bed %s -name $4 -fo %s"%(f"{indir}/03.peak_{min_peak_h}_{sample_min}.bed",f"{indir}/peak.fa.mid"))

            with open(f"{indir}/peak.fa.mid",'r') as p:
                region_dic={}
                for line in p:
                    line=line.strip()
                    if re.match(r"^>",line):
                        lst=line.split('::')
                        tmp=lst[0]
                        if not region_dic.get(lst[0]):
                            region_dic[lst[0]]=[gene_name,gene_id,trans_id,gene_direction,lst[1],""]
                        else:
                            region_dic[lst[0]][4]=f"{region_dic[lst[0]][4]}||{lst[1]}"
                    else:
                        if gene_direction=="-":
                            line=reverse_complement(line)
                        region_dic[tmp][5]=f"{region_dic[tmp][5]}{line}"
            
            with open(f"{indir}/04.peak_{min_peak_h}_{sample_min}.fa",'w') as out3:
                region_count=1 #记录基因名字重复的但是不同的片段
                target_region_lst=[]
                for key,value in region_dic.items():
                    if len(value[5])>=20: #将序列片段短于20bp过滤掉
                        if value[1] not in target_region_lst:
                            target_region_lst.append(value[1])
                            region_count=1
                            value[1]=f'{value[1]}::{region_count}'                         
                        else:
                            region_count=region_count+1
                            value[1]=f'{value[1]}::{region_count}'
                        length=len(value[5])
                        result=f">{value[0]}::{value[1]}::{value[2]}::{value[3]}::{value[4]}::length={length}\n{value[5]}"
                        print(result,file=out3)
        except:
            print("%s has no result"%indir)
