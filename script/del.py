import sys
import re
import os


merge_exon=sys.argv[1]
bams=sys.argv[2]  #文件形式输入，真实bam文件路径，以空格分割
target_gene=sys.argv[3]
sample_name=sys.argv[4] #字符串样本名，以","分割
gtf=sys.argv[5]
output=sys.argv[6]
threads=int(sys.argv[7])

with open(bams,'r') as b:
    bam_files=''
    for line in b:
        bam_files=bam_files + " " + line.strip()

sample_name='\t'.join(sample_name.split(','))

os.system("samtools depth -@ %d -a -b %s %s > %s/merge_depth_coverage.txt"%(threads,merge_exon,bam_files,output))

with open(target_gene,'r') as f1:
    target_gene_dic={}
    for line in f1:
        lst=line.strip().split('\t')
        target_gene_dic[lst[0]]=lst[1]


with open(gtf,'r') as f2:
    chr_dic={} #记录这个染色体中的目标基因的每个exon信息
    trans_exon_dic={} #记录目标转录本的外显子区间
    for line in f2:
        if re.search(r"^#",line):
            next
        else:
            lst=line.strip().split('\t')
            if lst[2]=="exon":
                gene_id=re.search(r"gene_id \"(.*?)\";",lst[-1]).group(1)
                gene_name=re.search(r"gene_name \"(.*?)\";",lst[-1]).group(1)
                trans_id=re.search(r"transcript_id \"(.*?)\";",lst[-1]).group(1)
                key_name=f'{gene_id}_{gene_name}/{trans_id}'
                if gene_id in target_gene_dic.keys():
                    line=line.strip()
                    value=f'{lst[0]}\t{lst[3]}\t{lst[4]}'
                    if not chr_dic.get(lst[0]):
                        chr_dic[lst[0]]=[line]
                    else:
                        chr_dic[lst[0]].append(line)
                    if not trans_exon_dic.get(key_name):    
                        trans_exon_dic[key_name]=[value]
                    else:
                        trans_exon_dic[key_name].append(value)
                                       


merge_depth_coverage=f'{output}/merge_depth_coverage.txt'
with open(merge_depth_coverage,'r') as f3:
    trans_dic={} #记录目标gene_trans为键，值为exon行记录
    for line in f3:
        line=line.strip()
        lst=line.split('\t')
        for exon_detail in chr_dic[lst[0]]:          
            exon_row_detail_lst=exon_detail.split('\t')
            gene_id=re.search(r"gene_id \"(.*?)\";",exon_row_detail_lst[-1]).group(1)
            gene_name=re.search(r"gene_name \"(.*?)\";",exon_row_detail_lst[-1]).group(1)
            trans_id=re.search(r"transcript_id \"(.*?)\";",exon_row_detail_lst[-1]).group(1)
            key_name=f'{gene_id}_{gene_name}/{trans_id}'
            if int(lst[1])>=int(exon_row_detail_lst[3]) and int(lst[1])<=int(exon_row_detail_lst[4]):
                if not trans_dic.get(key_name):
                    trans_dic[key_name]=[line]
                else:
                    trans_dic[key_name].append(line)

for key,value in trans_dic.items():
    os.system("mkdir -p %s/%s"%(output,key))
    result_value='\n'.join(value)
    trans_exon_value='\n'.join(trans_exon_dic[key])
    exon_depth_file=f'{output}/{key}/01.exon_depth.txt'
    trans_exon=f'{output}/{key}/02.trans_exon_detail.txt'
    with open(exon_depth_file,'w') as out,open(trans_exon,'w') as out2:
        col_names=f'chr\tpos\t{sample_name}'
        print(trans_exon_value,file=out2)
        print(col_names,file=out)
        print(result_value,sep='\t',file=out)
