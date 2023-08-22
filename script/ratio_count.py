import re
import sys

target_fa=sys.argv[1]
input=sys.argv[2] #config中的project_dir目录
all_sample=sys.argv[3]
outfile=sys.argv[4]
min_peak_h=sys.argv[5]
sample_min=sys.argv[6]


with open(target_fa,'r') as f1:
    region_dic={}
    for line in f1:
        if re.match(r">",line):
            lst=line.strip().split("::")
            gene_name=re.match(r">(.*?)::",line).group(1)
            gene_id=lst[1]
            number=lst[2]
            trans_id=lst[3]
            len=int(lst[-1].split('=')[1])+1 ##实际上真实read比区间少一位，因为当时用的是bed文件取得序列，然而bed文件是0开头得，这里就当多加一个位点
            indir=f'{input}/{gene_id}_{gene_name}/{trans_id}/'
            region=f'{gene_id}_{number}'
            region_dic[region]=[indir,lst[-2],len,gene_name]


with open(outfile,'w') as o1:
    print("Gene_ID","Gene_name","ratio",sep='\t',file=o1)
    for key,value in region_dic.items():
        pos_lst=[]
        count=0
        lst=value[1].split("||")
        for i in lst:
            start=int(i.split(":")[1].split("-")[0])
            end=int(i.split("-")[1])+1
            pos_lst.append([start,end])

        file_name=f'{region_dic[key][0]}/05.peak_{min_peak_h}_{sample_min}_detail.txt'
        with open(file_name,'r') as f2:            
            for line in f2:
                lst2=line.strip().split('\t')
                pos=int(lst2[1])
                for i in pos_lst:
                    if pos in range(i[0],i[1]):
                        count+=int(lst2[2])
        
        ratio=str(round(count/value[2]/int(all_sample)*100,2)) + '% (all_num=' + all_sample + ')'
        print(key,value[3],ratio,sep='\t',file=o1)
