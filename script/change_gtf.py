import re
import sys

peak_fa=sys.argv[1]
gtf=sys.argv[2] #这里我通常调用的是最长转录本信息的gtf，eg:/project/personal/chenzh/genome/hg38_p13/hg38_p13_longest_transcipt.gtf

with open(peak_fa,"r") as p:
    gene_detail={} #gene_id为键,重复的gene就多个标号，方向和bedlist,peak位置最小值peak位置最大值为值
    gene_lst=[]
    for line in p:
        if re.match(r">",line):
            lst=line.strip().split("::")
            gene_id=lst[1]
            gene_number=lst[2]
            gene_direction=lst[4]
            peak_bed=lst[5].split("||")
            peak_min=peak_bed[0].split(":")[1].split("-")[0]
            peak_max=peak_bed[-1].split(":")[1].split("-")[1]
            chr_id=peak_bed[0].split(":")[0]
            gene_lst.append(gene_id)
            gene_id=f'{gene_id}_{gene_number}'
            gene_detail[gene_id]=[gene_direction,peak_bed,peak_min,peak_max,chr_id]
            

with open(gtf,'r') as g:
    for line in g:
        line=line.strip()
        if re.match(r"^#",line):
            print(line)
        else:
            lst=line.split('\t')
            gene_id=re.search(r"gene_id \"(.*?)\";",lst[-1]).group(1)
            if gene_id not in gene_lst:
                print(line)
    for key,value in gene_detail.items():
        chr_id=value[-1]
        print(chr_id,"HAVANA","gene",value[2],value[3],".",value[0],".",f'gene_id "{key}";',sep="\t")
        print(chr_id,"HAVANA","transcript",value[2],value[3],".",value[0],".",f'gene_id "{key}";',sep="\t")
        for line in value[1]:
            peak_start=line.split(":")[1].split("-")[0]
            peak_end=line.split("-")[1].strip()             
            print(chr_id,"HAVANA","exon",peak_start,peak_end,".",value[0],".",f'gene_id "{key}";',sep="\t")
