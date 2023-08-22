import sys
import os
import pandas as pd

# eg：python3 feature_summary.py All_sample_featurecount.txt All_sample .

inpunt=sys.argv[1]
output_prex=sys.argv[2] #文件名前缀，eg：All_sample
outdir=sys.argv[3]

df=pd.read_csv(inpunt,sep='\t',header=1)

sample_list=[]
new_df=pd.DataFrame()
new_df['Gene_id']=df['Geneid']
for i in df.columns[6:]:
    sample_name=os.path.basename(i).split(".fix.sort.rdup.bam")[0]
    sample_list.append(sample_name)
    
    sum_number=df[i].sum()
    new_df[sample_name]=round(df[i].apply(lambda x: x*1000000000/sum_number)/df["Length"],2)
    df=df.rename(columns={i:sample_name})
df.rename(columns={"Geneid":"Gene_id"})
df2=df.iloc[:, [0] + list(range(6, len(df.columns)))]

out_count=f'{outdir}/{output_prex}_count.txt'
out_fpkm=f'{outdir}/{output_prex}_FPKM.txt'

df2.to_csv(out_count,sep='\t',index=False)
new_df.to_csv(out_fpkm,sep='\t',index=False)
