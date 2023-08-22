rule change_header:
    input:
        config['project_dir'] + '/03.EXP/All_sample_count.txt'
    output:
        temp(touch(config['project_dir'] + '/change_header.task'))
    shell:
        """
        sed -i '1d' {input} &&
        sed -i '1i\Gene_id\t{Treat_sample}\t{Control_sample}' {input}
        """


rule DEG:
    input:
        config['project_dir'] + '/change_header.task',
        config['project_dir'] + '/03.EXP/All_sample_count.txt'
    output:
        config['project_dir'] + "/04.DEG/" + config['DEG_info'] + "_raw_" + config['DEG_software'] + ".xls"
    params:
        scripts=config['script_dir'],
        DEG_software=config['DEG_software'],
        groups_samples_info=';'.join([f'{key}:{value}' for key, value in config['group'].items()]),
        outdir=config['project_dir'] + '/04.DEG',
        DEG_info=config['DEG_info']
    shell:
        """
        Rscript {params.scripts}/DEG_analysis.R {params.DEG_software} \
        "{params.groups_samples_info}" \
        {params.DEG_info} \
        {input[1]} \
        {params.outdir}
        """


rule DEG_summary:
    input:
        config['project_dir'] + "/04.DEG/" + config['DEG_info'] + "_raw_" + config['DEG_software'] + ".xls",
        config['project_dir'] + "/01.Peak/" + "peak_" + config['min_peak_h'] + "_" + config['sample_min'] + ".fa"
    output: 
        temp(config['project_dir'] + "/04.DEG/temp.txt")
    run:
        with open(input[1],'r') as f1,open(input[0],'r') as f2,open(output[0],'w') as o1:
            gene_id_lst=[]
            for line in f1:
                if re.match(r">",line):
                    lst=line.split('::')
                    gene_id_lst.append(lst[1])
            gene_id_set=set(gene_id_lst)       
            header=f2.readline().strip()
            print(header,file=o1,sep="\t")
            for line in f2:
                line=line.strip()
                lst=line.split('\t')
                for i in gene_id_set:
                    if re.match(i,lst[0]):
                        print(line,file=o1,sep='\t')



rule Gini_AUC:
    input:
        config['project_dir'] + "/04.DEG/temp.txt"
    output:
        temp(config['project_dir'] + "/04.DEG/temp2.txt")
    params:
        s_dir=config['script_dir'],
        T_sample=','.join(Treat_sample.split('\t')),
        C_sample=','.join(Control_sample.split('\t'))
    shell:
        """
        Rscript {params.s_dir}/Gini_AUC.R {input} {params.T_sample} {params.C_sample} {output}
        """



rule ratio_count:
    input:
        config['project_dir'] + '/01.Peak/peak_' + config['min_peak_h'] + '_' + config['sample_min'] + '.fa'
    output:
        temp(config['project_dir'] + '/05.Peak_rate/temp.txt')
    params:
        all_sample=config['all_sample_num'],
        input_dir=config['project_dir'],
        script=config['script_dir'],
        mp=config['min_peak_h'],
        n=config['sample_min']
    shell:
        """
        python3 {params.script}/ratio_count.py {input} {params.input_dir} {params.all_sample} {output} {params.mp} {params.n}
        """


rule rate_summary:
    input:
        config['project_dir'] + "/04.DEG/temp2.txt",
        config['project_dir'] + '/05.Peak_rate/temp.txt'
    output:
        config['project_dir'] + '/05.Peak_rate/result_' + config['min_peak_h'] + '_' + config['sample_min'] + '.txt'
    shell:
        """
        awk -F '\\t' 'NR==FNR{{a[$1]=$2"\\t"$3;next}}{{print $0"\\t"a[$1]}}' {input[1]} {input[0]} > {output}
        """