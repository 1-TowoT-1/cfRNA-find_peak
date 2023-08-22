rule cat_peak:
    input:
        config['project_dir'] + '/fp_complete.task',
        config['project_dir'] + '/longest_trans.txt'
    output:
        config['project_dir'] + '/01.Peak/peak_' + config['min_peak_h'] + '_' + config['sample_min'] + '.fa'
    params:
        mp=config['min_peak_h'],
        n=config['sample_min'],
        i=config['project_dir'],
        script=config['script_dir']
    run:
        with open(input[1],'r') as f1:
            dir_str=[]
            for line in f1:
                gene_id=re.search(r"gene_id \"(.*?)\"",line).group(1)
                gene_name=re.search(r"gene_name \"(.*?)\"",line).group(1)
                trans_id=re.search(r"transcript_id \"(.*?)\"",line).group(1)
                file_name=f'{params.i}/{gene_id}_{gene_name}/{trans_id}/04.peak_{params.mp}_{params.n}.fa'         
                if os.system("ls %s"%file_name)==0:
                    dir_str.append(f'{params.i}/{gene_id}_{gene_name}/{trans_id}/04.peak_{params.mp}_{params.n}.fa')
            dir_str=' '.join(dir_str)
            os.system("cat %s > %s"%(dir_str,output[0]))


rule change_gtf:
    input:
        config['project_dir'] + '/01.Peak/peak_' + config['min_peak_h'] + '_' + config['sample_min'] + '.fa'
    output:
        config['project_dir'] + '/02.GTF/new.gtf'
    params:
        script=config['script_dir'],
        gtf=config['input_gtf']
    shell:
        """
        python3 {params.script}/change_gtf.py {input} {params.gtf} > {output}
        """


rule feature_count:
    input:
        config['project_dir'] + '/02.GTF/new.gtf'
    output:
        config['project_dir'] + '/03.EXP/All_sample_featurecount.txt'
    params:
        cancer_bam=config['cancer_bam'],
        health_bam=config['health_bam']
    threads:
        4
    run:
        with open(params.cancer_bam,'r') as f1,open(params.health_bam,'r') as f2:
            cancer_file=''
            health_file=''
            for line in f1:
                line=line.strip()
                cancer_file=f'{cancer_file} {line}'
            for line in f2:
                line=line.strip()
                health_file=f'{health_file} {line}'
            os.system('featureCounts -p -T %d -t exon -a %s -g gene_id -O -o %s %s %s'%(threads,input[0],output[0],cancer_file,health_file))


rule count_summary:
    input:
        config['project_dir'] + '/03.EXP/All_sample_featurecount.txt'
    output:
        config['project_dir'] + '/03.EXP/All_sample_count.txt'
    params:
        script=config['script_dir'],
        outdir=config['project_dir'] + '/03.EXP'
    shell:
        """
        python {params.script}/feature_summary.py {input} All_sample {params.outdir}
        """