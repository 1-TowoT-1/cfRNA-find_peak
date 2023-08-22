rule del_bed:
    input:
        target_gene=config['target_gene'],
        gtf=config['input_gtf']
    output:
        temp(config['project_dir'] + "/raw_exon.bed"),
        temp(config['project_dir'] + "/raw_sort_exon.bed"),
        config['project_dir'] + "/merge_exon.bed"
    params:
        config['project_dir']
    shell:
        """
        cut -f 1 {input.target_gene} |while read i;do grep $i {input.gtf} ;done \
        |awk -F '\\t' '{{if($3=="exon")print $1"\\t"$4-1"\\t"$5}}' > {output[0]} && \
        bedtools sort -i {output[0]} > {output[1]} && \
        bedtools merge -i {output[1]} > {output[2]}
        """
 
rule count_depth:
    input:
        target_gene=config['target_gene'],
        gtf=config['input_gtf'],
        bed=config['project_dir'] + "/merge_exon.bed",
        cancer_bam=config['cancer_bam']
    output:
        config['project_dir'] + "/merge_depth_coverage.txt",
        config['project_dir'] + "/longest_trans.txt",
    threads:
        20
    params:
        scripts=config['script_dir'],
        outdir=config['project_dir'],
        sample_name=config['group'][Treat_name]
    shell:
        """
        python3 {params.scripts}/del.py {input.bed} \
        {input.cancer_bam} \
        {input.target_gene} \
        {params.sample_name} \
        {input.gtf} \
        {params.outdir} \
        {threads} &&\
        cut -f 1 {input.target_gene} |while read i;do grep $i {input.gtf};done |awk -F '\\t' '{{if($3=="transcript")print $9"\\t"$7}}' |awk -F '[\\t;]' '{{print $1"\\t"$4"::"$NF"\\t"$2}}' > {output[1]}
        """