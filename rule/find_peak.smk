rule find_peak:
    input:
        config['project_dir'] + "/longest_trans.txt"
    output:
        config['project_dir'] + '/fp_complete.task'
    params:
        min_peak_h=config['min_peak_h'],
        sample_min=config['sample_min'],
        scripts=config['script_dir'],
        outdir=config['project_dir']
    shell:
        """
        python3 {params.scripts}/find_peak.py -t {input} \
        -mp {params.min_peak_h} \
        -n {params.sample_min} \
        -i {params.outdir} > {output}
        """