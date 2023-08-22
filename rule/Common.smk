Treat_name=config['DEG_info'].split('_vs_')[0]
Control_name=config['DEG_info'].split('_vs_')[1]
Treat_sample='\t'.join(config['group'][Treat_name].split(','))
Control_sample='\t'.join(config['group'][Control_name].split(','))


if config['project_dir'][-1]=="/":
    config['project_dir']=config['project_dir'][:-1]


config['sample_min']=str(config['sample_min'])
config['min_peak_h']=str(config['min_peak_h'])


def get_final_output():
    final_output=[]
    if config['count_depth']:
        final_output.append(config['project_dir'] + "/merge_exon.bed")
        final_output.append(config['project_dir'] + "/merge_depth_coverage.txt")
        final_output.append(config['project_dir'] + "/longest_trans.txt")
    if config['find_peak']:
        final_output.append(config['project_dir'] + '/fp_complete.task')
    if config['exp_count']:
        final_output.append(config['project_dir'] + '/01.Peak/peak_' + config['min_peak_h'] + '_' + config['sample_min'] + '.fa')
        final_output.append(config['project_dir'] + '/02.GTF/new.gtf')
        final_output.append(config['project_dir'] + '/03.EXP/All_sample_featurecount.txt')
        final_output.append(config['project_dir'] + '/03.EXP/All_sample_count.txt')
    if config['peak_rate']:
        final_output.append(config['project_dir'] + "/04.DEG/" + config['DEG_info'] + "_raw_" + config['DEG_software'] + ".xls")
        final_output.append(config['project_dir'] + '/05.Peak_rate/result_' + config['min_peak_h'] + '_' + config['sample_min'] + '.txt')
    return final_output