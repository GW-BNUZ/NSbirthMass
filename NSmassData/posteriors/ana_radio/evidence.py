import bilby
import numpy as np
import csv
import pandas as pd
import os
import glob

sub_dir_name_list=['turn_on_pow','turn_on_pow_fix','turn_on_pow_G','pow', '2G','turn_on_pow_G_fixed_max','G','2G_fixed_max', 'sst', 'G_fixed_max', '2G_fixed_min','2G_fixed_max_min','logu','3G_fixed_max_min','lognorm','gamma', 'U']

dir_path=os.path.abspath( os.path.join(os.getcwd()) )
dir_name=os.path.basename(dir_path)

evidence={}
for i in range(len(sub_dir_name_list)):
    try:
        os.path.exists(dir_path+'/'+sub_dir_name_list[i]+"/hy_outdir/*.json")
        fnames= glob.glob(dir_path+'/'+sub_dir_name_list[i]+"/hy_outdir/*.json")
        fname=fnames[0]
        re=bilby.result.read_in_result(filename=fname)
        evidence[sub_dir_name_list[i]]=re.log_evidence
    except:
        evidence[sub_dir_name_list[i]]=0
print(evidence)
evidence=pd.DataFrame.from_dict(evidence,orient='index')

evidence.to_csv('{}_evidence.csv'.format(dir_name),header=0,index=sub_dir_name_list)
