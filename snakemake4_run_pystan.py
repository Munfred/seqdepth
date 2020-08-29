import numpy as np
import pandas as pd
import pystan as ps
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import plotly.express as px
import glob
import arviz
from tqdm import tqdm
import matplotlib
import os
import sys
import datetime
import io
import requests

# --output=%j.out
# snakemake -j 1 -s snakemake4_run_pystan.py --keep-going --rerun-incomplete -pn
# snakemake -j 100 -s snakemake4_run_pystan.py --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 500  -o ./logs/output.%a.pystan "
# snakemake -j 300 -s snakemake4_run_pystan.py --jobname "{jobid}.{wildcards.dataset_project_id}.{wildcards.dataset_sample_id}.pystan" --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 5000  --mem-per-cpu=20200 --output=./logs/snakemake4_run_pystan%j.logs"


url="https://docs.google.com/spreadsheets/d/"+ \
    "1-2bLIns8r8VRoDenHVk-cQE9feNDnXJXnGZNm70ROrA"+\
    "/export?gid="+\
    "0" + \
    "&format=csv"
metadatas=pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')))

pystan_success_indicator_files = []

for dataset_sample_id in metadatas[metadatas['pystan']==1]['dataset_sample_id']:
    dataset_project_id = metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['dataset_project_id']
    pystan_success_indicator_files.append('pystan_results/' + dataset_project_id + '-'+ dataset_sample_id + '-pystan_result_SUCCESS.txt')


rule all:
    input:
        pystan_success_indicator_files,
        
rule run_pystan:
    input: 
        scvi_final_summary_file='scvi_final_summaries/{dataset_project_id}-{dataset_sample_id}-final_summary.csv',
        
    output:
        pystan_success_indicator_file ='pystan_results/{dataset_project_id}-{dataset_sample_id}-pystan_result_SUCCESS.txt',
        
    params:
        results_folder='./pystan_results/',
    
    threads: 8
    run:
        ds = wildcards.dataset_sample_id
        dataset = ds
        project = wildcards.dataset_project_id
        
        print('   💝 💝 💝 💝 💝    PYSTAN PROCESSING DATASET: ', ds, ' PROJECT: ', project, '   💝 💝 💝 💝 💝    ')
        df = pd.read_csv(input.scvi_final_summary_file).sort_values(["sampled_cells", "total_UMIs"], ascending = (True, True))
        
        stan_model = ps.StanModel(file="transformed_piecewise_stan_model.stan", 
                          model_name = "transformed_piecewise_stan_model")

        os.environ['STAN_NUM_THREADS'] = "10"

        print(dataset)
        begin = datetime.datetime.now()
        print ("🧡 🧡 🧡 🧡 🧡 Start fit time : 🧡 🧡 🧡 🧡 🧡 ")
        print (begin.strftime("%Y-%m-%d %H:%M:%S"))

        data_dict = {"ncells": df["sampled_cells"], "umis_per_cell": df["UMIs_per_cell"], "validation_error": df["validation_error"], "N": len(df)}


        stan_fit = stan_model.sampling(data=data_dict,
                               iter=20000,
#                                 warmup = 15000,
                                n_jobs=8,
                                chains=8,
                                thin=5, # I was getting problems with autocorrelation with thin=1...
                                refresh = 200,
                                verbose=True,
                              control={'adapt_delta':1, 'max_treedepth': 20},
                                      )
#         print(stan_model.model_code)
        print ('💚 💚 💚 💚 💚 Finished ', project, '-', ds, '  fit time: 💚 💚 💚 💚 💚')
        now = datetime.datetime.now()
        print (now.strftime("%Y-%m-%d %H:%M:%S"))
        print('Time taken:')
        delta=now - begin
        print(str(delta))
        
        
        s = stan_fit.summary()
        summary = pd.DataFrame(s['summary'], columns=s['summary_colnames'], index=s['summary_rownames'])
        summary_head=pd.concat([summary.head(10),summary.iloc[-10:-1]]).copy()
        print(summary_head)


        arviz.plot_trace(stan_fit,['intercept',
                                   'bp',
                                   'bp_umis',
                                   'before_variance',
                                   'after_variance',
                                   'cell_slope_before_destandardized',
                                   'cell_slope_after_destandardized',
                                   'umi_slope_before_destandardized',
                                   'umi_slope_after_destandardized',
                                   'umi_slope_before',
                                   'umi_slope_after',
                                   'cell_slope_before',
                                   'cell_slope_after',                               
                                  ]
                        )
        plt.savefig(params.results_folder + project + '-' + dataset + str(now.strftime("+%Y-%m-%d_%H:%M:%S")) +'.png',format='png',dpi=80)


        full_stan_results = stan_fit.to_dataframe()

        summary_text = str(summary_head.round(3))

        extracted = stan_fit.extract()
        full_stan_results.to_csv(params.results_folder + project + '-' + dataset + str(now.strftime("+%Y-%m-%d_%H:%M:%S+")) + 'pystan_result_full.csv')
        
        summary.to_csv(params.results_folder + project + '-' + dataset + str(now.strftime("+%Y-%m-%d_%H:%M:%S+")) + 'pystan_result_summary.csv')

        print ("Done! Current date and time : ")
        print (now.strftime("%Y-%m-%d %H:%M:%S"))
        print('   💙 💙 💙 💙 💙    PROCESSED: ', project, '-', ds, '   💙 💙 💙 💙 💙  ')
        
        summary.to_csv(output.pystan_success_indicator_file)


