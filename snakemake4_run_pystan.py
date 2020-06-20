import numpy as np
import sklearn
import pandas as pd
from glob import glob
from tqdm import tqdm
import anndata
import torch
from scvi.dataset import GeneExpressionDataset
from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
from scvi.inference.posterior import Posterior
from openTSNE import TSNE
from openTSNE.callbacks import ErrorLogger
import os
import sys
from natsort import natsorted
import io 
import requests

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

# --output=%j.out
# snakemake -j 1 -s snakemake4_run_pystan.py --keep-going --rerun-incomplete -pn
# snakemake -j 100 -s snakemake4_run_pystan.py --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 500  -o ./logs/output.%a.out "
# snakemake -j 300 -s snakemake4_run_pystan.py --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 500   --output=./logs/snakemake4_run_pystan%j.logs"


url="https://docs.google.com/spreadsheets/d/"+ \
    "1-2bLIns8r8VRoDenHVk-cQE9feNDnXJXnGZNm70ROrA"+\
    "/export?gid="+\
    "0" + \
    "&format=csv"
metadatas=pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')))


scvi_final_summary_files = []
pystan_final_summary_files =[]
pystan_final_full_files = []
pystan_success_indicator_files = []

for dataset_sample_id in metadatas[metadatas['pystan']==1]['dataset_sample_id']:
    dataset_project_id = metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['dataset_project_id']
    scvi_final_summary_files.append('scvi_final_summaries/' + dataset_project_id + '-'+ dataset_sample_id + '-final_summary.csv')
    pystan_final_summary_files.append('pystan_results/' + dataset_project_id + '-'+ dataset_sample_id + '-pystan_result_summary.csv')
    pystan_final_full_files.append('pystan_results/' + dataset_project_id + '-'+ dataset_sample_id + '-pystan_result_full.csv')
    pystan_success_indicator_files.append('pystan_results/' + dataset_project_id + '-'+ dataset_sample_id + '-pystan_result_SUCCESS.txt')
    
# print( '💙💙💙💙💙💙💙💙💙💙    FINAL H5AD FILES 💙💙💙💙💙💙💙💙💙💙💙💙')
# print(final_stacked_h5ad_files)
# print( '🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺   MAKING  ',len(final_summary_files), ' H5AD FILES  🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺')



rule all:
    input:
        pystan_success_indicator_files
        
rule run_pystan:
    input: 
        scvi_final_summary_file='scvi_final_summaries/{dataset_project_id}-{dataset_sample_id}-final_summary.csv',
        
    output:
        pystan_success_indicator_file ='pystan_results/{dataset_project_id}-{dataset_sample_id}-pystan_result_SUCCESS.txt',
        
    params:
        results_folder='./pystan_results/',
    
    run:
        ds = wildcards.dataset_sample_id
        dataset = ds
        project = wildcards.dataset_project_id
        
        print('   💙 💙 💙 💙    PYSTAN PROCESSING DATASET: ', ds, ' PROJECT: ', project, '   💙 💙 💙 💙   ')
        df = pd.read_csv(input.scvi_final_summary_file).sort_values(["sampled_cells", "total_UMIs"], ascending = (True, True))


        
        stan_model = ps.StanModel(file="piecewise_stan_model.stan", 
                          model_name = "piecewise_stan_model")

        # for dataset in dfs:

        # set environmental variable STAN_NUM_THREADS
        # Use 4 cores per chain
        os.environ['STAN_NUM_THREADS'] = "10"

        print(dataset)
        begin = datetime.datetime.now()
        print ("Start fit time : ")
        print (begin.strftime("%Y-%m-%d %H:%M:%S"))

        data_dict = {"ncells": np.log2(df["sampled_cells"]), "umis_per_cell": np.log2(df["UMIs_per_cell"]), "validation_error": np.log2(df["validation_error"]), "N": len(df)}

        # use the default 4 chains == 4 parallel process
        # used cores = min(cpu_cores, 4*STAN_NUM_THREADS)

        stan_fit = stan_model.sampling(data=data_dict,
                               iter=25000,
                                warmup = 15000,
    #                             n_jobs=10,
                                chains=3,
                                refresh = 10,
                                verbose=True,
                              control={'adapt_delta':0.8, 'max_treedepth': 15},
                                      )
        print(stan_model.model_code)
        print ("Finished fit time: ")
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
                                   'bp_1d',
                                   'bp_umis',
                                   'before_variance',
                                   'after_variance',
                                   'umi_slope_before',
                                   'umi_slope_after',
                                   'cell_slope_before',
                                   'cell_slope_after',                               
                                    'cell_slope_difference',
                                    'cell_after_over_before',
                                    'cell_before_over_after',
                                    'umi_slope_difference', 
                                    'umi_after_over_before', 
                                    'umi_before_over_after', 
                                  ]
                        )
        plt.savefig(params.results_folder + dataset+'-'+ str(now.strftime("%Y-%m-%d_%H:%M:%S")) +'.png',format='png',dpi=80)


        full_stan_results = stan_fit.to_dataframe()

        summary_text = str(summary_head.round(3))

        extracted = stan_fit.extract()
        full_stan_results.to_csv((params.results_folder + project + '-' + dataset + str(now.strftime("_%Y-%m-%d_%H:%M:%S")) + 'pystan_result_full.csv')
        
        summary.to_csv(params.results_folder + project + '-' + dataset + str(now.strftime("_%Y-%m-%d_%H:%M:%S")) + 'pystan_result_summary.csv')

        print ("Done! Current date and time : ")
        print (now.strftime("%Y-%m-%d %H:%M:%S"))
        
        summary.to_csv(output.pystan_success_indicator_file)

