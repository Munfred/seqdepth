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


# load the 10xv3 results with 30x sampling for each cell/depth combination
dfs={}
for item in glob.glob('./10xv3_final_summaries/*'):
    dfs[item.split('/')[2].split('-final_summary.csv')[0]] = pd.read_csv(item).sort_values(["sampled_cells", "total_UMIs"], ascending = (True, True))

stan_model = ps.StanModel(file="./stan_models/seqdepth_2predictors_piecewise_v5.stan", 
                          model_name = "seqdepth_2s_piecewise")



results_folder='./results6/'
# for dataset in dfs:

# set environmental variable STAN_NUM_THREADS
# Use 4 cores per chain
os.environ['STAN_NUM_THREADS'] = "10"

for dataset in tqdm(['10x_genomics_data-neuron_1k_v3']):
    print(dataset)
    begin = datetime.datetime.now()
    print ("Start fit time : ")
    print (begin.strftime("%Y-%m-%d %H:%M:%S"))
    df = dfs[dataset]
    
    
    data_dict = {"ncells": np.log2(df["sampled_cells"]), "umis_per_cell": np.log2(df["UMIs_per_cell"]), "validation_error": np.log2(df["validation_error"]), "N": len(df)}

    # use the default 4 chains == 4 parallel process
    # used cores = min(cpu_cores, 4*STAN_NUM_THREADS)
    
    
    stan_fit = stan_model.sampling(data=data_dict,
                           iter=25000,
                            warmup = 15000,
#                             n_jobs=10,
                            chains=2,
                            refresh = 10,
                            verbose=True,
                          control={'adapt_delta':0.8, 'max_treedepth': 12},
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
    plt.savefig(results_folder + dataset+'-'+ str(now.strftime("%Y-%m-%d_%H:%M:%S")) +'.png',format='png',dpi=200)


    full_stan_results = stan_fit.to_dataframe()
#         full_stan_results.to_csv(results_folder + full_stan_' + dataset+'-'+str(ncells)+'.csv')
#         summary.to_csv(results_folder + summary_stan_' + dataset+'-'+str(ncells)+'.csv')
#     plt.show()

    summary_text = str(summary_head.round(3))

    extracted = stan_fit.extract()
    full_stan_results.to_csv(results_folder + 'full_stan_2predictors_' + dataset+'.csv')
    now = datetime.datetime.now()
    summary.to_csv(results_folder + 'summary_stan_2predictors_' + dataset + str(now.strftime("%Y-%m-%d_%H:%M:%S")) + '.csv')
    
now = datetime.datetime.now()
print ("Done! Current date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))