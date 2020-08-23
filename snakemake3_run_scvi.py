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
# how to call snakemake just to see if it works
# snakemake -j 1 -s snakemake3_run_scvi.py --keep-going --rerun-incomplete -pn

# how to call snakemake with a maximum of 100 cores used in parallel
# --mem-per-cpu=60200 requests 60GB RAM per CPU/core so that we avoid out of memory problems by having too many jobs in the same node

# snakemake -j 100 -s snakemake3_run_scvi.py --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 500 --mem-per-cpu=60200  --output=./logs/scvi_%j.scvi"


url="https://docs.google.com/spreadsheets/d/"+ \
    "1-2bLIns8r8VRoDenHVk-cQE9feNDnXJXnGZNm70ROrA"+\
    "/export?gid="+\
    "0" + \
    "&format=csv"
metadatas=pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')))

#how many times to run scvi
NUMBER_OF_SCVI_REPLICATES=30

final_summary_files = []

for dataset_sample_id in metadatas[metadatas['scvi']==1]['dataset_sample_id']:
    dataset_project_id = metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['dataset_project_id']
    final_summary_files.append('final_summaries/' + dataset_project_id + '-'+ dataset_sample_id + '-final_summary.csv')
    
# print( '💙💙💙💙💙💙💙💙💙💙    FINAL H5AD FILES 💙💙💙💙💙💙💙💙💙💙💙💙')
# print(final_stacked_h5ad_files)
# print( '🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺   MAKING  ',len(final_summary_files), ' H5AD FILES  🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺')

def make_partial_results_filenames(wildcards):
#     print('💙💙💙💙')
    STACKED_H5AD=f'stacked_h5ads/{wildcards.dataset_project_id}-{wildcards.dataset_sample_id}-stacked.h5ad'
    adata=anndata.read(STACKED_H5AD)
    total_cells = adata.n_obs
    n_retained_cells = int(0.85*total_cells)
    print('Total cells:', total_cells)
    print('Retained cells:', n_retained_cells)

    cells_sizes = []
    #initial number of sampled cells
    sampling_size = 500
    print('💙💙💙💙')
    print('total cells:', total_cells)
    print(sampling_size)
    print(n_retained_cells)
    print( '🔺🔺🔺')
    while sampling_size < n_retained_cells:
        cells_sizes.append(sampling_size)
        sampling_size = int(sampling_size*np.sqrt(2))
        
    ss_depths = [depth for depth in adata.layers.keys()]
    print('🥖🥖🥖🥖🥖🥖🥖')
    print(dataset_project_id)
    print(dataset_sample_id)
    print(cells_sizes)
    print(ss_depths)
    print('🥖🥖🥖🥖🥖🥖🥖')
        
    return expand(
    'scvi_output/partial_csvs/{{dataset_project_id}}-{{dataset_sample_id}}-c{ss_cells}-d{ss_depth}-SUCCESS.csv',
    ss_cells=cells_sizes, ss_depth=ss_depths    
    )

rule all:
    input:
        final_summary_files
        
rule run_scvi:
    input: 
        STACKED_H5AD='stacked_h5ads/{dataset_project_id}-{dataset_sample_id}-stacked.h5ad',
#         STACKED_H5AD='stacked_h5ads/'+ dataset_project_id + '-' + dataset_sample_id + '-stacked.h5ad',
        
    output:
        SCVI_PARTIAL_SUMMARY='scvi_output/partial_csvs/{dataset_project_id}-{dataset_sample_id}-c{ss_cells}-d{ss_depth}.csv',
        SCVI_PARTIAL_SUMMARY_SUCCESS='scvi_output/partial_csvs/{dataset_project_id}-{dataset_sample_id}-c{ss_cells}-d{ss_depth}-SUCCESS.csv',
    
    run:
        
        # try loading partial results of a previous run that might have failed, otherwise make a new dataframe
        try:
            df_full=pd.read_csv(output.SCVI_PARTIAL_SUMMARY, index=False)
        except:
            df_full=pd.DataFrame()
        
        current_numer_of_replicates=len(df_full)
        
        while current_numer_of_replicates < NUMBER_OF_SCVI_REPLICATES:
            ds = wildcards.dataset_sample_id
            ss_depth = int(wildcards.ss_depth)
            ss_cells = int(wildcards.ss_cells)

            print('============== PROCESSING DATASET: ', ds, 'subsampling depth: ', ss_depth, 'cells: ', ss_cells, '===============')
            adata = anndata.read(str(input.STACKED_H5AD))
            total_cells = adata.n_obs
            n_retained_cells = int(0.85*total_cells)
            print('Total cells:', total_cells)
            print('Retained cells:', n_retained_cells)

            cells_sizes = []
            #initial number of sampled cells
            sampling_size = 500
            while sampling_size < n_retained_cells:
                cells_sizes.append(sampling_size)
                sampling_size = int(sampling_size*np.sqrt(2))

            # cells_sizes = np.logspace(np.log2(500), np.log2(n_retained_cells), num=9, base=2).astype(int)
            print('Number of sampled cells for ', ds, cells_sizes)

            cells_dataset = GeneExpressionDataset()
            X_ = adata.layers['0']
            cells_dataset.populate_from_data(X_, gene_names=adata.var.index.values)

            #we subsambple to 1000 genes for speed and to prevent overfitting
            cells_dataset.filter_genes_by_count(per_batch=True)
            cells_dataset.subsample_genes(1000)
            sel_genes = cells_dataset.gene_names

            n_validation = adata.shape[0] - n_retained_cells
            print(ds, ' n_validation:', n_validation)

            validation_cells = np.random.choice(adata.obs.index, size=n_validation, replace=False)
            learning_cells = adata.obs.index.difference(validation_cells)

            val_adata = adata[validation_cells]
            lea_adata = adata[learning_cells]

            ne_cells = X_.sum(axis=1) > 0
            to_keep = np.where(ne_cells)[0]

            log_counts = np.log(X_[to_keep].sum(axis=1))
            local_means = (np.mean(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)
            local_vars = (np.var(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)

            batch_indices = np.ones((X_.shape[0], 1))
            labels = np.zeros_like(batch_indices)

            validation_cells_dataset = GeneExpressionDataset()
            X_ = val_adata[:, sel_genes].layers['0']
            validation_cells_dataset.populate_from_data(X_)

            # this will store all dataframes we also save as csv
            results_list = []

            print('                🍞  🥯  🧁  🍰  🥮  🥞  🥧  🎂  🍯  🥪  🥖  🥐  🥟   ')
            print('======== ', ds, ' ========  SUBSAMPLED CELLS:', ss_cells, ' ================== ')
            idx = np.random.choice(learning_cells, size=ss_cells, replace=False)
            tmp_adata = lea_adata[idx]
            n_epochs = int(27 * adata.shape[0] / ss_cells)
            print('Running ', n_epochs, 'epochs')

            print(ds, 'now running depth:', ss_depth)

            result_dict = {'ss_depth': ss_depth, 'ss_cells': ss_cells}

            print('🧁  🍰 🧁  🍰 🧁  🍰 🧁  🍰')
            print(str(input.STACKED_H5AD))
            print(tmp_adata[:, sel_genes].layers)
            print(tmp_adata[:, sel_genes].layers.keys)
            try:
                X_ = tmp_adata[:, sel_genes].layers[str(ss_depth)]
            except:
                pass
            try:
                X_ = tmp_adata[:, sel_genes].layers[ss_depth]
            except:
                pass
            try:
                X_ = tmp_adata[:, sel_genes].layers[wildcards.ss_depth]
            except:
                X_ = tmp_adata[:, sel_genes].layers[wildcards.ss_depth]
            ne_cells = X_.sum(axis=1) > 0
            to_keep = np.where(ne_cells)[0]

            log_counts = np.log(X_[to_keep].sum(axis=1))
            local_means = (np.mean(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)
            local_vars = (np.var(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)

            batch_indices = np.ones((X_.shape[0], 1))
            labels = np.zeros_like(batch_indices)

            cells_dataset = GeneExpressionDataset()
            cells_dataset.populate_from_data(X_, gene_names=sel_genes)

            vae = VAE(cells_dataset.nb_genes,
                      dispersion='gene',
                      n_layers=2,
                      n_hidden=128,
                      reconstruction_loss='nb')

            trainer = UnsupervisedTrainer(vae, cells_dataset)
            
            # very rarely the trainer fails with nans for the loss, 
            # simplest solution is to just try a few times before giving up
            # very unlikely it will fail 7 times in a row, if it does there's a real problem...
            try:
                trainer.train(n_epochs=n_epochs)
            except:
                try:
                    trainer.train(n_epochs=n_epochs)
                except:
                    try:
                        trainer.train(n_epochs=n_epochs)
                    except:
                        try:
                            trainer.train(n_epochs=n_epochs)          
                        except:
                            try:
                                trainer.train(n_epochs=n_epochs)
                            except:
                                try:
                                    trainer.train(n_epochs=n_epochs)
                                except:
                                    trainer.train(n_epochs=n_epochs)
                                    
            validation_posterior = Posterior(vae, validation_cells_dataset, use_cuda=False)
            print(X_.shape)
            result_dict['validation_error'] = validation_posterior.reconstruction_error()

            df=pd.DataFrame()
            df['ss_depth'] = [result_dict['ss_depth']]
            df['ss_cells'] = [result_dict['ss_cells']]
            df['validation_error'] = [result_dict['validation_error']]
            
            df_full=df_full.append(df)
            df_full.to_csv(output.SCVI_PARTIAL_SUMMARY, index=False)
            current_numer_of_replicates=len(df_full)
            
        # once the dataframe has the desired number of replicates we exit the while loop and save the final version
        df_full.to_csv(output.SCVI_PARTIAL_SUMMARY_SUCCESS, index=False)

rule make_final_summaries:
    input:
        SCVI_PARTIAL_SUMMARY_SUCCESS=make_partial_results_filenames,
        
    params:
        STACKED_H5AD='stacked_h5ads/{dataset_project_id}-{dataset_sample_id}-stacked.h5ad',
    output:
        FINAL_SUMMARY = 'final_summaries/{dataset_project_id}-{dataset_sample_id}-final_summary.csv'
    run:
        print('🥧   🥧   🥧   🥧 FINISHED DATASET: 🥧   🥧   🥧  🥧')
        print(input)
        print('🍰  🍰  🍰  🍰  🍰  🍰  🍰   🍰   🍰   🍰')
        dfs = []
        ds = wildcards.dataset_sample_id
        summary_results={}

        for fname in tqdm(input):
            df = pd.read_csv(fname)
            dfs.append(df)

        summary_results[ds] = pd.concat(dfs).reset_index(drop=True)

        #loads anndata object used for calculating UMIs seen
        
        adata = anndata.read(params.STACKED_H5AD)
        total_cells = adata.n_obs

        # this uses the anndata object to calculate the total UMIs for each subsampling depth
        for subsampled_depth in summary_results[ds]['ss_depth'].unique():
            total_UMIs = int(adata.layers[str(subsampled_depth)].sum())
            summary_results[ds].loc[summary_results[ds].ss_depth == subsampled_depth, 'total_UMIs'] = total_UMIs
        summary_results[ds]['total_UMIs'] = summary_results[ds]['total_UMIs'].astype('int')
        # rfull data depth from `0` to the total reads in the dataset,
        summary_results[ds]['ss_depth'] = summary_results[ds]['ss_depth'].map(lambda d: {0: int(metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['total_reads'].values)}.get(d, d))
        # we rename ss_depth to total_seqtk_reads to avoid abiguity
        summary_results[ds].rename(columns={'ss_depth':'total_seqtk_reads'}, inplace=True)
        summary_results[ds].rename(columns={'ss_cells':'sampled_cells'}, inplace=True)
        summary_results[ds]['validation_error'] = summary_results[ds]['validation_error'].round(1)

        summary_results[ds]['estimated_reads'] = (summary_results[ds]['sampled_cells'] / total_cells* summary_results[ds]['total_seqtk_reads']).astype('int')

        summary_results[ds]['estimated_UMIs'] = (summary_results[ds]['sampled_cells'] / total_cells * summary_results[ds]['total_UMIs']).astype('int')
        summary_results[ds]['reads_per_cell'] = (summary_results[ds]['total_seqtk_reads'] / total_cells).astype('int')

        summary_results[ds]['UMIs_per_cell'] = (summary_results[ds]['total_UMIs'] / total_cells).astype('int')
        summary_results[ds]['duplication_rate'] = (summary_results[ds]['total_seqtk_reads'] / summary_results[ds]['total_UMIs'] ).round(2)


#         summary_results[ds].drop(['index'], axis=1, inplace=True)
        summary_results[ds] = summary_results[ds].sort_values('total_seqtk_reads')
        summary_results[ds].index.rename(ds, inplace=True)
        summary_results[ds]['dataset'] = ds
        print(ds, 'summary results shape:' , summary_results[ds].shape)
        print('Number of points with missing validation errors: ', sum(summary_results[ds]['validation_error'].isna()) )


        summary_results[ds].to_csv(output.FINAL_SUMMARY, index=False)

