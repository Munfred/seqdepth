from glob import glob
import anndata
import pandas as pd
from tqdm import tqdm
import anndata
import os
import subprocess
import numpy as np
from scipy import sparse
from natsort import natsorted
import io 
import requests

# snakemake -j 1 -s snakemake_h5ads.py --keep-going --rerun-incomplete -pn
# snakemake -j 1 -s snakemake_h5ads.py --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 500"
# snakemake -j 100 -s snakemake_h5ads.py --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 500   --output=./logs/slurm_%j.snakemake_h5ads"
url="https://docs.google.com/spreadsheets/d/"+ \
    "1-2bLIns8r8VRoDenHVk-cQE9feNDnXJXnGZNm70ROrA"+\
    "/export?gid="+\
    "0" + \
    "&format=csv"
metadatas=pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')))

def fetch_subsampling_depths(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    subsampling_depths = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['subsampling_depths'].values[0]
    subsampling_depths = [int(x) for x in subsampling_depths.split(',')]
    return subsampling_depths


final_stacked_h5ad_files = []

for dataset_sample_id in metadatas[metadatas['h5ad']==1]['dataset_sample_id']:
    dataset_project_id = metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['dataset_project_id']
    final_stacked_h5ad_files.append('stacked_h5ads/' + dataset_project_id + '-'+ dataset_sample_id + '-stacked.h5ad')
    
# print( '💙💙💙💙💙💙💙💙💙💙    FINAL H5AD FILES 💙💙💙💙💙💙💙💙💙💙💙💙')
# print(final_stacked_h5ad_files)
print( '🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺   MAKING  ',len(final_stacked_h5ad_files), ' H5AD FILES  🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺🔺')

rule all:
    input:
        final_stacked_h5ad_files

rule stack_h5ad:
    input:
        FULL_H5AD = 'subsampling/{dataset_project_id}/{dataset_sample_id}/kbsub_0/output/counts_unfiltered/adata.h5ad',
        FULL_BUS_FILE = 'subsampling/{dataset_project_id}/{dataset_sample_id}/kbsub_0/output/output.unfiltered.bus',
    params:
        FILTERED_CELL_BARCODES_FILE =  'subsampling/{dataset_project_id}/{dataset_sample_id}/kbsub_0/output/filtered_cell_barcodes.txt',
    output:
        STACKED_H5AD='stacked_h5ads/{dataset_project_id}-{dataset_sample_id}-stacked.h5ad',
        
    run:
#         print('FULL_H5AD:', input.FULL_H5AD)
#         print('STACKED_H5AD:', output.STACKED_H5AD)
        ds = wildcards.dataset_sample_id
        project_id = wildcards.dataset_project_id
        print('Now processing dataset: ', ds)

        bus_whitelist_command = 'bustools whitelist -o ' + params.FILTERED_CELL_BARCODES_FILE + ' ' + input.FULL_BUS_FILE

        returned_value = subprocess.call(bus_whitelist_command , shell=True)  # returns the exit code in unix
        print('Bus whitelist return value (0 is good): ', returned_value)

        filtered_cell_barcodes = pd.read_csv(params.FILTERED_CELL_BARCODES_FILE, index_col=0, header = None)
        print('Number of selected cells:', len(filtered_cell_barcodes))

        base_ad = anndata.read(input.FULL_H5AD)

        # now keep only filtered barcodes
        base_ad = base_ad[base_ad.obs.index.isin(filtered_cell_barcodes.index)].copy()
        print('💚💚💚💚💚💚💚💚💚💚💚💚💚   ', wildcards.dataset_sample_id,'   💚💚💚💚💚💚💚💚💚💚💚💚💚💚💚')
        subdepths = [int(x) for x in metadatas[metadatas['dataset_sample_id']==wildcards.dataset_sample_id]['subsampling_depths'].values[0].split(',')]
        print(subdepths)
        sub_h5ads = ['subsampling/' + wildcards.dataset_project_id + '/' + wildcards.dataset_sample_id + '/' + 'kbsub_' + str(sub) + '/output/counts_unfiltered/adata.h5ad' for sub in subdepths] 
#         print(sub_h5ads)
        print('📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐📐')
#         print(natsorted(sub_h5ads))
        # sort files by increasing depth, but starting from the full file (label 0)
        for adata_file in tqdm(natsorted(sub_h5ads)):
            
            ss_depth = adata_file.split('kbsub_')[-1].split('/')[0]

            print(ss_depth)
            curr_adata = anndata.read(adata_file)

            shared_barcodes = curr_adata.obs_names.intersection(base_ad.obs_names)
            base_idx = base_ad.obs_names.get_indexer(shared_barcodes)
            curr_idx = curr_adata.obs_names.get_indexer(shared_barcodes)

            idxmtx = sparse.eye(base_ad.shape[0], format="csr")[:, base_idx]
            # For some reason this is way faster than equivalent indexing statement
            base_ad.layers[ss_depth] = idxmtx @ curr_adata.X[curr_idx]

        # This is a check to ensure the anndata layers are correct
        for depth in base_ad.layers.keys():
            print('UMIs: \t', int(base_ad.layers[depth].sum()),'\t Reads: \t', depth )
        print(base_ad)
        base_ad.write(output.STACKED_H5AD)