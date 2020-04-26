import glob
import os
import pandas as pd
import io
import requests
import ntpath
import numpy as np


# for file in ./*/fastqs/full/*S1_L001_R1_001.fastq.gz ; do echo $file; zcat $file | echo $((`wc -l`/4)); 
# snakemake -j 1 -s kbsub.py --keep-going --rerun-incomplete -pn
# snakemake -j 100 -s kbsub.py --keep-going --rerun-incomplete --latency-wait 50 --cluster "sbatch -A lpachter -t 500   --output=./logs/slurm13_%j.kb" --verbose
# squeue -u edaveiga | grep snakejob | awk '{print $1}' | xargs -n 1 scancel
# kb count -i ~/data/references/mus_musculus-ensembl-96/transcriptome.idx -g ~/data/references/mus_musculus-ensembl-96/transcripts_to_genes.txt -x 10xv3 -o out2 ~/data/10x_genomics_data/heart_1k_v3/heart_1k_v3_R1_concat_1.fastq.gz ~/data/10x_genomics_data/heart_1k_v3/heart_1k_v3_R2_concat_2.fastq.gz 

THREAD = 1
REF_PATH = '~/data/references'
# how to call kallisto
KALLISTO = 'kallisto'
#how to call bustools
BUSTOOLS = 'bustools'

# directory with fast r/w if there is one, to speed up kallisto
SCRATCH_DIR = '/central/scratchio/edaveiga/'

# SCRATCH_DIR = 'scratch/'

url="https://docs.google.com/spreadsheets/d/"+ \
    "1-2bLIns8r8VRoDenHVk-cQE9feNDnXJXnGZNm70ROrA"+\
    "/export?gid="+\
    "0" + \
    "&format=csv"
metadatas=pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')))

def make_t2g_path(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    species = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['species'].values[0]
    if species=='mouse':
        T2G = os.path.join(REF_PATH,'mus_musculus-ensembl-96/transcripts_to_genes.txt')
    if species=='human':
        T2G = os.path.join(REF_PATH,'homo_sapiens-ensembl-96/transcripts_to_genes.txt')
    if species=='human-mouse':
        T2G = os.path.join(REF_PATH,'mouse_human_mix_ensembl-96/transcripts_to_genes.txt')
    return T2G

def make_index_path(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    species = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['species'].values[0]
    if species=='mouse':
        INDEX = os.path.join(REF_PATH,'mus_musculus-ensembl-96/transcriptome.idx')
    if species=='human':
        INDEX = os.path.join(REF_PATH,'homo_sapiens-ensembl-96/transcriptome.idx')
    if species=='human-mouse':
        INDEX = os.path.join(REF_PATH,'mouse_human_mix_ensembl-96/transcriptome.idx')
    return INDEX

def make_whitelist_path(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    TECHNOLOGY = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['technology'].values[0]
    WHITELIST = 'nothing'
    if TECHNOLOGY=='10xv3':
        WHITELIST = os.path.join(REF_PATH,'10xv3_whitelist.txt')
    if TECHNOLOGY=='10xv2':
        WHITELIST = os.path.join(REF_PATH,'10xv2_whitelist.txt')
    return WHITELIST

def fetch_technology(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    TECHNOLOGY = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['technology'].values[0]
    return TECHNOLOGY

def fetch_read1_filepath(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    DATASET_SAMPLE_PATH = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['dataset_sample_path'].values[0]   
    READ1_FILES = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['concat_read1_file'].values[0].split(',')
    #remove trailing spaces
    READ1_FILES = [ read1_file.strip() for read1_file in READ1_FILES]
    
    READ1_FILEPATHS = [os.path.join(DATASET_SAMPLE_PATH, read_filename) for read_filename in READ1_FILES]    
    return READ1_FILEPATHS
    
def fetch_read2_filepath(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    DATASET_SAMPLE_PATH = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['dataset_sample_path'].values[0] 
    READ2_FILES = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['concat_read2_file'].values[0].split(',')
    #remove trailing spaces
    READ2_FILES = [ read2_file.strip() for read2_file in READ2_FILES]
    READ2_FILEPATHS = [os.path.join(DATASET_SAMPLE_PATH, read_filename) for read_filename in READ2_FILES]      
    return READ2_FILEPATHS

def fetch_subsampling_depths(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    subsampling_depths = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['subsampling_depths'].values[0]
    subsampling_depths = [int(x) for x in subsampling_depths.split(',')]
    return subsampling_depths

final_subsampled_fastqs = []
final_count_matrices = []
final_genecount_mtx = []
final_bus_txt =[]
final_files = []

for dataset_sample_id in metadatas[metadatas['process']==1]['dataset_sample_id']:
    dataset_project_id = metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['dataset_project_id']
    for sub_string in metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['subsampling_depths']:
        subsampling_depths = [int(x) for x in sub_string.split(',')]
        for subdepth in subsampling_depths:
            final_files.append(str('subsampling/' + dataset_project_id.values[0] + '/' + dataset_sample_id + '/kbsub_'+str(subdepth)+'/output/counts_unfiltered/adata.h5ad'))
print( '==========================================================')
print(final_files)
print('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
# print(str(dataset_project_id + '/' + dataset_sample_id + '/subsampled_'+str(sub)+'/counts_unfiltered/adata.h5ad'))
# print(final_h5ad[0])

rule all:
    input:
        final_files
        
rule subsample:
    input:
        READ1_FILEPATH = fetch_read1_filepath,
        READ2_FILEPATH = fetch_read2_filepath
    params:
        SUBSAMPLING=lambda wildcards: f'{wildcards.sub}',
        SUB_DIR_TMP = os.path.join(SCRATCH_DIR, '{dataset_project_id}/{dataset_sample_id}/tmp_{sub}/'),

    output:
        R1 = os.path.join(SCRATCH_DIR, '{dataset_project_id}/{dataset_sample_id}/tmp_{sub}/subsampledreads_{sub}_R1.fq'),
        R2 = os.path.join(SCRATCH_DIR, '{dataset_project_id}/{dataset_sample_id}/tmp_{sub}/subsampledreads_{sub}_R2.fq')
    shell:
        """       
        if (({params.SUBSAMPLING}==0)); then 
        cp {input.READ1_FILEPATH} {output.R1}
        cp {input.READ2_FILEPATH} {output.R2}
        fi
        
        if (({params.SUBSAMPLING}!=0)); then
        mkdir -p {params.SUB_DIR_TMP}
        cd {params.SUB_DIR_TMP}
        seqtk sample -2 -s100 {input.READ1_FILEPATH} {params.SUBSAMPLING} > {output.R1} && \
        seqtk sample -2 -s100 {input.READ2_FILEPATH} {params.SUBSAMPLING} > {output.R2}
        fi
        """   
        
rule run_kb:   
    input:
        R1 = os.path.join(SCRATCH_DIR, '{dataset_project_id}/{dataset_sample_id}/tmp_{sub}/subsampledreads_{sub}_R1.fq'),
        R2 = os.path.join(SCRATCH_DIR, '{dataset_project_id}/{dataset_sample_id}/tmp_{sub}/subsampledreads_{sub}_R2.fq'),
    params: 
        DATASET_SAMPLE_ID = '{dataset_sample_id}',
        INDEX=make_index_path,
        TECHNOLOGY = fetch_technology,
        WHITELIST=make_whitelist_path,
        T2G=make_t2g_path,
        KB_OUT_DIR=directory('subsampling/{dataset_project_id}/{dataset_sample_id}/kbsub_{sub}/'),
    output:
         INSPECT_FILE='subsampling/{dataset_project_id}/{dataset_sample_id}/kbsub_{sub}/output/counts_unfiltered/adata.h5ad',
         
    shell:
        """
        echo "deeeeeeeeeeeeerp"
        echo "kb count --h5ad -i {params.INDEX} -g {params.T2G} -x {params.TECHNOLOGY} --overwrite --verbose --keep-tmp -t {THREAD} -o output {input.R1} {input.R2}"
        echo "faaaaaaaaaart"
        mkdir -p {params.KB_OUT_DIR}
        cd  {params.KB_OUT_DIR}
        kb count --h5ad -i {params.INDEX} -g {params.T2G} -x {params.TECHNOLOGY} --overwrite --verbose --keep-tmp -t {THREAD} -o output {input.R1} {input.R2} && \
        rm {input.R1} && rm {input.R2}
        """ 
