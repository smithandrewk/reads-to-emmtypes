import pandas as pd
import subprocess
from tqdm import tqdm
import os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_yellow(str):
    print(f'{bcolors.WARNING}{str}{bcolors.ENDC}')

def print_green(str):
    print(f'{bcolors.OKGREEN}{str}{bcolors.ENDC}')

def print_on_start_on_end(func):
    import time
    def inner1(*args, **kwargs):
        begin = time.time()
        print_yellow(f'Starting execution of {func.__name__}')
        rets = func(*args, **kwargs)
        print_green(f'Execution completed of {func.__name__}')
        end = time.time()
        print(f'Total time taken in {func.__name__}: {end - begin} s')
        return rets
    return inner1

def print_execution_time(func):
    import time
    def inner1(*args, **kwargs):
        begin = time.time()
        rets = func(*args, **kwargs)
        end = time.time()
        print(f'Total time taken in {func.__name__} [{args}]: {end - begin} s')
        return rets
    return inner1

@print_execution_time
def get_new_sequence_ids():
    # TODO: automatic sequence id logic?
    ids = ["STREP22-0001","STREP22-0002","STREP22-0003","STREP22-0004"]
    print(ids)
    return ids

@print_execution_time
def download_sample(biosample_id):
    print_yellow(f'downloading {biosample_id}')
    command = f'bs list datasets --input-biosample={biosample_id} -f csv | cut -d , -f2 | tail -n +2'

    stdout = subprocess.check_output(command,encoding='utf8',shell=True)
    datasets = stdout.split()
    info = pd.DataFrame()

    for dataset in datasets:
        command = f'bs get dataset -i {dataset} -f csv'
        stdout = subprocess.check_output(command,encoding='utf8',shell=True)
        info = pd.concat([info,pd.DataFrame([interior.split(',') for interior in stdout.split('\n')[:2]]).T.set_index(0).T]) 
    info = info[['Id','QcStatus','Project.DateModified']]
    dataset_id = info[((info['QcStatus'] != 'QcFailed') & (~info['QcStatus'].isna()))].sort_values('Project.DateModified',ascending=False).sort_values('Project.DateModified',ascending=False).iloc[0,0]
    command = f'bs download dataset --id {dataset_id} -v --no-metadata -o data/0_raw'
    stdout = subprocess.check_output(command,encoding='utf8',shell=True)
    print_green(f'downloaded {biosample_id}')
    return dataset_id

@print_on_start_on_end
def download_samples(biosample_ids):
    for biosample_id in tqdm(biosample_ids):
        download_sample(biosample_id)

@print_on_start_on_end
def trim_raw_reads():
    source_dir = f'data/0_raw'
    target_dir = f'data/1_trimmed'
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    for file in tqdm(os.listdir(source_dir)):
        basespace_id = file.split("_L001_")[0]
        if os.path.exists(f'{target_dir}/{basespace_id}_R1.trimd.fastq'):
            print_green(f'{basespace_id} already trimmed, continuing')
            continue
        print_yellow(f'trimming {basespace_id}')
        command = f'cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 20 --minimum-length 50 -o {target_dir}/{basespace_id}_R1.temp.fastq -p {target_dir}/{basespace_id}_R2.temp.fastq -j 0 {source_dir}/{basespace_id}_L001_R1_001.fastq.gz {source_dir}/{basespace_id}_L001_R2_001.fastq.gz >> cutadapt.log'
        os.system(command)
        command = f'cutadapt -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -q 20 --minimum-length 50 -o {target_dir}/{basespace_id}_R1.trimd.fastq -p {target_dir}/{basespace_id}_R2.trimd.fastq -j 0 {target_dir}/{basespace_id}_R1.temp.fastq {target_dir}/{basespace_id}_R2.temp.fastq >> cutadapt.log'
        os.system(command)
        os.system(f'rm {target_dir}/{basespace_id}_R1.temp.fastq {target_dir}/{basespace_id}_R2.temp.fastq')

@print_on_start_on_end
def qc_trimmed_reads():
    source_dir = f'data/1_trimmed'
    target_dir = f'data/2_qc'
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    command = f'fastqc {source_dir}/* -o {target_dir}'
    print(command)
    os.system(command)

def rename_raw():
    import os
    dir = f'samples'
    for file in os.listdir(f'{dir}'):
        for subfile in os.listdir(f'{dir}/{file}'):
            prefix = subfile.split('_')[-2]
            command = f'mv {dir}/{file}/{subfile} {dir}/{file}/{prefix}.fastq.gz'
            os.system(command)
            print(command)