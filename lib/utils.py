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
        command = f'staphb-tk trimmomatic PE {source_dir}/{basespace_id}_L001_R1_001.fastq.gz {source_dir}/{basespace_id}_L001_R2_001.fastq.gz {target_dir}/{basespace_id}_R1.trimd.fastq {target_dir}/{basespace_id}_R2.trimd.fastq LEADING:20 TRAILING:20 MINLEN:50'
        os.system(command)

@print_on_start_on_end
def qc_trimmed_reads():
    source_dir = f'data/1_trimmed'
    target_dir = f'data/2_qc'
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    files = [file.replace(".fastq","") for file in os.listdir(source_dir)]
    for file in files:
        if os.path.exists(f'{target_dir}/{file}_fastqc.html'):
            print(f'{file} already qc, continue')
            continue
        command = f'staphb-tk fastqc {source_dir}/{file}.fastq -o {target_dir}'
        os.system(command)
    command = f'staphb-tk multiqc data/2_qc -o data/3_mqc --force'
    os.system(command)

def assemble_reads():
    import os
    source_dir = f'data/1_trimmed'
    target_dir = f'data/4_assembled'
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    files = [file.replace("_R1.trimd.fastq","") for file in os.listdir(source_dir) if file.endswith('R1.trimd.fastq')]
    for file in files:
        print(file)
        if os.path.isdir(f'{target_dir}/{file}'):
            print(f'already assembled {file}')
            continue
        command = f'staphb-tk spades -1 {source_dir}/{file}_R1.trimd.fastq -2 {source_dir}/{file}_R2.trimd.fastq -o {target_dir}/{file}'
        os.system(command)
    command = f'staphb-tk quast data/4_assembled/*/contigs.fasta -o data/5_quast -t 1'
    os.system(command)

def emmtype_assemblies():
    import os
    source_dir = f'data/4_assembled'
    target_dir = f'data/6_emmtyped'
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    for file in os.listdir(source_dir):
        print(file)
        if os.path.isdir(f'{target_dir}/{file}'):
            print(f'already emmtyped {file}')
            continue
        command = f'emmtyper {source_dir}/{file}/contigs.fasta -o {target_dir}/{file}.tsv'
        print(command)
        os.system(command)
    command = f'emmtyper {source_dir}/*/contigs.fasta -o {target_dir}/all_emmtypes.tsv'
    os.system(command)
