import pandas as pd
import subprocess
# bs list biosamples -f csv > biosamples.csv

def rename_raw():
    import os
    dir = f'samples'
    for file in os.listdir(f'{dir}'):
        for subfile in os.listdir(f'{dir}/{file}'):
            prefix = subfile.split('_')[-2]
            command = f'mv {dir}/{file}/{subfile} {dir}/{file}/{prefix}.fastq.gz'
            os.system(command)
            print(command)