"""
Downloads ENCODE datasets
"""

import sys
import pandas as pd
import numpy as np
import os
import tempfile

from datetime import datetime
import gzip
from itertools import islice
import time
import shutil

from urllib.request import urlretrieve
import hashlib
import multiprocessing

log = None
we_created_tempdir = False
if 'snakemake' in locals():

    log = open(str(snakemake.log), 'w')
    sys.stdout = sys.stderr = log

    temp_dir = snakemake.resources.tmpdir
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
        we_created_tempdir = True

    print(f"Setting temp dir to {temp_dir}")
    tempfile.tempdir = temp_dir

def md5_checksum(filename):
    # https://stackoverflow.com/a/3431838/171400
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

    
def batch(iterable, batch_size):
    # https://stackoverflow.com/a/62913856/171400
    iterator = iter(iterable)
    while batch := list(islice(iterator, batch_size)):
        yield batch

def download_file(url, checksum, output_location):

    # Inspired by SO:
    # https://stackoverflow.com/a/44378512/171400

    tmp_handle, tmp_loc = tempfile.mkstemp()
    # Close the temp file immediately to prevent "too many open files"
    os.close(tmp_handle)
    try:
        # Download file
        urlretrieve(url, tmp_loc)
        # Compute checksum
        checksum_of_downloaded = md5_checksum(tmp_loc)
        # Verify checksum
        if checksum_of_downloaded != checksum:
            raise Exception(f"Download for {url} failed checksum check {checksum_of_downloaded!r} != {checksum!r}")
        # Move file to correct location if checksum OK
        shutil.move(tmp_loc, output_location)
    finally:
        # Cleanup
        if os.path.isfile(tmp_loc):
            os.unlink(tmp_loc)


def _download_file(args):
    download_file(*args)

def main(input_tsv, output_dir, output_tsv, n_simultaneous=1, batch_size=100):

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    df_full = pd.read_csv(input_tsv, sep='\t')

    df = df_full[['File download URL', 'md5sum', 'Filename']].drop_duplicates()

    tasks = []

    for ix, row in df.iterrows():
        
        url = row['File download URL']
        checksum = row['md5sum']
        out = row['Filename']

        assert os.path.samefile(output_dir, os.path.dirname(out))
        tasks.append((url, checksum, out))

    n_tasks = len(tasks)
    if n_tasks > 0:
        first_task = tasks[0][2]
        last_task = tasks[-1][2]
    else:
        first_task = last_task = None
    
    start_time = datetime.now()
    print("{}: Processing {:,} tasks. First task: {}, Last task: {}".format(start_time.isoformat(), n_tasks, first_task, last_task))

    n_done = 0
    with multiprocessing.Pool(n_simultaneous) as pool:

        for task_batch in batch(tasks, batch_size=batch_size):
            
            for __ in pool.imap_unordered(_download_file, task_batch):
                n_done += 1
            
            print("{}: Done {:,}/{:,} ({:.2%})".format(datetime.now().isoformat(), n_done, n_tasks, n_done/n_tasks))
            if log is not None:
                # Flush python and os buffers to write instantly
                log.flush()
                os.fsync(log.fileno())
    
    end_time = datetime.now()
    print("{}: Done (took {} s)".format(end_time.isoformat(), (end_time-start_time).total_seconds()))
    
    # We're done, copy the full TSV into output file
    print(f"Copying input tsv to {output_tsv}")
    df_full.to_csv(output_tsv, sep='\t', index=False)

if __name__ == '__main__' and 'snakemake' in locals():

    try:
        main(
            snakemake.input.tsv,
            snakemake.output.files,
            snakemake.output.tsv,
            snakemake.threads,
            snakemake.params.get('batch_size', 100),
        )
    finally:
        if log is not None:
            # Flush python and os buffers to write instantly
            log.flush()
            os.fsync(log.fileno())

        # Cleanup the temp directory if it's empty
        if we_created_tempdir or not any(os.scandir(temp_dir)):
            shutil.rmtree(temp_dir)

        time.sleep(60)