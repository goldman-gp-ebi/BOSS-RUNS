import multiprocessing
from itertools import repeat
import os

# non-std lib
import numpy as np


'''
sample reads from a large fastq file
- used in simulations
- superseeded by data sampler module

'''



def parallel_batches(fq, bsize, workers):
    '''
    simulations: sample random reads from a large fastq file in parallel

    Parameters
    ----------
    fq: str
        path to fastq file for sampling
    bsize: int
        number of reads to sample
    workers: int
        number of parallel workers/threads

    Returns
    -------
    read_lengths: dict
        dict of read_id: len
    read_sequences: dict
        dict of read_id: sequence
    basesTOTAL: int
        total number of bases in this batch, used for recording time

    '''
    # initialise pool for parallelisation
    pool = multiprocessing.Pool(processes=workers)
    # prep arg lists
    fqs = [fq] * workers
    seeds = np.random.randint(10000, size=workers)
    # launch workers
    batches = pool.starmap(_get_random_batch, zip(fqs, repeat(int(bsize / workers)), seeds))
    # join the results from all workers
    batches_concat = ''.join(batches)
    # parse the batch (long string) into dicts
    read_lengths, read_sequences, basesTOTAL = _parse_batch(batch_string=batches_concat)
    return read_lengths, read_sequences, basesTOTAL


def _get_random_batch(fq, bsize, seed):
    batch = ''
    np.random.seed(seed)

    for b in range(bsize):
        read = _get_random_read(filepath=fq)
        batch += read
    return batch


def _get_random_read(filepath):
    file_size = os.path.getsize(filepath)
    fq_lines = b''
    proper_line = b''
    # open file in binary mode
    with open(filepath, 'rb') as f:
        while True:
            pos = np.random.randint(1, file_size)
            f.seek(pos)  # seek to random position
            f.readline()  # skip possibly incomplete line
            # skip lines until we are at a fastq header
            while not proper_line.decode().startswith('@'):
                proper_line = f.readline()
            # add the fastq header to the read
            fq_lines += proper_line
            for i in range(3): # add the other 3 fastq lines as well
                fq_lines += f.readline()
            if fq_lines:
                return fq_lines.decode()


def _parse_batch(batch_string):
    # take a batch in string format, parse it into containers
    # Imitates reading from actual fq
    read_lengths = {}
    read_sequences = {}

    batch_lines = batch_string.split('\n')
    n_lines = len(batch_lines)

    i = 0
    # since we increment i by 4 (lines for read in fq), loop until n_lines - 4
    while i < (n_lines - 4):
        # grab the name of the read. split on space, take first element, trim the @
        name = batch_lines[i].split(' ')[0][1:]
        seq = batch_lines[i + 1]
        read_len = len(seq)
        # update the containers
        read_lengths[name] = read_len
        read_sequences[name] = seq
        i += 4
    # get the total length of reads in this batch
    basesTOTAL = np.sum(np.array(list(read_lengths.values())))
    return read_lengths, read_sequences, basesTOTAL

