import pickle
import os
import logging
from itertools import repeat
import sys
from concurrent.futures import ThreadPoolExecutor as TPE

# non-std lib
import numpy as np



class Sampler:

    def __init__(self, fq, paf_full, paf_trunc, bsize, workers):
        self.fq = fq
        self.paf_full = paf_full
        self.paf_trunc = paf_trunc
        self.bsize = bsize
        self.workers = workers

        # load offsets of fq
        offsets_fq = self._load_offsets_fq(path=fq)
        self.offsets_fq = offsets_fq

        # load offsets for mappings
        offsets_full = self._load_offsets_paf(path=paf_full)
        self.offsets_full = offsets_full

        offsets_trunc = self._load_offsets_paf(path=paf_trunc)
        self.offsets_trunc = offsets_trunc


    def _load_offsets_paf(self, path):
        # check if offsets present
        if not os.path.exists(f'{path}.offsets'):
            print("no paf offsets found. exiting")
            return

        with open(f'{path}.offsets', "rb") as p:
            offsets = pickle.load(p)
        logging.info(f"loaded offsets for {path.split('/')[-1]}")
        return offsets


    def _load_offsets_fq(self, path):
        # check if offsets present
        if not os.path.exists(f'{path}.offsets.npy'):
            print("no fq offsets found. exiting")
            return

        offsets = np.load(f'{path}.offsets.npy')
        logging.info(f"loaded offsets for {path.split('/')[-1]}")
        return offsets


    def _grab_batch(self, pos):
        # helper func for multiprocc
        batch = ''

        for n in range(len(pos)):
            read = self._grab_read(pos[n])
            batch += read
        return batch


    def _grab_read(self, p):
        # get read from specific position
        fq_lines = b''
        proper_line = b''

        with open(self.fq, 'rb') as f:
            while True:
                f.seek(p)  # seek to random position
                f.readline()  # skip possibly incomplete line
                # skip lines until we are at a fastq header
                while not proper_line.decode().startswith('@'):
                    proper_line = f.readline()
                # add the fastq header to the read
                fq_lines += proper_line
                for i in range(3):  # add the other 3 fastq lines as well
                    fq_lines += f.readline()
                if fq_lines:
                    return fq_lines.decode()


    @staticmethod
    def _grab_paf_lines(args):
        # helper for multithread
        path, pos = args

        paf = b''
        with open(path, 'rb') as f:
            for p in pos:
                f.seek(p)  # seek to position
                paf += f.readline()

        return paf.decode()


    def parallel_batch(self, batch):
        # random sample from offsets
        # positions = np.sort(np.random.choice(a=self.offsets_fq, size=self.bsize, replace=False))
        # instead of random we can also just slice
        used_reads = batch * self.bsize
        positions = self.offsets_fq[used_reads : used_reads + self.bsize]
        if positions.shape[0] == 0:
            logging.info("fq offsets empty. all reads used up")
            sys.exit()

        worker_positions = np.array_split(positions, self.workers)
        # distribute work
        with TPE(max_workers=self.workers) as executor:
            batches = executor.map(self._grab_batch, worker_positions)

        # join results from workers
        batches_concat = ''.join(batches)

        # parse the batch, which is just a long string into dicts
        read_lengths, read_sequences, read_qualities, basesTOTAL = parse_batch(batch_string=batches_concat)
        # self.read_ids = list(read_lengths.keys())
        return read_lengths, read_sequences, read_qualities, basesTOTAL


    def grab_mappings(self, read_ids):
        # given some sampled read ids, grab their full and trunc mappings
        # first get their offsets
        rid_offsets_full = [self.offsets_full[rid] for rid in read_ids]
        unmapped = np.sum([1 for i in rid_offsets_full if len(i) == 0])
        mapped = np.sum([1 for i in rid_offsets_full if len(i) > 0])
        logging.info(f"full: mapped {mapped}, unmapped {unmapped}")
        # transform offsets into array
        offsets_arr_full = np.concatenate([np.array(i, dtype=int) for i in rid_offsets_full])

        rid_offsets_trunc = [self.offsets_trunc[rid] for rid in read_ids]
        unmapped = np.sum([1 for i in rid_offsets_trunc if len(i) == 0])
        mapped = np.sum([1 for i in rid_offsets_trunc if len(i) > 0])
        logging.info(f"trunc: mapped {mapped}, unmapped {unmapped}")
        # transform offsets into array
        offsets_arr_trunc = np.concatenate([np.array(i, dtype=int) for i in rid_offsets_trunc])

        # grab the pafs
        worker_offsets = np.array_split(offsets_arr_full, self.workers)
        with TPE(max_workers=self.workers) as executor:
            pafs = executor.map(self._grab_paf_lines,
                                zip(repeat(self.paf_full),
                                    worker_offsets))
        paf_full = ''.join(pafs)

        # paf_full = self._grab_paf_lines(self.paf_full, pos=offsets_arr_full)
        worker_offsets = np.array_split(offsets_arr_trunc, self.workers)
        with TPE(max_workers=self.workers) as executor:
            pafs = executor.map(self._grab_paf_lines,
                                zip(repeat(self.paf_trunc),
                                    worker_offsets))
        paf_trunc = ''.join(pafs)

        # paf_trunc = self._grab_paf_lines(self.paf_trunc, pos=offsets_arr_trunc)
        return paf_full, paf_trunc



def parse_batch(batch_string):
    # parse batch from string format into containers
    read_lengths = {}
    read_sequences = {}
    read_qualities = {}

    batch_lines = batch_string.split('\n')
    n_lines = len(batch_lines)

    i = 0
    # since we increment i by 4 (lines for read in fq), loop until n_lines - 4
    while i < (n_lines - 4):
        # grab the name of the read. split on space, take first element, trim the @
        name = batch_lines[i].split(' ')[0][1:]
        seq = batch_lines[i + 1]
        qual = batch_lines[i + 3]
        read_len = len(seq)
        # update the containers
        read_lengths[name] = read_len
        read_sequences[name] = seq
        read_qualities[name] = qual
        i += 4
    # get the total length of reads in this batch
    basesTOTAL = np.sum(np.array(list(read_lengths.values())))
    return read_lengths, read_sequences, read_qualities, basesTOTAL


