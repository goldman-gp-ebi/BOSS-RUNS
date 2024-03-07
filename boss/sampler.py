import logging
import time
import mmap
import gzip
import pickle
from pathlib import Path
from typing import Any
from itertools import repeat
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor as TPE

import numpy as np
from numpy.typing import NDArray


# module for random read sampling in BOSS


class Sampler:

    def __init__(self, source: str, paf_full: str = None, paf_trunc: str = None, workers: int = 4, **kwargs) -> None:
        """
        Wrapper to sample sequencing data (and their mappings) given source files.
        Paf files and offsets can be generated with a script in ../scripts

        :param source: Fastq files with sequencing data to sample from
        :param paf_full: Full-length mappings of the same reads
        :param paf_trunc: Truncated (mu-sized) mappings of the same reads
        :param workers: Max number of threads for mapping retrieval
        :param kwargs: Pass-through args to the FastqStream
        """
        self.fq_stream = FastqStream_mmap(source=source, **kwargs)
        # stream for mapping data is only initialised for BOSS-RUNS, not AEONS
        self.pafs = True if paf_full and paf_trunc else False
        if self.pafs:
            self.paf_stream = PafStream(paf_full=paf_full, paf_trunc=paf_trunc, workers=workers)



    def sample(self) -> tuple[dict[str, str], dict[str, str], str, str]:
        """
        Generate a new batch of data

        :return: Tuple of read sequences, qualities and (optional) raw string output of mappings
        """
        self.fq_stream.read_batch()
        if self.pafs:
            paf_f, paf_t = self.paf_stream.grab_mappings(read_ids=self.fq_stream.read_ids)
        else:
            paf_f = defaultdict(list)
            paf_t = defaultdict(list)
        return self.fq_stream.read_sequences, self.fq_stream.read_qualities, paf_f, paf_t



class FastqStream_mmap:

    def __init__(self, source: str, seed: int = 1, shuffle: bool = False, batchsize: int = 1, maxbatch: int = 1):
        """
        Stream reads from a fastq file (4 lines per read) during simulations
        This implementation won't sample the same read twice

        :param source: Path to the file-like object.
        :param seed: Seed for random number generation.
        :param shuffle: Whether to shuffle the offsets.
        :param batchsize: Batch size.
        :param maxbatch: Maximum number of batches.
        """
        self.source = source
        if not self.source:
            raise AssertionError("No source file provided")
        # check if file is gzipped. Not very good
        suffix = source.split('.')[-1]
        if suffix == 'gz':
            self.gzipped = True
        else:
            self.gzipped = False

        self.log_each = int(int(1e5))  # log frequency of scanning for offsets
        self.filesize = int(Path(source).stat().st_size)
        logging.info(f"Representing {self.filesize / 1e6} Mbytes of data from source: {self.source}")
        # scan the offsets if the file does not exist
        if not Path(f'{source}.offsets.npy').exists():
            logging.info("scanning offsets")
            self._scan_offsets()
        logging.info("loading offsets")
        self._load_offsets(seed=seed, shuffle=shuffle, batchsize=batchsize, maxbatch=maxbatch)
        self.batch = 0



    def _scan_offsets(self, k: int = 4, limit: float = 1e9) -> None:
        """
        Scan the file to find byte offsets of each sequencing read.

        :param k: Number of lines per chunk, default chink size is 4 for fastq
        :param limit: Maximum number of reads to scan.
        :return: None.
        """
        tic = time.time()
        tmp_offsets = []
        read_num = 0

        with open(self.source, 'rb') as f:
            k_tmp = 1
            # memory-map the file; lazy eval-on-demand via POSIX filesystem
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            if self.gzipped:
                mm = gzip.GzipFile(mode="rb", fileobj=mm)

            for _ in iter(mm.readline, b''):
                if k_tmp % k == 0:
                    pos = mm.tell()
                    tmp_offsets.append(pos)
                    k_tmp = 1
                    read_num += 1
                    # status update
                    if read_num % self.log_each == 0:
                        logging.info(f"{read_num} reads scanned")
                else:
                    k_tmp += 1

                if read_num >= limit:
                    break

        toc = time.time()
        offsets = np.asarray(tmp_offsets, dtype='uint64')
        del tmp_offsets
        # write the offsets to a file
        np.save(f'{self.source}.offsets', offsets)
        logging.info(f"DONE scanning {read_num} reads")
        logging.info(f'wrote {len(offsets)} offsets to {self.source}.offsets.npy')
        logging.info(f"{round(toc - tic, 4)} seconds elapsed scanning file for offsets")



    def _load_offsets(self, seed: int = 1, shuffle: bool = False, batchsize: int = 1, maxbatch: int = 1) -> None:
        """
        Load offsets of sequencing reads in a fastq file.

        :param seed: Seed for random number generation.
        :param shuffle: Whether to shuffle the offsets.
        :param batchsize: Batch size.
        :param maxbatch: Maximum number of batches.
        :return: None.
        """
        if seed == 0:
            seed = np.random.randint(low=0, high=int(1e6))
        np.random.seed(seed)
        offsets = np.load(f'{self.source}.offsets.npy')
        maxbatch = maxbatch + 1

        if shuffle:
            np.random.shuffle(offsets)
            logging.info(f"offsets shuffled using random seed: {seed}")

        # shorten the offsets to number of reads we need
        len_offsets = len(offsets)
        logging.info(f"available batches: {len_offsets / batchsize}")
        n_reads = batchsize * maxbatch
        if n_reads < len_offsets:
            offsets = offsets[: n_reads]
        else:
            raise ValueError("Requested more reads than there are available in the fastq")

        # restructure the offsets into 2D array to represent batches (rows)
        offsets = offsets.reshape((maxbatch, batchsize))
        self.offsets = offsets



    def _parse_batch(self, batch_string: str) -> None:
        """
        Parse a batch in string format into dictionaries.

        :param batch_string: Batch in string format.
        :return: Tuple containing dictionaries of read lengths, sequences, sources, and the total number of bases.
        """
        read_lengths = {}
        read_sequences = {}
        read_qualities = {}
        read_sources = {}

        batch_lines = batch_string.split('\n')
        n_lines = len(batch_lines)

        i = 0
        # since we increment i by 4 (lines for read in fq), loop until n_lines - 4
        while i < (n_lines - 4):
            # grab the name of the read. split on space, take first element, trim the @
            desc = batch_lines[i].split(' ')
            name = str(desc[0][1:])
            source = str(desc[-1])
            seq = batch_lines[i + 1]
            qual = batch_lines[i + 3]
            read_len = len(seq)
            # update the containers
            read_lengths[name] = read_len
            read_sequences[name] = seq
            read_sources[name] = source
            read_qualities[name] = qual
            i += 4
        # get the total length of reads in this batch
        total_bases = int(np.sum(np.array(list(read_lengths.values()))))
        # assign attributes
        self.read_ids = set(read_sequences.keys())
        self.read_lengths = read_lengths
        self.read_sequences = read_sequences
        self.read_qualities = read_qualities
        self.read_sources = read_sources
        self.total_bases = total_bases



    @staticmethod
    def _get_single_read(mm: Any, offset: int) -> str:
        """
        Return 4 lines from a memory-mapped fastq file given a byte-wise position.

        :param mm: The memory-mapped fastq file.
        :param offset: The byte-wise position of the read.
        :return: The read chunk.
        """
        mm.seek(offset)
        chunk_size = 4
        chunk = b''
        # read 4 lines of the fastq entry
        for _ in range(chunk_size):
            chunk += mm.readline()
        chunk = chunk.decode("utf-8")
        return chunk



    def read_batch(self, delete: bool = True) -> None:
        """
        Return a batch of reads from the fastq file.

        :param delete: Whether to delete the batch offsets after retrieval.
        :return:
        """
        # check if offsets are empty
        if self.offsets.shape[0] == 0:
            raise ValueError("No more reads left to sample. Exiting")

        with open(self.source, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            if self.gzipped:
                mm = gzip.GzipFile(mode="rb", fileobj=mm)
            # the first row of the offsets are the next batch
            batch_offsets = self.offsets[0, :]
            # initialise list instead of string concat
            batch = [''] * len(batch_offsets)
            # possibly LINUX specific and not POSIX
            # here we tell the kernel to preload pages of the mapped memory
            # magically makes it even faster to access
            # pagesize is a LINUX (system)-specific constant of 4096 bytes per "page"
            pagesize = 4096
            # the start of WILLNEED needs to be a multiple of pagesize
            # we take the modulo and move the start of the offset a bit earlier if needed
            new_offsets = batch_offsets - (batch_offsets % pagesize)

            if not self.gzipped:
                for new_offset in new_offsets:
                    # we preload 20 pages of data following each read start
                    # 20 pages = 80 kbytes (read of up to ~40 kbases, I think)
                    mm.madvise(mmap.MADV_RANDOM)
                    mm.madvise(mmap.MADV_WILLNEED, int(new_offset), 20)

            batch_offsets = np.sort(batch_offsets)
            for i in range(len(batch_offsets)):
                try:
                    # jump to position in file and return the next 4 lines
                    chunk = self._get_single_read(mm=mm, offset=int(batch_offsets[i]))
                    batch[i] = chunk
                except:
                    logging.info(f"Error at location: {batch_offsets[i]}")
                    continue
                if len(chunk) == 0:
                    continue

            # add call to close memory map, only file itself is under with()
            mm.close()

        if not batch[0].startswith('@') and not batch[0].startswith('>'):
            raise ValueError("The batch of reads is broken")

        if delete:
            # remove the row from the offsets, so it does not get sampled again
            new_offsets = np.delete(self.offsets, 0, 0)
            self.offsets = new_offsets
        # parse the batch (long string) into dicts
        batch_string = ''.join(batch)
        self._parse_batch(batch_string=batch_string)
        self.batch += 1
        logging.info(f'got new batch of {len(self.read_sequences)} reads')





class PafStream:


    def __init__(self, paf_full: str, paf_trunc: str, workers: int = 4):
        """
        Initialise a PafStream object to quickly grab paf entries given sampled read IDs

        :param paf_full: Path to paf file with mappings of full-length reads
        :param paf_trunc: Path to paf file with mappings of truncated reads
        :param workers: Max threads to spawn
        """

        self.paf_full = paf_full
        self.paf_trunc = paf_trunc
        self.workers = workers

        # scan the offsets if the file does not exist
        if not Path(f'{paf_full}.offsets').exists():
            logging.info("scanning full-length PAF offsets")
            self._scan_offsets_paf(path=self.paf_full)
        self.offsets_full = self._load_offsets_paf(path=paf_full)

        # same for truncated mappings
        if not Path(f'{paf_trunc}.offsets').exists():
            logging.info("scanning truncated PAF offsets")
            self._scan_offsets_paf(path=self.paf_trunc)
        self.offsets_trunc = self._load_offsets_paf(path=paf_trunc)
        self.batch = 0



    @staticmethod
    def _load_offsets_paf(path: str) -> defaultdict[list]:
        """
        Load offsets from file

        :param path: Path to offsets file
        :return: Defaultdict with list of mappings of reads
        """
        # check if offsets present
        if not Path(f'{path}.offsets').is_file():
            print("no paf offsets found. exiting")
            exit(1)

        with open(f'{path}.offsets', "rb") as p:
            offsets = pickle.load(p)
        logging.info(f"loaded offsets for {path.split('/')[-1]}")
        return offsets



    @staticmethod
    def _scan_offsets_paf(path: str, limit: int = 1e9) -> None:
        """
        Scan file to find byte offsets of mappings

        :param path: Path to paf mapping file
        :param limit: maximum lines to scan (debugging)
        :return:
        """
        mapping_num = 0
        offsets = defaultdict(list)
        pos = 0

        with open(path, 'rb') as paf:
            for line in paf:
                ll = line.split(b"\t")
                read_id = ll[0].decode()
                offsets[read_id].append(pos)
                pos = paf.tell()
                mapping_num += 1

                if mapping_num % 1e5 == 0:
                    print(f"{mapping_num} mappings scanned")

                if mapping_num >= limit:
                    break

        # write the offsets to a file
        with open(f'{path}.offsets', 'wb') as p:
            pickle.dump(offsets, p)
        print(f"DONE scanning {mapping_num} mapping records")



    @staticmethod
    def _grab_paf_lines(args: tuple[str, NDArray]) -> str:
        """
        Helper func for multithreaded grabbing of paf lines

        :param args: tuple of file path and array of byte offsets
        :return:
        """
        path, pos = args
        paf = b''
        with open(path, 'rb') as f:
            for p in pos:
                f.seek(p)  # seek to position
                paf += f.readline()
        return paf.decode()



    @staticmethod
    def _get_offset_array(offsets: defaultdict[list], read_ids: set) -> NDArray:
        """
        Given a set of read IDs, get the array of offsets of the mappings of those reads

        :param offsets: Defaultdict of offsets of mappings
        :param read_ids: Set of reads from which to retrieve offsets
        :return: Array of offsets of mappings of requested reads
        """
        rid_offsets = [offsets[rid] for rid in read_ids]
        unmapped = np.sum([1 for i in rid_offsets if len(i) == 0])
        mapped = np.sum([1 for i in rid_offsets if len(i) > 0])
        logging.info(f"mapped {mapped}, unmapped {unmapped}")
        # transform offsets into array
        offsets_arr = np.concatenate([np.array(i, dtype=int) for i in rid_offsets])
        return offsets_arr


    def _grab_pafs(self, offsets: NDArray, path: str) -> str:
        """
        Multithreaded retrieval of paf mappings given an array of offsets

        :param offsets: Array of offsets of a specific set of reads
        :param path: Path to the original file to grab from
        :return: Concatenated mappings of interest
        """
        worker_offsets = np.array_split(offsets, self.workers)
        with TPE(max_workers=self.workers) as executor:
            pafs = executor.map(self._grab_paf_lines,
                                zip(repeat(path),
                                    worker_offsets))
        paf_out = ''.join(pafs)
        return paf_out


    def grab_mappings(self, read_ids: set) -> tuple[str, str]:
        """
        Given a set of read IDs, get their offsets and grab the mappings data from these reads

        :param read_ids: Set of read IDs
        :return: Tuple of mappings for full and truncated mappings of the reads
        """
        # given some sampled read ids, grab their full and trunc mappings
        offsets_arr_full = self._get_offset_array(
            offsets=self.offsets_full, read_ids=read_ids)
        offsets_arr_trunc = self._get_offset_array(
            offsets=self.offsets_trunc, read_ids=read_ids)
        # grab the actual data
        paf_f = self._grab_pafs(
            offsets=offsets_arr_full, path=self.paf_full)
        paf_t = self._grab_pafs(
            offsets=offsets_arr_trunc, path=self.paf_trunc)
        return paf_f, paf_t

