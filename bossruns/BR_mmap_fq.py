import logging
import time
import os
import gzip
import mmap

# non-std lib
import numpy as np

# custom imports
from BR_sample_fq import _parse_batch




class FastqStream:
    """
    Stream reads from a fastq file (4 lines per read).
    used for random sampling in simulations
    """
    def __init__(self, source, log_each=int(1e5)):
        self.source = source  # path to file-like object

        # check if file is gzipped
        suffix = source.split('.')[-1]
        if suffix == 'gz':
            self.gzipped = True
        else:
            self.gzipped = False

        self.log_each = int(log_each)  # defining logging frequency of scanning for offsets
        self.filesize = int(os.stat(source).st_size)
        print(f"Representing {self.filesize / 1e6} Mbytes of data from source: {self.source}")


    def scan_offsets(self, k=4, limit=1e9):
        '''
        Scan file to find byte offsets of reads.
        Offsets are created for chunks of k lines each (4 for fastq)

        Parameters
        ----------
        k: int
            number of lines per read in fq format
        limit: int
            maximum number of offsets to scan from a file

        Returns
        -------

        '''
        tic = time.time()
        tmp_offsets = []
        read_num = 0

        with open(self.source, 'rb') as f:
            k_tmp = 1
            # memory-map the file; lazy eval-on-demand via POSIX filesystem
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            # handle gzipped
            if self.gzipped:
                mm = gzip.GzipFile(mode="rb", fileobj=mm)

            for _ in iter(mm.readline, b''):
                if k_tmp % k == 0:
                    pos = mm.tell()
                    tmp_offsets.append(pos)
                    k_tmp = 1
                    read_num += 1
                    # status update in case there are many reads
                    if read_num % self.log_each == 0:
                        print(f"{read_num} reads scanned")
                else:
                    k_tmp += 1

                if read_num >= limit:
                    break

        toc = time.time()
        # uint32 enough for most files
        offsets = np.asarray(tmp_offsets, dtype='uint64')
        del tmp_offsets
        # write the offsets to a file
        np.save(f'{self.source}.offsets', offsets)
        print(f"DONE scanning {read_num} reads")
        print(f'wrote {len(offsets)} offsets to {self.source}.offsets.npy')
        print(f"{round(toc - tic, 4)} seconds elapsed scanning file for offsets")


    def load_offsets(self, seed, shuffle=False, batchsize=1, maxbatch=1):
        if seed == 0:
            seed = np.random.randint(low=0, high=int(1e6))
        np.random.seed(seed)
        offsets = np.load(f'{self.source}.offsets.npy')
        # add one batch for initialising length dist
        maxbatch = maxbatch + 1

        if shuffle:
            np.random.shuffle(offsets)
            print(f"offsets shuffled using random seed: {seed}")

        # shorten the offsets to number of reads we need
        len_offsets = len(offsets)
        n_reads = batchsize * maxbatch
        if n_reads < len_offsets:
            offsets = offsets[: n_reads]
        else:
            print("requested more reads than there are available in the fastq")
            # sys.exit()

        # restructure the offsets into 2D array to represent batches (rows)
        offsets = offsets.reshape((maxbatch, batchsize))
        self.offsets = offsets


    def get_batch(self):
        '''
        return a batch of reads from the fastq file

        Returns
        -------
        read_lengths: dict
            dict of read_ids: lengths
        read_sequences: dict
            dict of read_ids: sequences
        basesTOTAL: int
            total number of sequenced bases, used for time keeping

        '''
        batch = ''
        # check if offsets are empty
        if self.offsets.shape[0] == 0:
            logging.info("no more batches left")
            # sys.exit()

        with open(self.source, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            # handle gzipped files
            if self.gzipped:
                mm = gzip.GzipFile(mode="rb", fileobj=mm)
            # the first row of the offsets are the next batch
            batch_offsets = self.offsets[0, :]
            # this is probably LINUX specific and not POSIX
            # tell the kernel what to do with the mapped memory
            # preload specific pages of the file
            # which magically makes it even faster to access
            # pagesize is LINUX constant of 4096 bytes per page
            pagesize = 4096 #
            # the start of WILLNEED needs to be a multiple of pagesize
            # so we take the modulo and move the start of the offset a little bit earlier if needed
            new_offsets = batch_offsets - (batch_offsets % pagesize)

            for new_offset in new_offsets:
                # preload 20 pages of data following each read start
                # 20 pages = 80 kbytes (covers read of up to ~40 kbases)
                mm.madvise(mmap.MADV_RANDOM)
                mm.madvise(mmap.MADV_WILLNEED, int(new_offset), 20)

            batch_offsets = np.sort(batch_offsets)
            for offset in batch_offsets:
                try:
                    # Use offset to jump to position in the file
                    # return the next 4 lines, i.e. one read
                    chunk = FastqStream.get_single_read(mm=mm, offset=offset)
                    # append the fastq entry to the batch
                    batch += chunk
                except:
                    print(f"Error at location: {offset}")
                    continue
                if len(chunk) == 0:
                    continue

            # add call to close memory map, only file itself is under with()
            mm.close()

        if not batch.startswith('@'):
            print("The batch is broken")


        # remove the row from the offsets so it does not get sampled again
        new_offsets = np.delete(self.offsets, 0, 0)
        self.offsets = new_offsets
        # parse the batch, which is just a long string into dicts
        read_lengths, read_sequences, basesTOTAL = _parse_batch(batch_string=batch)
        return read_lengths, read_sequences, basesTOTAL


    @staticmethod
    def get_single_read(mm, offset):
        # return 4 lines from a memory-mapped fastq file given a byte-wise position
        mm.seek(offset)
        chunk_size = 4
        chunk = b''
        # read the 4 lines of the fastq entry
        for _ in range(chunk_size):
            chunk += mm.readline()

        chunk = chunk.decode("utf-8")
        return chunk


