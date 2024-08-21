import logging
import re
import os
from pathlib import Path

import mappy
import numpy as np

from boss.utils import empty_file, execute, random_id
from boss.paf import paf_dict_type



class FastqBatch:

    def __init__(self, fq_files: list[str | Path], channels: set = None):
        """
        Initialise a new batch of sequencing reads using their filepaths

        :param fq_files: Paths of new files
        :param channels: Set of channels of the BOSS region
        """
        self.fq_files = fq_files
        self.channels = channels
        # attributes to be filled
        self.read_sequences = {}
        self.read_qualities = {}
        self.read_ids = {}
        self.read_lengths = {}
        self.total_bases = 0
        self._read_batch()



    def _read_batch(self) -> None:
        """
        Read sequencing data from all new files

        :return:
        """
        read_sequences = {}
        read_qualities = {}
        for fq in self.fq_files:
            rseq, rqual = self._read_single_batch(fastq_file=fq)
            read_sequences.update(rseq)
            read_qualities.update(rqual)
        # set attributes of the batch
        self.read_sequences = read_sequences
        self.read_qualities = read_qualities
        self.read_ids = set(read_sequences.keys())
        self.read_lengths = {rid: len(seq) for rid, seq in read_sequences.items()}
        self.total_bases = np.sum(list(self.read_lengths.values()))
        logging.info(f'total new reads: {len(read_sequences)}')



    def _read_single_batch(self, fastq_file: str) -> tuple[dict[str, str], dict[str, str]]:
        """
        Get the reads from a single fq file and put into dictionary

        :param fastq_file: A new fastq files
        :return: Dictionary of new sequences as header, sequence pairs
        """
        logging.info(f"Reading file: {fastq_file}")
        read_sequences = {}
        read_qualities = {}
        # to make sure it's a string
        if not isinstance(fastq_file, str):
            raise TypeError('Fastq file must be a string')
        # loop over all reads in the fastq file
        # if we consider all channels
        if not self.channels:
            for name, seq, qual, desc in mappy.fastx_read(fastq_file, read_comment=True):
                read_sequences[str(name)] = seq
                read_qualities[str(name)] = qual
        else:
            # consider source channel
            for name, seq, qual, desc in mappy.fastx_read(fastq_file, read_comment=True):
                try:
                    # regex to get the channel number from the header
                    # \s=whitespace followed by 'ch=' and then any amount of numeric characters
                    curr_channel = re.search("\sch=[0-9]*", desc).group()
                    ch_num = int(curr_channel.split('=')[1])
                except AttributeError:
                    # if the pattern is not in the header, skip the read
                    logging.info("ch= not found in header of fastq read")
                    continue
                # check if read comes from a BOSS channel
                if ch_num in self.channels:
                    read_sequences[str(name)] = seq
                    read_qualities[str(name)] = qual
        return read_sequences, read_qualities



class ReadCache:

    def __init__(self, batchsize: int, dumptime: int, alpha: int = 300, rho: int = 300, mu: int = 400):
        """
        Object for simulations only. Collects the pseudotime of different conditions and
        the reads seen by each condition. Then writes the observed reads to file in
        cumulative batches of data for downstream analysis. I.e. the files produced by this
        object can be used as timepoints in downstream analyses

        :param batchsize: How many reads are sampled in each batch
        :param dumptime: Time interval of writing data to file
        :param alpha: Acquisition time of model
        :param rho: Rejection time cost in model
        :param mu: Length of initial bit in model
        """
        self.alpha, self.rho, self.mu = alpha, rho, mu
        self.batchsize = batchsize
        # pseudotime counter
        self.time_boss = 0
        self.time_control = 0
        # for writing reads to file
        self.cache_control = dict()
        self.cache_boss = dict()
        # after how much time should sequences be written to file
        # dump time is incremented every time a batch is written, which happens once that is overcome
        self.dumptime = dumptime
        self.dump_n_control = 1
        self.dump_n_boss = 1
        # for storing the batches of reads for snakemake analysis
        if not os.path.exists('./00_reads'):
            os.mkdir('./00_reads')
        empty_file(f'00_reads/control_0.fa')
        empty_file(f'00_reads/boss_0.fa')




    def update_times_runs(self, total_bases: int, paf_dict: paf_dict_type, n_unmapped: int, n_reject: int) -> None:
        """
        Increment the pseudotime for control and boss regions on flowcell
        depending on the observed data and decisions made on them

        :param total_bases: Total observed bases in batch
        :param paf_dict: Dict of mappings
        :param n_unmapped: Number of unmapped reads
        :param n_reject: Number of rejected reads
        :return:
        """
        # Control: total number of bases + alpha for each read
        self.time_control += total_bases
        self.time_control += (self.batchsize * self.alpha)
        # BR: all accepted bases and ((mu + rho) * number of unmapped/rejected reads) + (alpha * batch_size)
        bases_br = np.sum([r[0].qlen for r in paf_dict.values()])
        self.time_boss += bases_br
        self.time_boss += (n_unmapped * (self.mu + self.rho))
        self.time_boss += (n_reject * (self.mu + self.rho))
        self.time_boss += (self.batchsize * self.alpha)
        logging.info(f"time control: {self.time_control}")
        logging.info(f"time boss-runs: {self.time_boss}")



    def update_times_aeons(self, read_sequences: dict[str, str], reads_decision: dict[str, str]) -> None:
        """
        Increment pseudotime for control and boss regions
        on a flowcell during simulations

        :param read_sequences: Dict of raw sequences, accepting everything
        :param reads_decision: Dict of processed sequences after decisions
        :return:
        """
        # for control: all reads as they come out of the sequencer
        # total bases + (#reads * alpha)
        bases_total = np.sum([len(seq) for seq in read_sequences.values()])
        acquisition = self.batchsize * self.alpha
        self.time_control += (bases_total + acquisition)
        logging.info(f"time control: {self.time_control}")

        # for aeons: bases of the fully sequenced reads (accepted & unmapped) and of the truncated reads
        read_lengths_decision = np.array([len(seq) for seq in reads_decision.values()])
        n_reject = np.sum(np.where(read_lengths_decision == self.mu, 1, 0))
        bases_aeons = np.sum(read_lengths_decision)
        rejection_cost = n_reject * self.rho
        self.time_boss += (bases_aeons + acquisition + rejection_cost)
        logging.info(f"time boss-aeons: {self.time_boss}")



    def fill_cache(self, read_sequences: dict[str, str], reads_decision: dict[str, str]) -> None:
        """
        Write batches of simulated reads for convenient processing after the experiment
        Either only adds reads to a cache, or dumps them to file if it's time

        :param read_sequences: Dict of raw sequences
        :param reads_decision: Dict of sequences after BOSS decisions
        :return:
        """
        # helper function for both conditions
        def add_to_cache(seqs, cache):
            for rid, seq in seqs.items():
                cache[rid] = seq
        # add the current sequences to the cache
        add_to_cache(seqs=read_sequences, cache=self.cache_control)
        add_to_cache(seqs=reads_decision, cache=self.cache_boss)
        # check if time to dump and execute
        self._prep_dump(cond='control')
        self._prep_dump(cond='boss')



    def _prep_dump(self, cond: str) -> None:
        """
        Check if it's time to write cache to file, and execute if it is

        :param cond: String to determine which attributes to grab
        :return:
        """
        # grab the attributes of the condition
        curr_time = getattr(self, f'time_{cond}')
        dump_number = getattr(self, f'dump_n_{cond}')
        cache = getattr(self, f'cache_{cond}')
        # check if it's time to write out the next file
        if curr_time > (self.dumptime * dump_number):
            self._execute_dump(cond=cond, dump_number=dump_number, cache=cache)



    def _execute_dump(self, cond: str, dump_number: int, cache: dict[str, str]) -> None:
        """
        Write out the next cumulative batch file from the cache.
        This is for simulations only

        :param cond: String indicating which region on the flowcell
        :param dump_number: Running number of the dump to execute
        :param cache: Cache of read sequences
        :return:
        """
        logging.info(f'dump {cond} #{dump_number}. # of reads {len(list(cache.keys()))}')
        filename = f'00_reads/{cond}_{dump_number}.fa'
        # copy previous file to make cumulative
        previous_filename = f'00_reads/{cond}_{dump_number - 1}.fa'
        try:
            execute(f"cp {previous_filename} {filename}")
        except FileNotFoundError:
            # at the first batch, create empty 0th and copy to 1st file
            # to make sure we don't append to the same file multiple times
            # otherwise we have duplicate reads causing issues
            empty_file(previous_filename)
            execute(f"cp {previous_filename} {filename}")
        # writing operation
        with open(filename, "a") as f:
            for rid, seq in cache.items():
                r = random_id()
                fa_line = f'>{rid}.{r}\n{seq}\n'
                f.write(fa_line)
        # increment dump counter
        setattr(self, f'dump_n_{cond}', dump_number + 1)
        # reset cache
        setattr(self, f'cache_{cond}', dict())






