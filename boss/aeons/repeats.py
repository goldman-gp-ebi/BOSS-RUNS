import logging

import numpy as np

from boss.mapper import Mapper
from boss.utils import find_blocks_ge
from boss.aeons.sequences import SequencePool




class RepeatFilter:
    """
    A class to:
     - extract repeats
     - generate a fasta
     - map batches to repeat library for filtering

    """

    def __init__(self, name: str, seqpool: SequencePool):
        """
        Initialise a RepeatFilter, incl indexing the sequences, chopping and mapping them

        :param name: Name of experiment
        :param seqpool: Sequence pool of initial data
        """
        self.seqpool = seqpool
        self.name = name
        self.library = f'{name}.repeat_lib.fa'
        # initialise a mapper against the long reads
        SequencePool.write_seq_dict(seqpool.seqdict(), f'{name}.seqs.fa')
        lr_mapper = Mapper(ref=f'{name}.seqs.fa')
        # chop and map all sequences
        little_seqs = self._chop_seqs()
        mappings = lr_mapper._mappy_batch(sequences=little_seqs)
        covs = self._count_cov(mappings)
        self.covs = covs
        # # DEBUG write mappings to file
        # with open("little_seqs_mappings.paf", 'w') as lsm:
        #     lsm.write(mappings)
        # import pickle
        # with open("coverage.pkl", 'wb') as covpkl:
        #     pickle.dump(covs, covpkl)
        # find the min cov for repeats
        self._find_limit()
        # find repeats
        repeat_blocks = self._identify_repeat_sites()
        # write to library file
        self._write_repeat_seqs(repeat_blocks)


    def _chop_seqs(self, window: int = 100, step: int = 100) -> dict:
        """
        Chop the current sequences into smaller bits using a sliding window

        :param window: Size of sliding window
        :param step: Stepsize of sliding windows
        :return: Dictionary of chopped sequences to map
        """
        seqs = self.seqpool.seqdict()
        little_seqs = {}
        for header, seq in seqs.items():
            i = 0
            while i < len(seq):
                little_seqs[f'{header}-{i:010}'] = seq[i: i + window]
                i += step
        with open('little_seqs.fa', 'w') as fh:
            for header, seq in little_seqs.items():
                fh.write(f'>{header}\n')
                fh.write(f'{seq}\n')
        return little_seqs


    @staticmethod
    def _count_cov(mappings: str) -> dict:
        """
        loop through mappings and count the coverage of all targets

        :param mappings: String of a PAF file
        :return: Dict of coverages
        """
        covs = {}
        for line in mappings.split('\n'):
            rec = line.split('\t')
            # check if target array exists
            if not rec[5] in covs.keys():
                covs[rec[5]] = np.zeros(shape=int(rec[6]))
            # grab the cov
            c = covs[rec[5]]
            c[int(rec[7]): int(rec[8])] += 1
        return covs



    def _find_limit(self) -> None:
        """
        Find the limit of coverage to consider repeat sequence
        :return:
        """
        # find the max value first
        maximum = 0
        for c in self.covs.values():
            cmax = np.max(c)
            if cmax > maximum:
                maximum = cmax

        # count all coverage values
        bcounts = np.zeros(int(np.max(maximum) + 1), dtype="int")
        for c in self.covs.values():
            c[0] = 0  # make sure count starts at 0
            bcounts_arr = np.bincount(c.astype('int'))
            for i in range(len(bcounts_arr)):
                bcounts[i] += bcounts_arr[i]

        # limit
        lim = np.quantile(np.repeat(np.arange(len(bcounts)), repeats=bcounts), 0.999)
        if lim < 3:
            lim = 3.0
        self.lim = lim


    def _identify_repeat_sites(self) -> dict:
        """
        find positions where supposed repeats are in sequences
        :return: Dict of repeat blocks for each header
        """

        repeat_blocks = {}
        for header, cov in self.covs.items():
            blocks = find_blocks_ge(cov, self.lim, min_len=100)
            if len(blocks) > 0:
                repeat_blocks[header] = blocks
        return repeat_blocks


    def _write_repeat_seqs(self, repeat_blocks: dict) -> None:
        """
        Write the repeat seqs to file and collect them in dictionary

        :param repeat_blocks: Dict of repeat block coordinates
        :return:
        """
        n_seqs = 0
        repeats = {}
        with open(self.library, 'w') as fh:
            for header, blocks in repeat_blocks.items():
                for start, end in blocks:
                    r = Repeat(header, start, end)
                    r.get_sequence(seqpool=self.seqpool.sequences)
                    fa = r.fasta()
                    fh.write(fa)
                    repeats[r.header] = r.seq
                    n_seqs += 1
        self.repeats = repeats


    @staticmethod
    def _check_coverage(rep_cov: dict, window: int = 500) -> set:
        """
        Check whether a read has a potential repeat on either end

        :param rep_cov: Dict of coverage counts
        :param window: size of end windows
        :return: set of read ids with potential repeats at end
        """
        danger = set()
        for header, rcov in rep_cov.items():
            beginning = rcov[: window]
            if np.sum(beginning) > 5:
                danger.add(header)
            ending = rcov[window:]
            if np.sum(ending) > 5:
                danger.add(header)
        return danger



    def filter_batch(self, seq_dict: dict) -> dict:
        """
        Check a dict of input sequences against a repeat library

        :param seq_dict: Dict of input sequences
        :return: Dict of sequences with potential repeat-seqs removed
        """
        logging.info("repeat filtering batch of reads")
        # write batch to file
        bfile = f'{self.name}.batch.fa'
        with open(bfile, 'w') as fh:
            for header, seq in seq_dict.items():
                fa = f'>{header}\n{seq}\n'
                fh.write(fa)
        # initialise a LinearMapper object
        lm = Mapper(ref=bfile)
        # first map them to the library
        mappings = lm._mappy_batch(self.repeats)
        rep_cov = self._count_cov(mappings)
        danger_ids = self._check_coverage(rep_cov)
        filtered_seqs = {h: s for h, s in seq_dict.items() if h not in danger_ids}
        return filtered_seqs




class Repeat:

    def __init__(self, rid: str = None, start: int = 0, end: int = -1):
        """
        Initialise a repeat object

        :param rid: ID of source sequence
        :param start: Start pos on source sequence
        :param end: End pos on source sequence
        """
        self.rid = rid
        self.start = start
        self.end = end
        self.seq = ''


    def get_sequence(self, seqpool: dict):
        """
        Get sequence from Sequencepool dict using ID and coordinates

        :param seqpool: Dict of SequencePool, i.e. SequencePool.sequences
        :return:
        """
        # index into seqpool and trim
        try:
            self.seq = seqpool[self.rid].seq[self.start: self.end]
        except KeyError:
            logging.info(f'{self.rid} not found in seqpool')
            return


    def fasta(self) -> str:
        """
        Generate fasta representation of itself

        :return: string representation in FASTA format
        """
        if not self.seq:
            return ""
        # construct fasta entry
        self.header = f'{self.rid}-{self.start}:{self.end}'
        fa = f'>{self.header}\n{self.seq}\n'
        return fa






