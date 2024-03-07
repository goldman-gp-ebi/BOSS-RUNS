from itertools import product
from math import sqrt
from typing import Any

import numpy as np
from numpy.typing import NDArray

from boss.utils import reverse_complement



class KmerCounter:

    def __init__(self):
        """
        Initialize the KmerCounter object.
        Specifically for di-, tri-, and tetramers
        General purpose class that can be initialized once, then used for many sequences
        """
        # initialise kmer_tables for n = {2, 3, 4}
        self.t2, self.k2 = self._kmer_table(2)
        self.t3, self.k3 = self._kmer_table(3)
        self.t4, self.k4 = self._kmer_table(4)


    @staticmethod
    def _kmer_table(k: int) -> tuple[NDArray, list[str]]:
        """
        Generate a kmer table.

        :param k: The length of the kmers.
        :return: The kmer table and the list of kmers as strings.
        """
        # initialise nucleotides and their ordinants
        NUC = ['A', 'C', 'G', 'T']
        ordinants = [ord(c) for c in NUC]
        # generate all possible kmers for k
        kmers_str = [''.join(nu) for nu in product(NUC, repeat=k)]
        # the respective integer kmers
        # numbers modified to reduce table size
        kmers = np.array(list(product(ordinants, repeat=k))) - 64
        kmers = np.clip(kmers, 1, 8)
        # create multi-dim array for indexing
        ishape = [np.max(kmers) + 1] * k
        table = np.zeros(shape=ishape, dtype=np.uint8)
        # populate table
        kmer_ind = np.arange(len(kmers))
        for j in range(len(kmers)):
            # using tuple to NOT trigger fancy indexing
            table[tuple(kmers[j])] = kmer_ind[j]
        return table, kmers_str


    @staticmethod
    def _integer_seq(seq: str) -> NDArray:
        """
        Convert a sequence into an integer array.

        :param seq: The input sequence.
        :return: The integer array representing the sequence.
        """
        # prepare a sequence with its reverse complement in byte array
        pseq = np.array([seq + reverse_complement(seq)], dtype=bytes)
        # use view to cast to integer array
        # apply same transformation to reduce table size
        iseq = pseq.view('|S1').view(np.uint8) - 64
        iseq = np.clip(iseq, 1, 8)
        return iseq


    @staticmethod
    def _translate_into_kmer_indices(iseq: NDArray, k: int, table: NDArray) -> NDArray:
        """
        Translate the integer sequence into kmer indices using the kmer table.

        :param iseq: The integer array representing the sequence.
        :param k: The length of the kmers.
        :param table: The kmer table.
        :return: The array of kmer indices.
        """
        # use indexing to get unique numbers for imers
        ituple = tuple(iseq[i: -(k-i+1)] for i in range(k))
        imer_t = table[ituple]
        return imer_t


    def count(self, seq: str, k: int) -> dict[str, int]:
        """
        Count the occurrences of kmers in the given sequence.

        :param seq: The input sequence.
        :param k: The length of the kmers.
        :return: The dictionary of kmer counts.
        """
        # str to integer array
        iseq = self._integer_seq(seq)
        # grab the correct pre-computed table
        table = getattr(self, f't{k}')
        # translate to kmer indices
        tseq = self._translate_into_kmer_indices(iseq, k, table)
        # fast counting with bincount
        counts = np.bincount(tseq)
        # make a dict for results
        count_dict = dict(zip(getattr(self, f'k{k}'), counts))
        return count_dict


    @staticmethod
    def _expected_tetramer_frequencies(km: list[dict[str, int]]) -> dict[str, float]:
        """
        Calculate the expected tetramer frequencies.

        :param km: The list of kmer count dictionaries for di-, tri-, and tetramers.
        :return: The dictionary of expected tetramer frequencies.
        """
        tetra_exp = {}
        tetra = km[2]
        for tet in tetra:
            tetra_exp[tet] = (1.0 * km[1][tet[:3]] * km[1][tet[1:]] / km[0][tet[1:3]])
        return tetra_exp


    @staticmethod
    def _estimate_zscores(km: list[dict[str, int]], tetra_exp: dict[str, float]) -> dict[str, float]:
        """
        Estimate the z-scores for tetramers.

        :param km: list of kmer count dictionaries for di-, tri-, and tetramers.
        :param tetra_exp: dictionary of expected tetramer frequencies.
        :return: dictionary of tetramer z-scores.
        """
        tetra_sd = {}
        for tet, exp in list(tetra_exp.items()):
            den = km[0][tet[1:3]]
            tetra_sd[tet] = sqrt(exp * (den - km[1][tet[:3]]) * (den - km[1][tet[1:]]) / (den * den))

        tetra_z_num = np.array([km[2][tet] - exp for tet, exp in tetra_exp.items()])
        tetra_z_denom = np.array([tetra_sd[tet] for tet in tetra_exp.keys()])
        tetra_z_arr = np.divide(tetra_z_num, tetra_z_denom, where=tetra_z_denom > 0)
        tetra_z = dict(zip(tetra_exp.keys(), tetra_z_arr))
        return tetra_z


    def tetra_zscores(self, seq: str) -> dict[str, float]:
        """
        WRAPPER to calculate the z-scores for tetramers in the given sequence.

        :param seq: The input sequence.
        :return: The dictionary of tetramer z-scores.
        """
        # first count 2,3,4-mers
        kmers = [dict()] * 3
        kmers[0] = self.count(seq, 2)
        kmers[1] = self.count(seq, 3)
        kmers[2] = self.count(seq, 4)
        # calculate expected frequencies
        tetramer_exp = self._expected_tetramer_frequencies(kmers)
        tm_zscores = self._estimate_zscores(km=kmers, tetra_exp=tetramer_exp)
        return tm_zscores



class TetramerDist:

    def __init__(self):
        """
        This object takes sequence objects, not raw sequences.
        To avoid recomputing kmers.
        """
        self.kmc = KmerCounter()


    def euclidean_dist(self, seqo1: Any, seqo2: Any) -> float:
        """
        Calculate the Euclidean distance between two sequence objects.

        :param seqo1: first sequence object.
        :param seqo2: second sequence object.
        :return: Euclidean distance.
        """
        # check if kmers have already been counted
        t1 = getattr(seqo1, 'tmers', None)
        t2 = getattr(seqo2, 'tmers', None)
        # count them if not
        if not t1:
            seqo1.tmers = self.kmc.count(seqo1.seq, 4)
        if not t2:
            seqo2.tmers = self.kmc.count(seqo2.seq, 4)
        # grab the tetramers
        t1 = seqo1.tmers
        t2 = seqo2.tmers
        # intersection to get all 4mers present in both sequences
        tetramers = set(sorted(t1.keys())) & set(sorted(t2.keys()))
        c1 = np.array([t1[t] for t in tetramers])
        c2 = np.array([t2[t] for t in tetramers])
        # normalize
        n1 = c1 / np.sum(c1)
        n2 = c2 / np.sum(c2)
        # difference to calc Euclidean distance
        diff = n1 - n2
        euc = np.sqrt(np.sum(diff * diff))
        return euc




class IntraProb:

    def __init__(self):
        """
        This class calculates the intra-probability and threshold values.
        """
        self.mean = 0
        self.std = 0.01037897 / 2
        self.mean2 = 0.0676654
        self.std2 = 0.03419337
        # self.t = self.calc_threshold()
        self.t = 0.036  # hard-coded to avoid scipy dependency


    # def intra_prob(self, e: NDArray) -> float:
    #     """
    #     Calculate the intra-probability.
    #
    #     :param e: input value
    #     :return: intra-probability.
    #     """
    #     from scipy.stats import norm
    #     self.norm_intra = norm(self.mean, self.std)
    #     self.norm_inter = norm(self.mean2, self.std2)
    #     prob = self.norm_intra.pdf(e) / (self.norm_inter.pdf(e) + self.norm_intra.pdf(e))
    #     return prob
    #
    #
    # def calc_threshold(self) -> float:
    #     """
    #     Calculate the threshold value based on the empirical parameters.
    #
    #     :return: threshold value.
    #     """
    #     e = np.arange(0, 0.1, 0.001)
    #     ip = self.intra_prob(e)
    #     t = e[np.where(ip > 1e-10)][-1]
    #     return t


# instantiate to make available when imported
kmc = KmerCounter()
count_kmers = kmc.count                                 # basic counting function
tetramer_zscores = kmc.tetra_zscores                    # wrapper for tetramer z-scores

tdist = TetramerDist()
euclidean_dist = tdist.euclidean_dist                   # distance function for sequence OBJECTS
# pearson_dist = tdist.pearson_cor                      # distance function for sequence OBJECTS
euclidean_threshold = IntraProb().t                     # constant from empirical parameters











