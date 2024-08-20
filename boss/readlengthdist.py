import logging

import numpy as np
from numpy.typing import NDArray


class ReadlengthDist:

    def __init__(self, mu: int = 400, sd: int = 4000, lam: int = 6000, eta: int = 11):
        """
        Initializes a truncated normal distribution of read lengths.
        Used as prior for read length distribution, gets updated throughout the experiment

        :param mu: Length of anchor bases
        :param lam: Lambda parameter of the distribution.
        :param sd: Standard deviation of read lengths.
        :param eta: Number of partitions in the piece-wise approximation.
        """
        self.sd = sd
        self.lam = lam
        self.eta = eta
        self.mu = mu
        self.read_lengths = np.zeros(shape=int(1e6), dtype='uint16')
        # get the maximum read length
        longest_read = int(lam + 10 * sd)
        # prob density of normal distribution
        x = np.arange(longest_read, dtype='int')
        L = np.exp(-((x - lam + 1) ** 2) / (2 * (sd ** 2))) / (sd * np.sqrt(2 * np.pi))
        L /= sum(L)           # normalise
        self.L = L
        # stepwise approx
        self.approx_ccl = self.ccl_approx_constant()



    def update(self, read_lengths: dict) -> None:
        """
        Updates the distribution of read lengths.

        :param read_lengths: Dictionary containing read IDs and their lengths.
        """
        # loop through reads and record their length in array
        for rid, length in read_lengths.items():
            # ignore rejected reads for the read length dist
            # might overestimate the length slightly
            if length > self.mu * 2:
                # reads longer than 1M are counted as 1M
                if length >= 1e6:
                    length = int(1e6) - 1
                self.read_lengths[length] += 1
            else:
                continue

        # calc current stats of read lengths
        observed_read_lengths = np.nonzero(self.read_lengths)
        # if this function was called before any reads are observed, we just return
        if len(observed_read_lengths[0]) == 0:
            logging.info('Attempted update of read lengths before observing any reads')
            return
        length_sum = np.sum(observed_read_lengths * self.read_lengths[observed_read_lengths])
        self.lam = length_sum / np.sum(self.read_lengths[observed_read_lengths])
        self.longest_read = np.max(np.where(self.read_lengths))
        self.L = np.copy(self.read_lengths[:self.longest_read + 1]).astype('float64')
        self.L /= sum(self.L)
        self.approx_ccl = self.ccl_approx_constant()
        logging.info(f'rld: {self.approx_ccl}')
        # lambda - mu - rho
        self.time_cost = self.lam - 400 - 300



    def ccl_approx_constant(self) -> NDArray:
        """
        CCL is 1-CL (cumulative distribution) of read lengths (L).
        Starts at 1 and decreases with increasing length.
        CCL[i] is the probability that a read is at least i+1 long.
        Plus approximated by piece-wise constant function with eta pieces

        :return: approximated pieces of complementary cumulative read length distribution
        """
        # complement of cumulative distribtuion of read lengths
        ccl = np.zeros(len(self.L) + 1)
        ccl[0] = 1
        ccl[1:] = 1 - np.concatenate((self.L[1:].cumsum(), np.ones(1)))
        # cut distribution off to reduce complexity
        ccl[ccl < 1e-6] = 0
        ccl = np.concatenate((np.trim_zeros(ccl, trim='b'), np.zeros(1)))
        self.ccl = ccl
        # approx. with piecewise constant function
        approx_ccl = np.zeros(self.eta - 1, dtype='int32')
        i = 0
        for part in range(self.eta - 1):
            prob = 1 - (part + 0.5) / (self.eta - 1)
            while (ccl[i] > prob) and (len(ccl) > i):
                i += 1
            approx_ccl[part] = i
        return approx_ccl



