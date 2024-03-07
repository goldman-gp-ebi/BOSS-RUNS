from collections import defaultdict

import numpy as np
from numpy.typing import NDArray
from scipy.special import betaln

from boss.paf import Paf



class ReadStartDist:

    def __init__(self, contigs: dict, window_size: int = 2000, alpha: float = 1.0, p0: float = 0.1):
        """
        Initialise the probability distribution of fragment start sites

        :param contigs: Dictionary of contig objects
        :param window_size: Window size of counting read starts
        :param alpha: prior for alpha (hyperparameter for dirichlet prior)
        :param p0: Prior for sites with 0 observations
        """
        self.alpha = alpha
        self.p0 = p0
        self.window_size = window_size
        # track read start positions (forward and rev) in windows
        self.read_starts = {cname : np.zeros(shape=(int(c.length / window_size), 2)) for cname, c in contigs.items()}
        # fhat exists only in its merged form, i.e. for use in updating on a merged array
        self.total_len = np.sum([a.shape[0] for a in self.read_starts.values()])
        self.target_size = int(np.sum([c.length for c in contigs.values()]) // 100)
        self.on_target = 1   # TODO
        self.fhat = self.update_f_pointmass()



    def merge(self) -> NDArray:
        """
        :return: Concatenated array of read starting sites across all contigs
        """
        return np.concatenate(list(self.read_starts.values()))



    def count_read_starts(self, paf_dict: dict[str, list]) -> None:
        """
        Keep track of the read starting positions C_{i,o}, used to update F
        Read starts are saved in non-overlapping windows of windowsize

        :param paf_dict: Dictionary of read mappings
        :return:
        """
        # collect all starting positions
        starts_fwd = defaultdict(list)
        starts_rev = defaultdict(list)

        for rid in paf_dict.keys():
            rec = paf_dict[rid]
            # choose the highest ranked mapping
            if len(rec) > 1:
                rec = Paf.choose_best_mapper(rec)[0]
            else:
                rec = rec[0]

            if rec.rev:
                starts_rev[rec.tname].append(rec.tend)
            else:
                starts_fwd[rec.tname].append(rec.tstart)

        for cname, r_starts in self.read_starts.items():
            # count the number of read starts in windows
            n_windows = int(r_starts.shape[0])
            bins_fwd = np.histogram(
                starts_fwd[cname],
                bins=n_windows,
                range=(0, self.window_size * n_windows))[0].astype(dtype='float')
            bins_rev = np.histogram(
                starts_rev[cname],
                bins=n_windows,
                range=(0, self.window_size * n_windows))[0].astype(dtype='float')

            # add new counts to the array
            self.read_starts[cname][:, 0] += bins_fwd
            self.read_starts[cname][:, 1] += bins_rev



    def update_f_pointmass(self) -> NDArray:
        """
        Bayesian approach to update posterior values of f_hat using counts of read starting positions.
        Version with point mass at 0 (for sites with C == 0)

        :return: Array of read starting probabilities for all positions
        """
        # concatenate arrays
        merged = self.merge()
        n_windows = merged.shape[0]
        fhat = np.zeros(shape=merged.shape)
        # First, sites with C > 0
        nonzero_indices = np.nonzero(merged)
        nonzero = merged[nonzero_indices]
        num = np.add(self.alpha, nonzero)
        Csum = np.sum(nonzero)
        denom = 2 * n_windows * self.alpha + Csum
        fhat[nonzero_indices] = np.divide(num, denom)
        # then sites with C == 0
        rhs = (self.alpha / (2 * n_windows * self.alpha + Csum))
        beta_num = np.exp(betaln(self.alpha, ((2 * n_windows - 1) * self.alpha + Csum)))
        beta_denom = np.exp(betaln(self.alpha, ((2 * n_windows - 1) * self.alpha))) or 1e-20
        p0_bit = self.p0 / (self.p0 + (1 - self.p0))
        lhs = 1 - p0_bit * (beta_num / beta_denom)
        expectedPost = lhs * rhs
        # mask for the zero count sites - derived from nonzero indices
        zero_indices = np.ones(shape=fhat.shape, dtype="bool")
        zero_indices[nonzero_indices] = 0
        fhat[zero_indices] = expectedPost
        # expand from downsampled size
        fhat_exp = self._expand_fhat(fhat)
        return fhat_exp



    def _expand_fhat(self, fhat, downsample_window: int = 100) -> NDArray:
        """
        Expand and normalise Fhat from the compacted length of genome_length / windowSize

        :param fhat: probabilities of read starting positions
        :param downsample_window: length of windows used in generating decision strategies
        :return: read start probs expanded from genome_length / windowSize
        """
        nrep = int(self.window_size // downsample_window)
        # expand to the downsampled size
        fhat_exp = np.repeat(fhat, nrep, axis=0)
        # correct for small length difference
        lendiff = self.target_size - fhat_exp.shape[0]
        assert lendiff < self.window_size
        if lendiff > 0:
            # genome is longer than Fhat
            fhat_exp = np.append(fhat_exp, fhat_exp[-lendiff:], axis=0)
        elif lendiff < 0:
            # genome is shorter than Fhat
            fhat_exp = fhat_exp[:-abs(lendiff)]
        else:
            pass

        # normalise only if not empty
        fhat_sum = np.sum(fhat_exp)
        if fhat_sum != 0:
            # normalise for ratio of on/off target reads
            # uses estimated on-target proportion from initially accepted sites
            normalizer = self.on_target / fhat_sum     # TODO get the on-target estimator
            fhat_exp = np.multiply(fhat_exp, normalizer)
        return fhat_exp



    def estimate_priors(self) -> tuple[float, float]:
        """
        Estimate alpha, the concentration hyperparameter of the dirichlet prior,
        by equating the variance of Fhat with the variance of the dirichlet.

        :return: Prior of dirichlet and proportion of gap sites
        """
        merged = self.merge()
        n_windows = merged.shape[0]
        # filter readStartCounts for positions with 0
        zeroCounts = np.count_nonzero(merged == 0)
        p0 = zeroCounts / (n_windows * 2)
        # get sum of mapped reads
        Csum = np.sum(merged) or 1e-30
        # simple estimator of alpha
        Fhat = np.divide(merged, Csum)
        # variance of Fhat to equate with var of D(alpha)
        Vhat = np.var(Fhat, ddof=0) or 1e-30
        # if Vhat is very small, alpha becomes very large
        lhs = (2 * n_windows - 1) / (Vhat * 8 * (n_windows ** 3))
        rhs = 1 / (2 * n_windows)
        alpha = lhs - rhs
        return alpha, p0



