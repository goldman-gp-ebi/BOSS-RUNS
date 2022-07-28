import logging

# non-std lib
import numpy as np
from scipy.special import betaln

# custom imports
from .BR_batch import CurrentBatch
choose_multimapper = CurrentBatch.choose_multimapper




class MergedGenome:

    def __init__(self, otu, windowSize, alpha, p0):
        self.otu = otu
        self.timeCost = otu.timeCost
        # for fhat expansion
        self.target_size = int(otu.ch_cum_sum[-1] // otu.args.window) + 1
        self.on_target = otu.on_target
        self.readStartCounts = None
        self.fhat = None
        self.windowSize = windowSize
        self.alpha = alpha
        self.p0 = p0
        # initialise read start probs
        self._init_f()


    def _init_f(self):
        '''
        initialise the probability distribution of fragment start sites
        requires prior for alpha (hyperparameter for dirichlet prior) and p0

        Returns
        -------

        '''
        logging.info(f"init read start distribution Fhat")
        # keep track of read start positions (forward and rev) in windows of length windowSize
        readStartCounts = np.zeros(shape=(int(self.otu.ch_cum_sum[-1] / self.windowSize), 2), dtype='float64')
        self.readStartCounts = readStartCounts
        # updating function can be used to initialise
        self.fhat = self.updateFpointmass()
        self.windowSize = self.windowSize


    def countReadStarts(self, cigars_full):
        '''
        Keep track of the read starting positions C_{i,o}
        Used to update F during sequencing
        Read starts are saved in non-overlapping windows of windowsize

        Parameters
        ----------
        cigars_full: dict
            dict of mappings of full length reads

        Returns
        -------

        '''
        # collect all starting positions
        startsForward = []
        startsReverse = []
        windowNumber = self.readStartCounts.shape[0]

        for read_id in cigars_full.keys():
            # get cigars
            cigars = cigars_full[read_id]
            # choose highest ranked mapping
            if len(cigars) > 1:
                cigars = choose_multimapper(cigars)

            for cigar in cigars:
                try:
                    chrom_index = self.otu.chromosome_indices[str(cigar.target_name)]
                except KeyError:
                    continue
                chrom_offset = self.otu.ch_cum_sum[chrom_index]
                cigar_start = chrom_offset + cigar.target_start

                # add a count to list
                if cigar.strand == '+':
                    startsForward.append(cigar_start)
                else:
                    startsReverse.append(cigar_start)

        # count the number of read starts in windows
        binsForward = np.histogram(
            startsForward,
            bins=windowNumber,
            range=(0, self.windowSize * windowNumber))[0].astype(dtype='float')
        binsReverse = np.histogram(
            startsReverse,
            bins=windowNumber,
            range=(0, self.windowSize * windowNumber))[0].astype(dtype='float')

        # add new counts to the array
        self.readStartCounts[:, 0] += binsForward
        self.readStartCounts[:, 1] += binsReverse


    def updateFpointmass(self):
        '''
        Bayesian approach to update posterior values of Fhat
        using counts of read starting positions.
        version with point mass at 0 (for sites with C == 0)

        Returns
        -------
        Fhat: np.array
            probabilities of read starting positions given observed counts

        '''
        windowNumber = self.readStartCounts.shape[0]
        Fhat = np.zeros(shape=self.readStartCounts.shape)
        # First, sites with C > 0
        nonzero_indices = np.nonzero(self.readStartCounts)
        nonzero = self.readStartCounts[nonzero_indices]
        num = np.add(self.alpha, nonzero)
        Csum = np.sum(nonzero)
        denom = 2 * windowNumber * self.alpha + Csum
        Fhat[nonzero_indices] = np.divide(num, denom)
        # then sites with C == 0
        rhs = (self.alpha / (2 * windowNumber * self.alpha + Csum))
        beta_num = np.exp(betaln(self.alpha, ((2 * windowNumber - 1) * self.alpha + Csum)))
        beta_denom = np.exp(betaln(self.alpha, ((2 * windowNumber - 1) * self.alpha))) or 1e-20
        p0_bit = self.p0 / (self.p0 + (1 - self.p0))
        lhs = 1 - p0_bit * (beta_num / beta_denom)
        expectedPost = lhs * rhs
        # mask for the zero count sites - derived from nonzero indices
        zero_indices = np.ones(shape=Fhat.shape, dtype="bool")
        zero_indices[nonzero_indices] = 0
        Fhat[zero_indices] = expectedPost
        return Fhat


    def expandF(self, fhat, downsample_window):
        '''
        Expand and normalise Fhat from the compacted length of genome_length / windowSize
        takes both windowing schemes into account: fhat windows and downsampling windows

        Parameters
        ----------
        fhat: np.array
            probabilities of read starting positions
        downsample_window: int
            length of windows used in generating decision strategies

        Returns
        -------
        fhat_exp: np.array
            read start probs expanded from genome_length / windowSize

        '''
        nrep = int(self.windowSize // downsample_window)
        # expand to the downsampled size
        fhat_exp = np.repeat(fhat, nrep, axis=0)
        # correct for small length difference
        lendiff = self.target_size - fhat_exp.shape[0]
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
            normalizer = self.on_target / fhat_sum
            fhat_exp = np.multiply(fhat_exp, normalizer)
        return fhat_exp


    def estimateFhatPriors(self):
        '''
        Estimate alpha, the concentration hyperparameter of the dirichlet prior,
        by equating the variance of Fhat with the variance of the dirichlet.
        And p0, of the extension to the bayesian approach

        Returns
        -------
        alpha: float
            prior for dirichlet
        p0: float
            proportion of gap sites

        '''
        windowNum = self.readStartCounts.shape[0]
        # filter readStartCounts for positions with 0
        zeroCounts = np.count_nonzero(self.readStartCounts == 0)
        p0 = zeroCounts / (windowNum * 2)
        # get sum of mapped reads
        Csum = np.sum(self.readStartCounts) or 1e-30
        # simple estimator of alpha
        Fhat = np.divide(self.readStartCounts, Csum)
        # variance of Fhat to equate with var of D(alpha)
        Vhat = np.var(Fhat, ddof=0) or 1e-30
        # if Vhat is very small, alpha becomes very large
        lhs = (2 * windowNum - 1) / (Vhat * 8 * (windowNum ** 3))
        rhs = 1 / (2 * windowNum)
        alpha = lhs - rhs
        return alpha, p0



