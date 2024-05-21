import logging
import re
from collections import defaultdict
from itertools import permutations, product
from concurrent.futures import ThreadPoolExecutor as TPexe
from typing import Any

import numpy as np
from numpy.typing import NDArray
from numba import njit

from boss.utils import reverse_complement, binc, adjust_length
from boss.paf import Paf, paf_dict_type



class Priors:

    def __init__(self, ploidy: int = 1):
        """
        Initialise a prior container

        :param ploidy: Only haploid and diploid implemented at the moment
        """
        # set ploidy to either haploid or diploid
        if int(ploidy) == 1:
            self.diploid = False
        elif int(ploidy) == 2:
            self.diploid = True
        else:
            raise ValueError("Given ploidy is not defined")
        # initialise phi and priors
        self.len_b, self.len_g, self.phi = self._generate_phi(diploid=self.diploid)
        self._init_phi_stored()
        self._init_prior()




    @staticmethod
    def _generate_phi(
        diploid: bool = False,
        deletion_error: float = 0.03,
        err_missed_deletion: float = 0.1,
        substitution_error: float = 0.04,
    ) -> tuple[int, int, NDArray]:
        """
        Generate phi. For haploid & diploid. Either with deletions as genotype or ignoring them

        :param diploid: diploid or haploid
        :param deletion_error: probability that a base in the genome is skipped in a read (observed - while in truth it's not)
        :param err_missed_deletion: Probability that an actual deletion is missed in a read. e+
        :param substitution_error: Probability that wrong base is read.
        :return: Tuple of number of bases, genotypes, and phi array
        """
        # two distinctions: haploid/diploid & deletions
        if not diploid:
            if not deletion_error:
                # simplest case: haploid, no deletion genotype
                bases = ['A', 'C', 'G', 'T']
                genotypes = bases
                lenB = len(bases)
                lenG = len(genotypes)
                phi = np.zeros((lenB, lenG))
                for b in range(lenB):
                    for g in range(lenG):
                        if b == g:
                            phi[b][g] = 1 - substitution_error
                        else:
                            phi[b][g] = substitution_error / (lenB - 1)

            else:
                # haploid but with deletions
                bases = ["A", "C", "G", "T", "-"]
                genotypes = bases
                lenB = len(bases)
                lenG = len(genotypes)
                phi = np.zeros((lenB, lenG))
                for b in range(lenB - 1):
                    for g in range(lenG - 1):
                        if b == g:
                            phi[b][g] = 1.0 - (substitution_error + deletion_error)
                        else:
                            phi[b][g] = substitution_error / (lenB - 2)

                        # prob that a base is true when deletion is detected
                        phi[lenB - 1][g] = deletion_error

                    # prob that a deletion is true when base is detected
                    phi[b][lenG - 1] = err_missed_deletion / (lenB - 1)

                # prob that deletion is true when deletion is detected
                phi[lenB - 1][lenG - 1] = 1.0 - err_missed_deletion

        else:
            if not deletion_error:
                # diploid, no deletions
                bases = 'ACGT'
                genotypes = ['AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT']
                lenB = len(bases)
                lenG = len(genotypes)
                phi = np.zeros((lenB, lenG))
                for b in range(lenB):
                    for g in range(lenG):
                        # countBase indicates homozygous, heterozygous or diff to reference
                        countBase = genotypes[g].count(bases[b])
                        if countBase == 2:
                            phi[b][g] = 1 - substitution_error
                        elif countBase == 1:
                            phi[b][g] = ((1 - substitution_error) / 2) + (substitution_error / (2 * (lenB - 1)))
                        elif countBase == 0:
                            phi[b][g] = substitution_error / (lenB - 1)

            else:
                # most complex phi: diploid, with deletions
                bases = 'ACGT-'
                genotypes = ['AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT',
                             'A-', 'C-', 'G-', 'T-', '--']
                lenB = len(bases)
                lenG = len(genotypes)
                phi = np.zeros((lenB, lenG))
                prob_no_error = 1.0 - (substitution_error + deletion_error)
                for b in range(lenB - 1):
                    for g in range(lenG - 5):
                        # the standard base genotypes
                        countBase = genotypes[g].count(bases[b])
                        if countBase == 2:
                            phi[b][g] = prob_no_error  # homozygous
                        elif countBase == 1:
                            phi[b][g] = prob_no_error / 2 + substitution_error / (2 * (lenB - 2))  # heterozygous
                        elif countBase == 0:
                            phi[b][g] = substitution_error / (lenB - 2)  # total difference

                    # treatment of deletion genotypes
                    for g in range(10, 14):
                        countBase = genotypes[g].count(bases[b])
                        if countBase == 1:
                            # heterozygous del if base true
                            phi[b][g] = prob_no_error / 2 + err_missed_deletion / (2 * (lenB - 1))
                        elif countBase == 0:
                            # heterozygous del if base false
                            phi[b][g] = substitution_error / (2 * (lenB - 2)) + err_missed_deletion / (2 * (lenB - 1))

                    # homozygous del if base actually true
                    phi[b][lenG - 1] = err_missed_deletion / (lenB - 1)

                # probs when deletion is truth
                for g in range(lenG):
                    countGap = genotypes[g].count(bases[lenB - 1])
                    if countGap == 2:
                        phi[lenB - 1][g] = 1.0 - err_missed_deletion  # homozygous del
                    elif countGap == 1:
                        phi[lenB - 1][g] = (1.0 - err_missed_deletion) / 2 + deletion_error / 2  # heterozygous del
                    elif countGap == 0:
                        phi[lenB - 1][g] = deletion_error  # no del observed

        return lenB, lenG, phi



    def _init_phi_stored(self) -> None:
        """
        initialise the stored version of phi
        :return:
        """
        phi_stored = np.full((self.len_b, self.len_g, 1000), 1.0)
        for i in range(self.len_b):
            for j in range(self.len_g):
                phi_stored[i, j, :] = self.phi[i, j] ** np.arange(1000)
        self.phi_stored = phi_stored



    def _init_prior(self) -> None:
        """
        Initialize matrix of priors for all possible genotypes
        :return:
        """
        if not self.diploid:
            self.priors = self._haploid_priors()
        else:
            self.priors = self._diploid_priors()
        self.prior_dist = np.array([self.priors[0]])



    @staticmethod
    def _haploid_priors(
        deletion_error: float = 0.03,
        theta: float = 0.01,
        del_subs_ratio: float = 0.4
    ) -> NDArray:
        """
        Initialise priors for haploid genotypes

        :param deletion_error: probability that a base in the genome is skipped in a read (observed - while in truth it's not)
        :param theta: Population diversity
        :param del_subs_ratio: Ratio of deletions and substitutions
        :return:
        """
        if not deletion_error:
            # haploid, no deletions
            bases = ['A', 'C', 'G', 'T']
            genotypes = bases
            len_b = len(bases)
            len_g = len(genotypes)
            priors = np.zeros((len_b, len_g))

            for i in range(len_b):
                for j in range(len_g):
                    if i == j:
                        # prior for observing ref
                        priors[i][j] = 1.0 - theta
                    else:
                        # prior for observing base different to ref
                        priors[i][j] = theta / 3


        else:
            # haploid with deletions
            bases = ["A", "C", "G", "T", "-"]
            genotypes = bases
            len_b = len(bases)
            len_g = len(genotypes)
            priors = np.zeros((len_b - 1, len_g))

            for i in range(len_b - 1):
                for j in range(len_g - 1):
                    if i == j:
                        # prior for observing ref
                        priors[i][j] = 1.0 - (theta * (1.0 + del_subs_ratio))
                    else:
                        # prior for observing base different to ref
                        priors[i][j] = theta / 3

            # prior for observing a deletion
            if del_subs_ratio > 0.0001:
                priors[:, -1] = theta * del_subs_ratio
        return priors



    @staticmethod
    def _diploid_priors(
        deletion_error: float = 0.03,
        theta: float = 0.01,
        del_subs_ratio: float = 0.4
    ) -> NDArray:
        """
        Initialise priors for diploid genotypes

        :param deletion_error: probability that a base in the genome is skipped in a read (observed - while in truth it's not)
        :param theta: Population diversity
        :param del_subs_ratio: Ratio of deletions and substitutions
        :return:
        """
        # The prior proportion of diploid genome sites different from the reference expected to be homozygous
        # depends only on population size
        popsize = 1000
        homo = 0.0
        hetero = 0.0
        aN = np.sum(1.0 / (np.arange(1, popsize + 1)))
        for i in range(popsize):
            homo += (1.0 / ((i + 1) * aN)) * ((i + 1) * float(i + 1) / (popsize ** 2))
            hetero += (1.0 / ((i + 1) * aN)) * 2 * ((popsize - (i + 1)) * float(i + 1) / (popsize ** 2))
        p_homo = homo / (homo + hetero)

        if not deletion_error:
            # diploid, no deletions
            bases = 'ACGT'
            genotypes = ['AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT']
            len_b = len(bases)
            len_g = len(genotypes)
            priors = np.zeros((len_b, len_g))

            for b in range(len_b):
                for g in range(len_g):
                    # homo- or heterozygous observations in g compared to ref b
                    countBase = genotypes[g].count(bases[b])
                    if countBase == 2:
                        priors[b][g] = 1 - theta
                    elif countBase == 1:
                        priors[b][g] = ((1 - p_homo) * theta) / 3
                    elif countBase == 0:
                        priors[b][g] = (p_homo * theta) / 3

        else:
            # diploid with deletions
            bases = 'ACGT-'
            genotypes = ['AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT',
                         'A-', 'C-', 'G-', 'T-', '--']
            len_b = len(bases)
            len_g = len(genotypes)
            priors = np.zeros((len_b - 1, len_g))

            for b in range(len_b - 1):
                for g in range(len_g - 5):
                    # standard base situation
                    countBase = genotypes[g].count(bases[b])
                    if countBase == 2:
                        priors[b][g] = 1 - theta * (1 + del_subs_ratio)  # homozygous
                    elif countBase == 1:
                        priors[b][g] = ((1 - p_homo) * theta) / 3  # heterozygous
                    elif countBase == 0:
                        priors[b][g] = (p_homo * theta) / 3  # total difference

                # treatment of deletion genotypes
                for g in range(10, 14):
                    # prior for heterozygous deletion
                    # should we have a distinction between g2 = b_r and g2 != b_r?
                    priors[b][g] = (1 - p_homo) * del_subs_ratio * theta

                # prior for homozygous deletion
                priors[b][len_g - 1] = p_homo * del_subs_ratio * theta
        return priors




    def uniform_priors(self) -> None:
        """
        change the priors of genotype probabilities to be uniform
        :return:
        """
        len_g = self.priors.shape[1]
        uniform_p = 1 / len_g
        self.priors.fill(uniform_p)
        self.prior_dist = np.array([self.priors[0]])






class Scoring:

    def __init__(self, ploidy: int = 1):
        """
        Initialise a scoring container
        :param ploidy: Ploidy of organism in experiment
        """
        self.priors = Priors(ploidy=ploidy)
        self.n_ref = 4  # n of distinct reference nuc
        self.score0, self.ent0 = self.calc_score(scores=np.zeros(1), pos_posterior=self.priors.prior_dist[0:1])




    def init_score_array(self, max_cov: int = 40, size: int = 4) -> None:
        """
        Function to precompute the posterior and scores
        of different coverage patterns. The pattern generating code uses
        itertools.product() to produce all possible combinations of the ranges
        (all orderings, with repeated elements), then all permutations of the
        combinations are made with itertools.permutations()
        (all possible orderings, no repeats). This makes sure the patterns
        are produced for all possible reference bases.
        Posteriors and scores are then calculated and placed in their locations of
        a multidimensional array for fast indexing.

        :param max_cov: Max size of precomputed array
        :param size: defines sparsity of precomputed array
        :return:
        """
        # create cartesian product
        # if deletions are not accounted for, array will be the same but deletion dimension will not
        # influence the score
        p = np.array(
            list(
                product(range(int(max_cov)),
                        range(int(max_cov / size)),
                        range(int(max_cov / (size + 7))),
                        range(int(max_cov / (size + 9))),
                        range(int(max_cov / (size + 11))))
                ), dtype='uint8'
        )
        # all permutations per element
        cov_patterns = set()
        for i in range(len(p)):
            q = set(list(permutations(p[i])))
            cov_patterns |= q

        # transform generated coverage patterns to array
        cov_patterns = np.array(list(cov_patterns), dtype='uint16')
        n = len(cov_patterns)
        # calculate posterior and scores for coverage patterns
        entropies, scores = self.calc_posterior_and_scores(cov_patterns=cov_patterns)
        # construct lookup table with multidimensional array
        self.score_arr = np.zeros(shape=(max_cov, max_cov, max_cov, max_cov, max_cov, self.n_ref))
        self.entropy_arr = np.zeros(shape=(max_cov, max_cov, max_cov, max_cov, max_cov, self.n_ref))
        ind = np.swapaxes(cov_patterns, 0, 1)
        # assign scores and entropies to their locations
        for i in range(n):
            self.score_arr[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = scores[:, i]
            self.entropy_arr[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = entropies[:, i]




    def update_scores(self, contig) -> tuple[NDArray, NDArray]:
        """
        Update the scores of a contig and return the updated scores
        This is a method of the scoring class to make the large score array persistent without reassigning

        :param contig: Contig object for which to update scores
        :return: Arrays of scores and entropy to reassign
        """
        # grab the correct arrays
        scores = contig.scores
        entropy = contig.entropy
        coverage = contig.coverage
        # this is run after coverage has been updated
        cm = contig.change_mask
        # set deletions to 0 if not counting them
        if contig.len_b == 4:
            coverage[:, 4] = 0
        # indices where array_limit is reached
        maxed_indices = np.where(coverage.sum(axis=1) >= 30)[0]
        cm[maxed_indices] = False
        # indices of change
        cc_pos = np.nonzero(cm)[0]
        # target coverage patterns
        cc = coverage[cc_pos]
        # grab base of changed positions
        ref_bases = contig.seq_int[cc_pos]
        # main operation: indexing in pre-computed values using coverage patterns
        scores[cc_pos] = self.score_arr[cc[:, 0], cc[:, 1], cc[:, 2], cc[:, 3], cc[:, 4], ref_bases]
        # prevent maxed sites from being recalculated
        scores[maxed_indices] = np.finfo(float).tiny
        # check positions where scores were not found in the scoreArray
        # i.e. if they remain 0.0
        missing = np.argwhere(scores == 0.0).flatten()
        nmiss = missing.shape[0]
        # if any positions without scores, calc and add to array
        if nmiss != 0:
            missing_patterns = coverage[missing]
            # calculate new scores and entropies
            miss_entropies, miss_scores = self.calc_posterior_and_scores(cov_patterns=missing_patterns)
            # place into pre-computed array
            bases = contig.seq_int[missing]
            scores[missing] = miss_scores[bases, np.arange(nmiss)]
            entropy[missing] = miss_entropies[bases, np.arange(nmiss)]
            ind = np.swapaxes(missing_patterns, 0, 1)
            # assign the calculated scores and entropies to the large array
            for i in range(ind.shape[1]):
                self.score_arr[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = miss_scores[:, i]
                self.entropy_arr[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = miss_entropies[:, i]
        # grab the entropy values of the changed sites
        entropy[cc_pos] = self.entropy_arr[cc[:, 0], cc[:, 1], cc[:, 2], cc[:, 3], cc[:, 4], ref_bases]
        return scores, entropy




    def calc_posterior_and_scores(self, cov_patterns: NDArray) -> tuple[NDArray, NDArray]:
        """
        Calculate posterior and positional benefits. Operates on any
        number of given coverage patterns and return the values for all
        possible reference bases. Used for initialising the score_arr
        and to calculate missing values.

        :param cov_patterns: Array with different patterns of coverage for which to calc the posterior and scores
        :return: tuple of entropy and scores of input patterns
        """
        n = len(cov_patterns)
        # container for scores with same length as cov_patterns
        scores = np.repeat(self.score0, repeats=n, axis=0)
        pattern_posterior = self.calc_posterior(target_cov=cov_patterns)
        # calculate scores from posteriors
        pattern_score = np.zeros((self.n_ref, n))
        pattern_entropy = np.zeros((self.n_ref, n))
        for i in range(self.n_ref):
            tmpScore, tmpEnt = self.calc_score(scores=scores, pos_posterior=pattern_posterior[i, :, :])
            pattern_score[i] = tmpScore
            pattern_entropy[i] = tmpEnt
        return pattern_entropy, pattern_score



    def calc_posterior(self, target_cov: NDArray) -> NDArray:
        """
        Calculate posterior for site patterns.

        :param target_cov: Taget coverage patterns
        :return: Posteriors of coverage patterns shape: (4, 5, n)
        """
        # prevent overflow when indexing into phi_stored
        target_cov[target_cov > 990] = 990
        len_b = self.priors.phi_stored.shape[0]
        len_g = self.priors.phi_stored.shape[1]
        # how many patterns of coverage are input
        n = target_cov.shape[0]
        # repeat the priors for all input patterns
        coverage_post = np.repeat(self.priors.priors[:, np.newaxis], repeats=n, axis=1)
        # init empty phi, since we will use the stored version
        phis = np.full(n, 1.0)
        # calculate the posteriors for all coverage patterns
        for j in range(len_g):
            if j > 0:
                phis.fill(1.0)
            for i in range(len_b):
                phis *= self.priors.phi_stored[i, j, target_cov[:, i]]
            for h in range(4):
                coverage_post[h, :, j] *= phis
        # sum_b is Z_i(D) (Normalizing factor)
        # for the posterior (the likelihood of the new data at the position).
        for _h in range(4):
            sum_b = np.sum(coverage_post[_h, :, :], axis=1)
            sum_b[sum_b < 1e-300] = 1e-300  # prevent the next line to divide by zero
            coverage_post[_h, :, :] /= sum_b[:, np.newaxis]
        return coverage_post



    def calc_score(self, scores: NDArray, pos_posterior: NDArray) -> tuple[NDArray, NDArray]:
        """
        Uses the current posterior probabilities to find the benefit score of a given position.

        :param scores: Array of scores to update
        :param pos_posterior: posterior probabilites to calculate the scores
        :return: Updated arrays of scores and entropies
        """
        len_s = len(scores)
        # Shannon's entropy of the current posterior distribution.
        logs = np.log(pos_posterior, where=pos_posterior > 0.0)
        logs[logs == np.nan] = 0
        entropy = np.sum(-pos_posterior * logs, axis=1)
        # Expected entropy after a new batch of reads
        new_entropy = np.zeros(len_s)
        # probability of observing a new base of each type at the considered position
        observation_probabilities = np.zeros(len_s)
        # posterior probabilities after reading a certain base i
        new_post = np.zeros((len_s, self.priors.len_g))
        for i in range(self.priors.len_b):  # new read base
            np.multiply(pos_posterior, self.priors.phi[i], out=new_post)
            np.sum(new_post, axis=1, out=observation_probabilities)
            observation_probabilities[observation_probabilities == 0] = 1e-300
            new_post /= observation_probabilities[:, np.newaxis]
            np.log(new_post, where=new_post > 0.0, out=logs)
            logs[logs == np.nan] = 0
            for j in range(self.priors.len_g):  # genotype
                new_entropy -= observation_probabilities * new_post[:, j] * logs[:, j]
        scores[:] = entropy - new_entropy
        return scores, entropy


    @staticmethod
    def merge_benefit(contigs) -> tuple[NDArray, NDArray]:
        """
        :return: Concatenated array of read starting sites across all contigs
        """
        c_benefits = [c.additional_benefit for c in contigs.values()]
        c_smu = [c.smu for c in contigs.values()]
        return np.concatenate(c_benefits), np.concatenate(c_smu)




    @staticmethod
    def find_strat_thread(benefit: NDArray, smu: NDArray, fhat: NDArray, time_cost: int) -> tuple[NDArray, float]:
        """
        Finding approximate decision strategy from the current read benefits

        :param benefit: Additional benefit at each position
        :param smu: S_mu at each position
        :param fhat: Read starting probability at each position
        :param time_cost: Additional time cost of sequencing
        :return: Tuple of strategy and threshold of acceptance
        """
        # take downsampling into consideration
        window = 100
        alpha = 300 // window
        rho = 300 // window
        mu = 400 // window
        tc = time_cost // window
        # group benefit into bins of similar values
        # using binary exponent
        benefit_flat = benefit.flatten('F')
        benefit_nz_ind = np.nonzero(benefit_flat)
        benefit_flat_nz = benefit_flat[benefit_nz_ind]
        # to make binary exponents work, normalise benefit values
        normaliser = np.max(benefit_flat_nz)
        benefit_flat_norm = benefit_flat_nz / normaliser
        mantissa, benefit_exponents = np.frexp(benefit_flat_norm)
        # count how often each exponent is present
        # absolute value because counting positive integers is quicker
        benefit_exponents_pos = np.abs(benefit_exponents)
        # multi-thread counting of exponents
        exponent_arrays = np.array_split(benefit_exponents_pos, 12)
        with TPexe(max_workers=12) as executor:
            exponent_counts = executor.map(np.bincount, exponent_arrays)
        exponent_counts = list(exponent_counts)
        # aggregate results from threads
        # target array needs to have the largest shape of the thread results
        max_exp = np.max([e.shape[0] for e in exponent_counts])
        bincounts = np.zeros(shape=max_exp, dtype='int')
        # sum up results from individual threads
        for exp in exponent_counts:
            bincounts[0:exp.shape[0]] += exp

        exponents_unique = np.nonzero(bincounts)[0]  # filter empty bins
        counts = bincounts[exponents_unique]  # counts of the existing benefit exponents
        # perform weighted bincount in multiple threads
        exponent_arrays = np.array_split(benefit_exponents_pos, 12)
        fhat_arrays = np.array_split(fhat.flatten('F')[benefit_nz_ind], 12)
        arguments = zip(exponent_arrays, fhat_arrays)

        with TPexe(max_workers=12) as executor:
            fgs = executor.map(binc, arguments)
        fgs = list(fgs)
        # aggregates results of threads
        max_fg = np.max([f.shape[0] for f in fgs])
        f_grid = np.zeros(shape=max_fg, dtype='float')
        # aggregate results
        for fg in fgs:
            f_grid[0: fg.shape[0]] += fg

        f_grid = f_grid[exponents_unique]  # filter empty bins
        f_grid_mean = f_grid / counts  # mean fhat for exponent bin
        # use exponents to rebuild benefit values
        benefit_bin = np.power(2.0, -exponents_unique) * normaliser
        # average benefit of strategy in the case that all fragments are rejected
        ubar0 = np.sum(fhat * smu)
        tbar0 = alpha + rho + mu
        # cumsum of the benefit (bins multiplied by how many sites are in the bin)
        # weighted by probability of a read coming from that bin
        cs_u = np.cumsum(benefit_bin * f_grid_mean * counts) + ubar0
        cs_t = np.cumsum(tc * counts * f_grid_mean) + tbar0
        peak = cs_u / cs_t
        strat_size = np.argmax(peak) + 1
        # plt.plot(cs_u)
        # plt.plot(cs_t)
        # plt.plot(peak)
        # plt.show()

        # calculate threshold exponent and where values are geq
        try:
            threshold = benefit_bin[strat_size]
        except IndexError:
            threshold = benefit_bin[-1]

        strat = np.where(benefit >= threshold, True, False)
        return strat, threshold







class CoverageConverter:

    def __init__(self, qt: int = 0):
        """
        Class that handles the conversion of read mappings to coverage counts of nucleotides

        :param qt: Minimum quality threshold to count observations
        """
        # set up translation dict for bases
        transDict = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}
        self.base2int = str.maketrans(transDict)
        # translation dict for qualities (phred encoded)
        transDict_qual = {chr(min(126, q + 33)): f'{q},' for q in range(94)}
        self.qual2int = str.maketrans(transDict_qual)
        # and a translation dict for the cigar ops
        transDict_cig = {'M': '0', 'D': '1', 'I': '2', 'S': '3'}
        self.cig2int = str.maketrans(transDict_cig)
        # compile regex for cigar ops
        self.cigar_regex = re.compile("(\d+)([MIDNSHP=XB])")
        # quality threshold
        self.qt = qt



    def convert_records(self, paf_dict: paf_dict_type, seqs: dict[str, str], quals: dict[str, str]
                        ) -> defaultdict[Any, list]:
        """
        Convert mappings to coverage counts

        :param paf_dict: Dict of mappings
        :param seqs: Dict of read sequences
        :param quals: Dict of read qualities
        :return: Defaultdict of lists of coverage counts per target seq
        """
        # container to collect increments per contig
        increments = defaultdict(list)
        # loop through all reads
        read_ids = list(paf_dict.keys())
        for rid in read_ids:
            # grab paf record
            rec = paf_dict[rid]
            if len(rec) > 1:
                rec = Paf.choose_best_mapper(rec)[0]
            else:
                rec = rec[0]

            # range of coverage of the read
            start = min(rec.tstart, rec.tend)
            end = max(rec.tstart, rec.tend)

            # handle strands
            if rec.rev:   # reverse
                seq = reverse_complement(seqs[rec.qname])
                qual = quals[rec.qname][::-1]
                offcut = rec.qlen - rec.qend
            else:
                seq = seqs[rec.qname]
                qual = quals[rec.qname]
                offcut = rec.qstart

            # get array of counts for each base and indels
            query_arr, qual_arr = self._parse_cigar(
                cigar_string=rec.cigar,
                read_string=seq,
                qual_string=qual,
                offcut=offcut
            )
            # make sure the arrays are the correct length
            assert (end - start) == query_arr.shape[0]
            assert (end - start) == qual_arr.shape[0]
            # eliminate low qual coverage
            addition = np.ones(shape=query_arr.shape[0])
            addition[np.where(qual_arr < self.qt)] = 0
            # collect increments
            increments[rec.tname].append((start, end, query_arr, addition))
        return increments




    def _parse_cigar(self, cigar_string: str, read_string: str, qual_string: str, offcut: int) -> tuple[NDArray, NDArray]:
        """
        Parse cigar string to integer base observations

        :param cigar_string: A CIGAR string
        :param read_string: Read that the CIGAR string describes
        :param qual_string: string of base qualities (phred)
        :param offcut: piece of the read at the start, which did not map. Offset for the query_pos
        :return: tuple of integerised base and quality arrays of query at each pos
        """
        # translate sequence and qualities to integers
        read_integer = read_string.translate(self.base2int)
        int_seq = np.frombuffer(read_integer.encode(), 'u1') - ord('0')
        qual_integer = qual_string.translate(self.qual2int)
        int_qual = np.array(qual_integer.split(',')[:-1], dtype="uint8")
        # split up the cigar
        lengths, ops = self._prep_cigar(cigar_string)
        # loop cigar operators
        query_array, qual_array = self._loop_cigar(
            cigar_lengths=lengths,
            cigar_ops=ops,
            int_seq=int_seq,
            int_qual=int_qual,
            qpos=offcut
        )
        return np.array(query_array), np.array(qual_array)



    def _prep_cigar(self, cigar_string: str) -> tuple[NDArray, NDArray]:
        """
        Get arrays for lengths and codes of cigar ops

        :param cigar_string: CIGAR string of an alignment
        :return: Tuple of arrays for lengths and codes of cigar operations
        """
        # takes a cigar string and returns arrays for lengths and opcodes
        # first split cigar into tuples of length, opcode
        parts = self.cigar_regex.findall(cigar_string)
        # cast into lists
        lengths, ops = zip(*parts)
        # convert to array
        lengths_arr = np.array(lengths, dtype=np.uint32)
        # for the opcodes, first put them into a string and translate to integers
        ops_seq = ''.join(ops)
        ops_tr = ops_seq.translate(self.cig2int)
        opsSeq = np.frombuffer(ops_tr.encode(), 'u1') - ord('0')
        return lengths_arr, opsSeq



    #@njit
    def _loop_cigar(self, cigar_lengths: NDArray, cigar_ops: NDArray, int_seq: NDArray, int_qual: NDArray, qpos: int) -> tuple[list, list]:
        """
        Parse a cigar string given it's deconstructed components

        :param cigar_lengths: Array of length of operations
        :param cigar_ops: Array of cigar operations
        :param int_seq: Interget sequence of read
        :param int_qual: Integer quality of read
        :param qpos: Offcut length
        :return: Tuple of query and quality lists
        """
        query = []
        qual = []
        for count, operation in zip(cigar_lengths, cigar_ops):
            # MATCH
            if operation == 0:
                query.extend(int_seq[qpos: qpos + count])
                qual.extend(int_qual[qpos: qpos + count])
                qpos += count
            # GAP
            elif operation == 1:
                query.extend([4] * count)
                qual.extend([20] * count)
            # INSERTION
            elif operation == 2:
                # query_array[array_pos, 5] += 1  # UNUSED
                # query_array[array_pos, 6] += count # UNUSED
                qpos += count
            # OTHER, e.g. softclip
            elif operation == 3:
                # soft clipped bases are skipped
                # minimap2 does not report softclipped reads anyway
                qpos += count
            else:
                logging.info(f'Invalid cigar operation {operation}')
        return query, qual








