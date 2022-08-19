import logging
import os
import glob
import time
import sys
from itertools import product, permutations
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor as TPexe

# non-std lib
import numpy as np
import pandas as pd
from natsort import natsorted
from minknow_api.manager import Manager
from minknow_api import __version__ as minknow_api_version


# custom imports
from .BR_utils import execute, read_fa, integer_genome, append_row,\
    grab_br_channels, init_logger, binc, window_sum, adjust_length
from .BR_mapper import Mapper
from .BR_merged_genome import MergedGenome
from .BR_batch import CurrentBatch
from .BR_reference import ReferenceWGS, ReferenceVCF
from .BR_abundance_tracker import AbundanceTracker


class Constants:

    def __init__(self):
        self.eta = 11
        self.windowSize = 2000
        self.alphaPrior = 1
        self.p0 = 0.1
        self.mu = 400
        self.rho = 300
        self.alpha = 300
        self.flank = 10000
        self.base_quality_threshold = 1
        self.substitution_error = 0.06
        self.deletion_substitution_ratio = 0.4
        self.theta = 0.01
        self.deletion_error = 0.05
        self.err_missed_deletion = 0.1
        self.ref_priors = 1
        self.cov_until = 5
        self.window = 100   # windowsize for downsampling (1 == no downsampling)
        self.bucket_size = 20000    # bucket size for strategy switches (large value for single threshold)



class OTU:

    def __init__(self, args, ref, name='otu', ploidy=1):
        self.args = args
        self.name = name
        self.ref = ref
        self.switch = False  # switch whether strategy is on or not
        self.threshold = 0
        self.approx_ccl = 0
        # set ploidy to either haploid or diploid
        if int(ploidy) == 1:
            self.diploid = False
        elif int(ploidy) == 2:
            self.diploid = True
        else:
            logging.info("higher ploidy not defined yet")

        self.linear = True
        self.regions = True
        # for dynamic read length dist, we init an array that keeps track of lengths
        self.read_lengths = np.zeros(shape=int(1e6), dtype='uint16')
        self.frame = []

        # attributes that get overwritten by reference object
        self.chromosome_sequences = {}
        self.chromosome_lengths = {}
        self.chromosome_indices = {}
        self.chromosome_indices_considered = {}
        self.chNum = 0
        self.chArray = np.zeros(1)
        self.ch_cum_sum = np.zeros(1)
        self.genome = np.zeros(1)
        self.roi_length = 0
        self.genome2roi = defaultdict(dict)
        self.genome2roi_arr = {}
        self.chrom_rpos = {}
        self.chrom_rpos_ranges = {}
        self.roi_indices = np.zeros(1)
        # prep rejected refs
        if not self.args.reject_refs:
            self.reject_refs = set()
        else:
            self.reject_refs = set(self.args.reject_refs.split(','))



    def init_phi(self):
        '''
        Initialise \phi, observation probabilities and lookup table,
        given model parameters specified as constants

        Returns
        -------

        '''
        logging.info(f"initialising phi ")

        len_b, len_g, phi = OTU.generate_phi(
            deletion_error=self.args.deletion_error,
            err_missed_deletion=self.args.err_missed_deletion,
            substitution_error=self.args.substitution_error,
            diploid=self.diploid)

        self.len_b = len_b
        self.len_g = len_g
        self.phi = phi
        self.phi_stored = OTU.initPhiStored(phi=self.phi, len_b=self.len_b, len_g=self.len_g)



    @staticmethod
    def generate_phi(deletion_error, err_missed_deletion, substitution_error, diploid):
        """
        Generate phi. For haploid & diploid. Either with deletions as genotype or ignoring them
        Parameters
        ----------
        deletion_error: float
            probability that a base in the genome is skipped in a read (observed - while in truth it's not)
        err_missed_deletion: float
            Probability that an actual deletion is missed in a read. e+
        substitution_error: float
            Probability that wrong base is read.
        diploid: bool
            diploid or haploid

        Returns
        -------
        lenG: int
            Number of genotypes considered
        lenB: int
            Number of bases considered
        phi: np.ndarray
            2D numpy array containing phis
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



    @staticmethod
    def initPhiStored(phi, len_b, len_g):
        '''
        initialise the stored version of phi, up to 1000.

        Parameters
        ----------
        len_g: int
            Number of genotypes considered (A,C,G,T,D)
        len_b: int
            Number of bases considered (A,C,G,T,D)
        phi: np.ndarray
            5x5 2D numpy array containing our phis.

        Returns
        -------
        phi_stored: np.ndarray
            Array that stores phi for different counts of each base.
        '''
        phi_stored = np.full((len_b, len_g, 1000), 1.0)
        for i in range(len_b):
            for j in range(len_g):
                phi_stored[i, j, :] = phi[i, j] ** np.arange(1000)
        return phi_stored



    @staticmethod
    def prior_readlength_dist(mu=400, sd=4000, lam=6000, eta=11):
        '''
        As prior for the read length distribution we use a truncated normal distribution

        Parameters
        ----------
        mu: int
            exclude reads shorter than mu, no strategy needed for those
        sd: float
            sd of normal dist for read lengths
        lam: int
            mean read length
        eta: int
            number of steps in the step-wise approximation to the read length dist

        Returns
        -------
        L: np.array
            array holding distribution of read lengths
        lam: float
            mean read length from the generated truncated normal distribution
        lastl: int
            maximum read length of the distribution

        '''
        logging.info("initialising prior read length distribution")
        # get the maximum read length
        longest_read = int(lam + 10 * sd)
        # prob density of normal distribution
        x = np.arange(longest_read, dtype='int')
        L = np.exp(-((x - lam + 1) ** 2) / (2 * (sd ** 2))) / (sd * np.sqrt(2 * np.pi))
        # exclude reads shorter than mu
        L[:mu] = 0.0
        # normalise
        L /= sum(L)
        # transform lambda to read length from distribution mode
        mean_length = np.average(x, weights=L) + 1
        # get the stepwise approx
        approx_CCL = OTU.CCL_ApproxConstant(L=L, eta=eta)
        return mean_length, longest_read, L, approx_CCL



    @staticmethod
    def CCL_ApproxConstant(L, eta):
        '''
        CCL is 1-CL, with CL the cumulative distribution of read lengths (L).
        CCL starts at 1 and decreases with increased read lengths.
        CCL[i] represents the probability that a random fragment is at least i+1 long.
        eta determines number of pieces for piece-wise linear approx.
        approx_CCL[i]: probability of reads being at least approx_CCL[i] long,
         is approximated by probability 1 - i / (eta-1),
        while the probability of reads being at least approx_CCL[i]+1 long,
        is approximated by probability 1-(i+1)/(eta-1).
        I.e. approx_CCL[i] is the point of the i-th change of value for
        piece-wise constant function.

        Parameters
        ----------
        L: np.array
            array containing the distribution of read lengths
        eta: int
            # of pieces in the piece-wise approximation.

        Returns
        -------
        approx_CCL: np.array
            array of length eta - 1, with each value as size of that
            partition

        '''
        # complement of cumulative distribtuion of read lengths
        CCL = np.zeros(len(L) + 1)
        CCL[0] = 1
        # instead of a loop we use numpy.cumsum, which generates an array itself
        CCL = 1 - L[1:].cumsum()
        # correct numerical errors, that should be 0.0
        CCL[CCL < 1e-10] = 0
        # to approx. U with a piecewise constant function
        # more partitions = more accurate, but slower (see update_U())
        approx_CCL = np.zeros(eta - 1, dtype='int32')
        i = 0
        for part in range(eta - 1):
            prob = 1 - (part + 0.5) / (eta - 1)
            while CCL[i] > prob:
                i += 1
            approx_CCL[part] = i
        return approx_CCL



    def init_timecost(self, mu, rho, lam):
        '''
        Expected increase in time cost for continuing to sequence a read instead of rejecting it.
        Only the increase in cost, not alpha, mu or rho. (lam - mu - rho)
        '''
        tc = lam - mu - rho
        self.timeCost = tc



    def init_coverage(self, modes=False):
        '''
        initialise containers to hold the observed coverage values
        modes are used in simulations to compare conditions

        Parameters
        ----------
        modes: bool
            switch used during simulations to record coverage for other conditions

        Returns
        -------

        '''
        # number of considered sites
        cov_shape = (self.roi_length, 5)

        # container for BR
        self.coverageBR = np.zeros(shape=cov_shape, dtype="uint16")
        # container for other conditions
        if modes:
            self.coverageNaive = np.zeros(shape=cov_shape, dtype="uint16")
            self.coverageRU = np.zeros(shape=cov_shape, dtype="uint16")



    def init_prior(self):
        '''
        Initialize matrix of priors for all possible genotypes
        '''
        logging.info(f"initializing priors")
        # check that theta and del_subs_ratio have reasonable levels
        if self.args.theta * (1.0 + self.args.deletion_substitution_ratio) > 1.0:
            raise ValueError("too much diversity, most of the genome expected to be different from reference")

        # two distinction, like in generate_phi(): ploidy and deletion genotype
        if not self.diploid:
            priors = OTU.haploid_priors(self.args.deletion_error,
                                        self.args.theta,
                                        self.args.deletion_substitution_ratio)

        else:
            priors = OTU.diploid_priors(self.args.deletion_error,
                                         self.args.theta,
                                         self.args.deletion_substitution_ratio)

        # set attr priors
        self.priors = priors
        # init prior distribution; use the first value if not incorp SNP info
        self.prior_dist = np.array([priors[0]])
        # number of possible reference bases (calcPostAndScore)
        self.n_ref = 4



    @staticmethod
    def haploid_priors(deletion_error, theta, del_subs_ratio):
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
    def diploid_priors(deletion_error, theta, del_subs_ratio):
        # DIPLOID
        # Prior proportion of diploid genome sites different from the reference expected to be homozygous
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



    @staticmethod
    def uniform_priors(priors):
        '''
        change the priors of genotype probabilities to be uniform
        useful if ROIs are specifically variant sites

        Parameters
        ----------
        priors: np.array
            initialised priors

        Returns
        -------
        priors: np.array
            overwritten, uniform genotype priors
        prior_dist: np.array
            different shape used for computations

        '''
        len_g = priors.shape[1]
        uniform_p = 1 / len_g
        priors.fill(uniform_p)
        prior_dist = np.array([priors[0]])
        return priors, prior_dist


    def init_reference_vcf(self, vcf, ref=None, mmi=None, min_len=None, reject_refs=None):
        # create the reference object for the VCF scenario
        reference = ReferenceVCF(ref=ref, mmi=mmi, min_len=min_len, reject_refs=reject_refs)
        reference.load_reference(vcf)
        # transfer the attributes to otu instance
        self.__dict__.update(reference.__dict__)


    def init_reference_wgs(self, ref=None, mmi=None, min_len=None, reject_refs=None):
        # get reference object with necessary attributes
        reference = ReferenceWGS(ref=ref, mmi=mmi, min_len=min_len, reject_refs=reject_refs)
        reference.load_reference()
        # transfer all attributes to OTU instance
        self.__dict__.update(reference.__dict__)


    def init_buckets(self, size):
        # put all sites in buckets for strategy switches
        logging.info("initialising strategy switches")
        self.bucket_size = size
        num_roi_in_bucket = np.zeros(shape=int(self.ch_cum_sum[-1] // size) + 1)
        roi_indices = self.roi_indices // size
        # circumvent buffering
        np.add.at(num_roi_in_bucket, roi_indices, 1)
        self.num_roi_in_bucket = num_roi_in_bucket
        # switch for each bucket
        self.bucket_switches = np.zeros(shape=num_roi_in_bucket.shape[0], dtype="bool")


    def init_frames(self):
        self.frame_accept = defaultdict(list)
        self.frame_mean = defaultdict(list)
        self.frame_ent = defaultdict(list)
        self.frame_ent_chrom = defaultdict(list)
        self.frame_cov = defaultdict(list)
        self.frame_corr = defaultdict(list)


    def init_score_array(self, maxCov=40, size=4, model=0):
        '''
        Function to precompute the posterior and scores
        of different coverage patterns. The pattern generating code uses
        itertools.product() to produce all possible combinations of the ranges
        (all orderings, with repeated elements), then all permutations of the
        combinations are made with itertools.permutations()
        (all possible orderings, no repeats). This makes sure the patterns
        are produced for all possible reference bases.
        Posteriors and scores are then calculated and placed in their locations of
        a multidimensional array for fast indexing.

        Parameters
        ----------
        maxCov: int
            define up to which values the scores will be saved in the array.
            Large values can take up huge amounts of memory.
        size: int
            roughly guides how many values are precalculated, i.e. allele frequencies
        model: int
            indicator which model to use. Experimental

        Returns
        -------
        None
        '''
        # create cartesian product
        # if deletions are not accounted for, array will be the same but deletion dimension will not
        # influence the score
        if maxCov >= 20:
            p = list(product(range(int(maxCov)),
                         range(int(maxCov / size)),
                         range(int(maxCov / (size + 7))),
                         range(int(maxCov / (size + 9))),
                         range(int(maxCov / (size + 11)))))

        else:
            p = list(product(range(int(maxCov)),
                             range(1),
                             range(1),
                             range(1),
                             range(1)))

        p = np.array(p, dtype='uint8')
        coveragePatterns = set()

        # all permutations per element
        for i in range(len(p)):
            q = set(list(permutations(p[i])))
            coveragePatterns |= q

        # transform generated coverage patterns to array
        coveragePatterns = np.array(list(coveragePatterns), dtype='uint8')
        n = len(coveragePatterns)
        # calculate posterior and scores for coverage patterns
        entropies, scores = self.calcPostAndScore(coveragePatterns=coveragePatterns, model=model)

        # construct lookup table with multi dimensional array
        self.scoreArray = np.zeros(shape=(maxCov, maxCov, maxCov, maxCov, maxCov, self.n_ref))
        self.entropyArray = np.zeros(shape=(maxCov, maxCov, maxCov, maxCov, maxCov, self.n_ref))

        ind = np.swapaxes(coveragePatterns, 0, 1)

        # assign scores and entropies to their locations
        for i in range(n):
            self.scoreArray[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = scores[:, i]
            self.entropyArray[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = entropies[:, i]



    def init_scores(self, modes=False):
        '''
        Initialise the benefit scores at each site using the priors

        Parameters
        ----------
        modes: bool
            simulations: switch used to initialise additional scores for other conditions

        Returns
        -------

        '''
        logging.info("initialising positional scores")
        # start by calculating score of an empty position
        score0, ent0 = sitewise_score(scores=np.zeros(1),
                                      pos_posterior=self.prior_dist[0:1],
                                      phi=self.phi,
                                      _len_b=self.len_b, _len_g=self.len_g)


        # entropy is only initialised for the rois
        initial_ent = np.zeros(shape=self.roi_length)
        initial_ent.fill(ent0[0])
        entropy_sum = np.sum(initial_ent)

        # score array
        self.initial_scores = np.zeros(shape=self.roi_length)
        self.initial_scores.fill(score0[0])

        # create copies of the scores, entropy & sum of entropy
        self.S_BR = np.copy(self.initial_scores)
        self.sitewise_eBR = np.copy(initial_ent)
        self.eBR = entropy_sum

        # simulations: extra containers
        if modes:
            self.S_Naive = np.copy(self.initial_scores)
            self.sitewise_eNV = np.copy(initial_ent)
            self.eNV = entropy_sum
            self.S_RU = np.copy(self.initial_scores)
            self.sitewise_eRU = np.copy(initial_ent)
            self.eRU = entropy_sum

        # maximum benefit gain achievable
        totalscore = np.nansum(self.initial_scores)
        logging.info(f"total score is: {totalscore}")
        # only first value is required to calculate new scores
        self.initial_scores = self.initial_scores[0:2]



    def init_dummy_strats(self, flank, window, modes=False):
        '''
        Initialise the first strategy given the ROI mask

        Parameters
        ----------
        flank: int
            number of sites surrounding each ROI for read benefit calculation
        window: int
            size of downsampling window
        modes: bool
            SIMULATIONS: initialise strategies for other conditions (control, readfish)

        Returns
        -------

        '''
        logging.info("initialising fist strategy")
        # this will be a separate array for each chromosome
        self.strat_dict_BR = self._init_strategy_dict(window=window, method='BR')

        if self.args.wgs:
            on_target = OTU._find_initial_accept_wgs(strat_dict=self.strat_dict_BR,
                                                     chromosome_indices=self.chromosome_indices,
                                                     reject_refs=self.reject_refs)
        else:
            on_target = OTU._find_initial_accept(strat_dict=self.strat_dict_BR,
                                                 genome2roi=self.genome2roi,
                                                 flank=flank, window=window)
        # save the proportion of on target sites for fhat normalisation
        self.on_target = on_target


        # log number of accepted sites
        for c, c_ind in self.chromosome_indices.items():
            chrom_len = self.strat_dict_BR[c].shape[0]
            f_accept = np.count_nonzero(self.strat_dict_BR[c][:, 0])  # * window
            r_accept = np.count_nonzero(self.strat_dict_BR[c][:, 1])  # * window
            f_perc = f_accept / chrom_len
            r_perc = r_accept / chrom_len
            logging.info(f'initial acceptance rates BR {c}: {f_accept}, {r_accept}; {f_perc}, {r_perc}')

        if modes:
            self.strat_dict_NV = self._init_strategy_dict(window=window, method='NV')
            self.strat_dict_RU = self._init_strategy_dict(window=window, method='RU')

            if self.args.wgs:
                _ = OTU._find_initial_accept_wgs(strat_dict=self.strat_dict_RU,
                                                 chromosome_indices=self.chromosome_indices,
                                                 reject_refs=self.reject_refs)
            else:
                _ = OTU._find_initial_accept(strat_dict=self.strat_dict_RU,
                                             genome2roi=self.genome2roi,
                                             flank=flank, window=window)

        # write to disk at the end
        self.save_strat()



    def _init_strategy_dict(self, window, method):
        chrom_lengths = {c: clen // window for c, clen in self.chromosome_lengths.items()}
        strat_dict = dict()
        for c, clen in chrom_lengths.items():
            chrom_strat = np.zeros(dtype="bool", shape=(clen, 2))
            # set NV to always accept everything
            if method == 'NV':
                chrom_strat[:, :] = 1
            strat_dict[c] = chrom_strat
        return strat_dict



    @staticmethod
    def _find_initial_accept(strat_dict, genome2roi, flank, window):
        # for the initial strategy we set all ROIs and flanks to accept
        flank_size = flank // window
        # for fhat normalisation: estimate for on-target reads
        n_accept = 0
        n_total = 0

        for chrom, position_dict in genome2roi.items():
            try:
                chrom_strat = strat_dict[chrom]
            except KeyError:
                continue

            # adjust positions into their windows
            positions = np.array(list(position_dict.keys()))
            positions = positions // window
            positions = np.unique(positions)  # reduce size of next loop by factor of window

            for pos in positions:
                flank_left = pos - flank_size
                flank_right = pos + flank_size
                if flank_left < 0:
                    flank_left = 0
                if flank_right > chrom_strat.shape[0]:
                    flank_right = chrom_strat.shape[0]
                # forward
                chrom_strat[flank_left: pos, 0] = 1
                # reverse
                chrom_strat[pos: flank_right, 1] = 1

            chrom_strat_nz = np.count_nonzero(chrom_strat[:, 0])
            n_accept += chrom_strat_nz
            n_total += chrom_strat.shape[0]

        # calculate on-target proportion as average of chromosomes
        on_target = n_accept / n_total
        return on_target


    @staticmethod
    def _find_initial_accept_wgs(strat_dict, chromosome_indices, reject_refs):
        # whole genomes: set everything to accept
        for chrom in chromosome_indices.keys():
            try:
                chrom_strat = strat_dict[chrom]
            except KeyError:
                continue
            # forward & reverse
            if chrom not in reject_refs:
                chrom_strat.fill(1)
            else:
                chrom_strat.fill(0)
        # on-target proportion for fhat TODO
        on_target = 1
        return on_target



    def save_strat(self):
        '''
        write newly calculated decision strategies to file
        makes them available to readfish
        also place a marker file that tells readfish to reload the strategies

        Returns
        -------

        '''
        basename = f'{self.args.out_dir}/masks'
        for c, clen in self.chromosome_lengths.items():
            chrom_path = f'{basename}/{c}'
            strats = getattr(self, f'strat_dict_BR')
            np.save(f"{chrom_path}", strats[c])
        # place a marker that the strategies were updated
        execute(f"touch {self.args.out_dir}/masks/masks.updated")


    def add_coverage(self, covADD_naive=0, covADD_RU=0, covADD_BR=0):
        # after calculating scores of a new batch, add up coverage counts
        self.coverageBR = np.add(self.coverageBR, covADD_BR)
        if type(covADD_naive) != int:
            self.coverageNaive = np.add(self.coverageNaive, covADD_naive)
            self.coverageRU = np.add(self.coverageRU, covADD_RU)


    def calcPostAndScore(self, coveragePatterns, model=0):
        '''
        wrapper function to calculate the posterior (from scratch version,
        instead of updating) and the expected benefit. Both with their
        vectorized numpy implementation. This function operates on any
        number of given coverage patterns and return the values for all
        possible reference bases. It is used for initialising the scoreArray
        and to calculate values missing in the scoreArray while updating.

        Parameters
        ----------

        coveragePatterns: np.ndarray
            Array with different patterns of coverage for which to calc the posterior and scores
        model: int
            which posterior model to use

        Returns
        -------
        patternPost: np.array
            array of dimension (n, 4, 5), holds posteriors for all given
            coverage patterns, for all possible reference bases.
        patternScore: np.array
            similar array of dimension (n, 4), score for any pattern, for
            any reference base.

        '''
        n = len(coveragePatterns)
        # container for scores with same length as coveragePatterns
        initialScore = np.repeat(self.initial_scores[0], repeats=n, axis=0)
        n_ref = self.n_ref

        if model == 0:
            patternPost = calcCoveragePost(
                prior=self.priors,
                phi_stored=self.phi_stored,
                targetCoverage=coveragePatterns)

            # calculate scores from posteriors
            patternScore = np.zeros((n_ref, n))
            patternEntropy = np.zeros((n_ref, n))
            for i in range(n_ref):
                tmpScore, tmpEnt = sitewise_score(
                    scores=initialScore,
                    pos_posterior=patternPost[i, :, :],
                    _len_g=self.len_g,
                    _len_b=self.len_b,
                    phi=self.phi)
                patternScore[i] = tmpScore
                patternEntropy[i] = tmpEnt
            return patternEntropy, patternScore

        # elif model == 1:
        #     patternPost = calcCoveragePost_d2(prior=self.priors, phi=self.phi, targetCoverage=coveragePatterns)
        #
        #     # calculate the scores from the posteriors
        #     patternScore = np.zeros((1, n))
        #     patternEntropy = np.zeros((1, n))
        #
        #     tmpScore, tmpEnt = sitewise_score2(pos_posterior=patternPost[0, :, :], phi=self.phi)
        #     patternScore[0] = tmpScore
        #     patternEntropy[0] = tmpEnt
        #     return patternEntropy, patternScore

        else:
            logging.info("invalid model specified")
            exit()


    def update_S_mu_thread(self, scores, mu, window=1):
        '''
        Parallelised version to calculate positional score of initial mu bases

        Parameters
        ----------
        scores: np.array
            benefit scores for each site
        mu: int
            length of initial fragment used for mapping
        window: int
            downsampling window size

        Returns
        -------
        smu: np.array
            benefit score for mu bases from each site
        scores_placed: np.array
            expanded version of benefit scores

        '''
        # expand benefit scores from downsampled version
        scores_placed = np.zeros(shape=int(self.ch_cum_sum[-1] // window) + 1)
        # position of sites in downsampled array
        roi_indices = self.roi_indices // window
        # avoid buffering
        np.add.at(scores_placed, roi_indices, scores)
        # adjust for windowsize
        mu_ds = mu // window
        # prepare arg list for multithreading
        chrom_starts = (self.ch_cum_sum // window)[:-1]
        chrom_ends = (self.ch_cum_sum // window)[1:]

        arglist = []
        for i in range(self.chNum):
            # scores of one chromosome go to one worker
            s = scores_placed[chrom_starts[i]: chrom_ends[i]]
            # arguments unpacked inside function
            arglist.append((s, mu_ds))
            arglist.append((s[::-1], mu_ds))  # reverse

        with TPexe(max_workers=20) as executor:
            smus = executor.map(calc_smu, arglist)

        # aggregate worker results
        smu = np.zeros((scores_placed.shape[0], 2), dtype="float64")
        chrom = 0
        smu_result_list = list(smus)
        for i in range(len(smu_result_list)):
            smu_result, ch_st = smu_result_list[i]
            # place results either in forward or reverse dim
            if i % 2 == 0:
                smu[int(chrom_starts[chrom] + 1): chrom_ends[chrom], 0] = smu_result[3::2]
                smu[chrom_starts[chrom], 0] = smu_result[1]
            else:
                smu_result_rev = smu_result[::-1]
                smu[chrom_starts[chrom]: chrom_ends[chrom], 1] = smu_result_rev[1::2]
                chrom += 1
        return smu, scores_placed


    def update_U_thread(self, scores_placed, eta, window=1):
        '''
        Parallelised version to calculate read benefit

        Parameters
        ----------
        scores_placed: np.array
            benefit scores expanded from downsampling
        eta: int
            number of segments used in piece-wise approx. function
        window: int
            size of downsampling window

        Returns
        -------
        expected_benefit: np.array
            read benefit for each position

        '''
        # downsample read length distribution
        approx_ccl_ds = self.approx_ccl // window
        # prep arg list for multithreading
        chrom_starts = (self.ch_cum_sum // window)[:-1]
        chrom_ends = (self.ch_cum_sum // window)[1:]

        arglist = []
        for i in range(self.chNum):
            # scores for one chromosome go to one worker
            s = scores_placed[chrom_starts[i]: chrom_ends[i]]
            # args unpacked inside function
            targs = (update_U_forward, s, approx_ccl_ds, eta)
            arglist.append(targs)
            targs = (update_U_reverse, s, approx_ccl_ds, eta)
            arglist.append(targs)

        with TPexe(max_workers=20) as executor:
            benefits = executor.map(calc_u, arglist)

        # aggregate worker results
        expected_benefit = np.zeros((scores_placed.shape[0], 2), dtype="float64")
        chrom = 0
        benefits = list(benefits)
        for i in range(len(benefits)):
            # place results either in forward or reverse dim
            if i % 2 == 0:
                expected_benefit[chrom_starts[chrom]: chrom_ends[chrom], 0] = benefits[i]
            else:
                expected_benefit[chrom_starts[chrom]: chrom_ends[chrom], 1] = benefits[i]
                chrom += 1
        return expected_benefit


    def update_length_dist(self, batch, cigars_full):
        '''
        Record the lengths of the newly observed sequencing reads
        Periodically also update the read length distribution

        Parameters
        ----------
        batch: int
            update number
        cigars_full: dict
            dictionary of read mappings. Read_id: list of mappings

        Returns
        -------

        '''
        # in case there was no read in the batch
        if len(cigars_full) == 0:
            return

        # loop through reads to record lengths
        for read_id, record_list in cigars_full.items():
            read_len = 0
            for r in record_list:
                read_len = r.query_len
                break
            # increment this length in recording array
            # ignore rejected reads for read length dist
            # might overestimate the dist slightly
            if read_len > self.args.mu * 2:
                self.read_lengths[read_len] += 1
            else:
                continue

        # periodically calc current stats of read lengths
        if (batch) % 10 == 0:
            logging.info("updating read length distribution")
            observed_read_lengths = np.nonzero(self.read_lengths)
            length_sum = np.sum(observed_read_lengths * self.read_lengths[observed_read_lengths])
            self.lam = length_sum / np.sum(self.read_lengths[observed_read_lengths])
            self.longest_read = np.max(np.where(self.read_lengths))
            self.L = np.copy(self.read_lengths[:self.longest_read]).astype('float64')
            self.L /= sum(self.L)
            # update approx CCL
            self.approx_ccl = OTU.CCL_ApproxConstant(L=self.L, eta=self.args.eta)
            self.approx_CCL = self.approx_ccl  # fix case error
            logging.info("new approx ccl")
            logging.info(self.approx_ccl)
            # update time cost
            self.timeCost = self.lam - self.args.mu - self.args.rho


    def save_metrics_chroms(self, batch, timeNV, timeRU, timeBR, yield_NV, yield_RU, yield_BR):
        # SIMULATIONS
        # bases per second * 60s * num of channels (512 / 3),
        time_factor = 450 * 60 * 170

        # MEAN COVERAGE
        for c, chrom_rois in self.chrom_rpos_ranges.items():
            croi_start, croi_end = chrom_rois

            cov_sum_NV = self.coverageNaive[croi_start: croi_end].sum(axis=1)
            cov_sum_RU = self.coverageRU[croi_start: croi_end].sum(axis=1)
            cov_sum_BR = self.coverageBR[croi_start: croi_end].sum(axis=1)
            cov_mean_NV = cov_sum_NV.mean()
            cov_mean_RU = cov_sum_RU.mean()
            cov_mean_BR = cov_sum_BR.mean()

            logging.info(f'mean cov {c} BR: {cov_mean_BR}')

            # fraction of sites with cov < X
            x = 5
            low_cov_NV = np.where(cov_sum_NV < x)[0].shape[0] / cov_sum_NV.shape[0]
            low_cov_RU = np.where(cov_sum_RU < x)[0].shape[0] / cov_sum_RU.shape[0]
            low_cov_BR = np.where(cov_sum_BR < x)[0].shape[0] / cov_sum_BR.shape[0]

            self.frame_mean = append_row(self.frame_mean,
                                         {'mode': 'Control', 'batch': batch,
                                          'value': cov_mean_NV, 'lowcov': low_cov_NV,
                                          'seq_time': timeNV / time_factor, 'chrom': c, 'yield': yield_NV})

            self.frame_mean = append_row(self.frame_mean,
                                         {'mode': 'Readfish', 'batch': batch,
                                          'value': cov_mean_RU, 'lowcov': low_cov_RU,
                                          'seq_time': timeRU / time_factor, 'chrom': c, 'yield': yield_RU})

            self.frame_mean = append_row(self.frame_mean,
                                         {'mode': 'BOSS-RUNS', 'batch': batch,
                                          'value': cov_mean_BR, 'lowcov': low_cov_BR,
                                          'seq_time': timeBR / time_factor, 'chrom': c, 'yield': yield_BR})

        # write to file
        mean_df = pd.DataFrame(self.frame_mean)
        mean_csv = f"{self.args.out_dir}/metrics/{self.args.run_name}_mean.csv"
        with open(mean_csv, 'w'):
            pass
        mean_df.to_csv(mean_csv)


        # ENTROPY PER CHROMOSOME
        for c, chrom_rois in self.chrom_rpos_ranges.items():
            croi_start, croi_end = chrom_rois

            eNV_sum = self.sitewise_eNV[croi_start: croi_end].sum()
            eRU_sum = self.sitewise_eRU[croi_start: croi_end].sum()
            eBR_sum = self.sitewise_eBR[croi_start: croi_end].sum()

            eNV_mean = self.sitewise_eNV[croi_start: croi_end].mean()
            eRU_mean = self.sitewise_eRU[croi_start: croi_end].mean()
            eBR_mean = self.sitewise_eBR[croi_start: croi_end].mean()


            self.frame_ent_chrom = append_row(self.frame_ent_chrom,
                                              {'mode': 'Control', 'batch': batch,
                                               'value': eNV_sum, 'value_mean': eNV_mean,
                                               'chrom': c,
                                               'seq_time': timeNV / time_factor, 'yield': yield_NV})


            self.frame_ent_chrom = append_row(self.frame_ent_chrom,
                                              {'mode': 'Readfish', 'batch': batch,
                                               'value': eRU_sum, 'value_mean': eRU_mean,
                                               'chrom': c,
                                               'seq_time': timeRU / time_factor, 'yield': yield_RU})

            self.frame_ent_chrom = append_row(self.frame_ent_chrom,
                                              {'mode': 'BOSS-RUNS', 'batch': batch,
                                               'value': eBR_sum, 'value_mean': eBR_mean,
                                               'chrom': c,
                                               'seq_time': timeBR / time_factor, 'yield': yield_BR})

        ent_chrom_df = pd.DataFrame(self.frame_ent_chrom)
        ent_csv = f"{self.args.out_dir}/metrics/{self.args.run_name}_ent_chrom.csv"
        with open(ent_csv, 'w'):
            pass
        ent_chrom_df.to_csv(ent_csv)


        # ENTROPY TOTAL
        eNV_mean = self.sitewise_eNV.mean()
        eRU_mean = self.sitewise_eRU.mean()
        eBR_mean = self.sitewise_eBR.mean()

        self.frame_ent = append_row(self.frame_ent,
                                    {'mode': 'Control', 'batch': batch,
                                     'value': self.eNV, 'value_mean': eNV_mean,
                                     'seq_time': timeNV / time_factor, 'yield': yield_NV})

        self.frame_ent = append_row(self.frame_ent,
                                    {'mode': 'Readfish', 'batch': batch,
                                     'value': self.eRU, 'value_mean': eRU_mean,
                                     'seq_time': timeRU / time_factor, 'yield': yield_RU})

        self.frame_ent = append_row(self.frame_ent,
                                    {'mode': 'BOSS-RUNS', 'batch': batch,
                                     'value': self.eBR, 'value_mean': eBR_mean,
                                     'seq_time': timeBR / time_factor, 'yield': yield_BR})

        ent_df = pd.DataFrame(self.frame_ent)
        ent_csv = f"{self.args.out_dir}/metrics/{self.args.run_name}_ent.csv"
        with open(ent_csv, 'w'):
            pass
        ent_df.to_csv(ent_csv)

        # ACCEPTED SITES
        for c, clen in self.chromosome_lengths.items():
            forward_accept = np.count_nonzero(self.strat_dict_BR[c][:, 0])
            reverse_accept = np.count_nonzero(self.strat_dict_BR[c][:, 1])
            forward_perc = forward_accept / self.strat_dict_BR[c].shape[0]
            reverse_perc = reverse_accept / self.strat_dict_BR[c].shape[0]
            self.frame_accept = append_row(self.frame_accept,
                                           {'mode': 'BOSS-RUNS', 'batch': batch,
                                            'chrom': c, 'orient': '+', 'value': forward_perc,
                                            'seq_time': timeBR / time_factor, 'yield': yield_BR})

            self.frame_accept = append_row(self.frame_accept,
                                           {'mode': 'BOSS-RUNS', 'batch': batch,
                                            'chrom': c, 'orient': '-', 'value': reverse_perc,
                                            'seq_time': timeBR / time_factor, 'yield': yield_BR})

        accept_df = pd.DataFrame(self.frame_accept)
        acc_csv = f"{self.args.out_dir}/metrics/{self.args.run_name}_acc.csv"
        with open(acc_csv, 'w'):
            pass
        accept_df.to_csv(acc_csv)



    def save_metrics_acc(self, batch):
        # PLAYBACK: record accepted sites
        for c, clen in self.chromosome_lengths.items():
            forward_accept = np.count_nonzero(self.strat_dict_BR[c][:, 0])
            reverse_accept = np.count_nonzero(self.strat_dict_BR[c][:, 1])
            forward_perc = forward_accept / self.strat_dict_BR[c].shape[0]
            reverse_perc = reverse_accept / self.strat_dict_BR[c].shape[0]
            self.frame_accept = append_row(self.frame_accept,
                                           {'mode': 'BOSS-RUNS', 'batch': batch,
                                            'chrom': c, 'orient': '+', 'value': forward_perc})

            self.frame_accept = append_row(self.frame_accept,
                                           {'mode': 'BOSS-RUNS', 'batch': batch,
                                            'chrom': c, 'orient': '-', 'value': reverse_perc})

        accept_df = pd.DataFrame(self.frame_accept)
        acc_csv = f"{self.args.out_dir}/metrics/{self.args.run_name}_acc.csv"
        with open(acc_csv, 'w'):
            pass
        accept_df.to_csv(acc_csv)




class BossRun:

    def __init__(self, args):
        const = Constants()
        for c, cval in const.__dict__.items():
            setattr(args, c, cval)

        self.args = args
        self.batch = 0
        # to keep track of off and on target fragments for fhat
        self.num_mapped = 0
        self.num_unmapped = 0
        # total time (in sequencing units)
        self.timeNaive = 0
        self.timeRU = 0
        self.timeBR = 0
        # track yield of each condition
        self.yield_NV = 0
        self.yield_RU = 0
        self.yield_BR = 0
        # store filenames of processed files
        self.processed_files = []
        # for tracking benefits
        self.threshold = 0
        # create directories for logs etc.
        args.out_dir = f'./bossruns_{args.run_name}'
        if not os.path.exists(args.out_dir):
            os.mkdir(args.out_dir)
            os.mkdir(f'{args.out_dir}/masks')
            os.mkdir(f'{args.out_dir}/ckp')
            os.mkdir(f'{args.out_dir}/fq')
            os.mkdir(f'{args.out_dir}/fq/RU')
            os.mkdir(f'{args.out_dir}/fq/BR')
            os.mkdir(f'{args.out_dir}/fq/NV')
            os.mkdir(f'{args.out_dir}/metrics')
        # initialise logger
        init_logger(logfile=f'{args.out_dir}/{args.run_name}.bossruns.log')
        # if multiple conditions on flowcell, grab channels assigned to BR
        if args.channels is not None:
            self.channels = grab_br_channels(channels_toml=args.channels, run_name=args.run_name)
        else:
            self.channels = None
        # initialise the abundance tracker
        self.abundance_tracker = AbundanceTracker()
        # check which scenario
        args.wgs = 1  # default
        if args.vcf:
            args.wgs = 0
        # print args for debugging
        logging.info("BOSS-RUNS")
        for a, aval in self.args.__dict__.items():
            logging.info(f'{a} {aval}')
        logging.info('\n')





    def initialise_OTUs(self, modes=False):
        otu = OTU(args=self.args,
                  ref=self.args.ref,
                  ploidy=self.args.ploidy)
        self.otu = otu

        logging.info(f"starting init ------- ")
        # initialise the reference
        # tmp
        if self.args.wgs:
            otu.init_reference_wgs(ref=self.args.ref, mmi=self.args.ref_idx, reject_refs=self.otu.reject_refs)
        # this is for the case of VCF input
        elif self.args.vcf:
            otu.init_reference_vcf(ref=self.args.ref, mmi=self.args.ref_idx, vcf=self.args.vcf, reject_refs=self.otu.reject_refs)
        else:
            logging.info("Need either of arguments: --wgs --vcf")
            sys.exit()
        # initialise buckets and switches for the strategy
        otu.init_buckets(size=self.args.bucket_size)
        # init frames for saving metrics, only used for sims
        otu.init_frames()
        # Generate phi and phi_stored
        otu.init_phi()
        otu.init_prior()
        # alternative uniform priors
        if not self.args.wgs:
            otu.priors, otu.prior_dist = OTU.uniform_priors(priors=otu.priors)
        # initialise a prior for read length distribution
        otu.lam, otu.longest_read, otu.L, otu.approx_ccl = OTU.prior_readlength_dist()
        # init time cost of acceptance
        otu.init_timecost(mu=self.args.mu, rho=self.args.rho, lam=otu.lam)
        # initialise scores and pre-computed array for updating
        otu.init_scores(modes=modes)
        otu.init_score_array()
        otu.init_coverage(modes=modes)
        # initialise first strategies
        otu.init_dummy_strats(flank=self.args.flank, window=self.args.window, modes=modes)
        logging.info(f"finished init ------- \n")



    def initialise_merged(self, mappy=True):
        '''
        initialise the Mapper class for mapping reads with mappy
        this keeps a persistent index throughout the run
        can be deactivated for sims

        Parameters
        ----------
        mappy: bool
            whether to initialise the mappy object

        Returns
        -------

        '''
        if mappy:
            ref = self.args.ref_idx
            logging.info(f"loading mapping index: {ref}")
            self.mapper = Mapper(mu=self.args.mu)
            self.mapper.init_mappy(ref=ref)

        # initialise read start position probabilities
        self.merged_genome = MergedGenome(otu=self.otu, windowSize=self.args.windowSize,
                                          alpha=self.args.alphaPrior, p0=self.args.p0)




    def read_and_map(self, fastq_file=None, workers=8):
        '''
        SIMULATIONS: ingest and map new sequencing reads
        can use either individual fastq files or
        randomly sample from a single large file

        Parameters
        ----------
        fastq_file: str
            path to a fastq file, if individual files are used
        workers: int
            number of threads used for mapping reads

        Returns
        -------
        paf_full: str
            alignments of full length reads in PAF format
        paf_trunc: str
            alignments of truncated reads in PAF format (for decision making)

        '''
        self.current_batch = CurrentBatch(num=self.batch,
                                          qt=self.args.base_quality_threshold,
                                          window=self.args.window)

        if fastq_file is None:
            # get batch of randomly sampled reads
            self.current_batch.read_batch_parallel(fq=self.args.fq, bsize=self.args.batch_size)
        else:
            # use fq file and save file name
            self.current_batch.read_batch(fastq_file=fastq_file, channels=self.channels)
            self.fq_name = fastq_file.split("/")[-1]

        # map reads in full length and truncated
        paf_trunc = self.mapper.mappy_batch(
            sequences=self.current_batch.read_sequences, truncate=True, workers=workers)
        paf_full = self.mapper.mappy_batch(
            sequences=self.current_batch.read_sequences, truncate=False, workers=workers)
        return paf_full, paf_trunc




    def read_and_inject(self, sampler):
        '''
        SIMULATIONS: alternative to loading new sequencing reads
        here we sample reads and load their pre-computed alignments

        Parameters
        ----------
        sampler: obj
            data sampler object that handles sequences and their alignments

        Returns
        -------
        paf_full: str
            alignments of full length reads in PAF format
        paf_trunc: str
            alignments of truncated reads in PAF format (for decision making)

        '''
        # initialise current batch
        self.current_batch = CurrentBatch(num=self.batch,
                                          qt=self.args.base_quality_threshold,
                                          window=self.args.window)
        # sample batch of reads
        read_lengths, read_sequences, read_qualities, basesTOTAL = sampler.parallel_batch(batch=self.batch)
        self.current_batch.read_ids = set(read_sequences.keys())
        self.current_batch.read_lengths = read_lengths
        self.current_batch.read_sequences = read_sequences
        self.current_batch.read_qualities = read_qualities
        self.current_batch.basesTOTAL = basesTOTAL
        self.sampled_n = len(self.current_batch.read_ids)
        # grab mappings
        paf_full, paf_trunc = sampler.grab_mappings(read_ids=list(self.current_batch.read_lengths.keys()))
        return paf_full, paf_trunc




    def update_score(self, covADD_naive, covADD_RU, covADD_BR):
        '''
        SIMULATIONS: wrapper to update benefit scores of all conditions

        Parameters
        ----------
        covADD_naive: np.array
            newly observed coverage of control
        covADD_RU: np.array
            newly observed coverage of readfish
        covADD_BR: np.array
            newly observed coverage of BR

        Returns
        -------

        '''
        logging.info("updating scores")
        otu = self.otu
        # NAIVE
        otu.S_Naive, otu.sitewise_eNV = update_S_scoreArray(
            otu=otu,
            scores=otu.S_Naive,
            entropy=otu.sitewise_eNV,
            coverage=otu.coverageNaive,
            covADD=covADD_naive)
        # RU
        otu.S_RU, otu.sitewise_eRU = update_S_scoreArray(
            otu=otu,
            scores=otu.S_RU,
            entropy=otu.sitewise_eRU,
            coverage=otu.coverageRU,
            covADD=covADD_RU)
        # BR
        otu.S_BR, otu.sitewise_eBR = update_S_scoreArray(
            otu=otu,
            scores=otu.S_BR,
            entropy=otu.sitewise_eBR,
            coverage=otu.coverageBR,
            covADD=covADD_BR)
        # adjust scores for large deletions and missing regions
        otu.S_Naive, otu.sitewise_eNV, otu.eNV = modify_scores(
            cov=otu.coverageNaive,
            scores=otu.S_Naive,
            ent=otu.sitewise_eNV,
            chrom_rpos=otu.chrom_rpos_ranges)

        otu.S_RU, otu.sitewise_eRU, otu.eRU = modify_scores(
            cov=otu.coverageRU,
            scores=otu.S_RU,
            ent=otu.sitewise_eRU,
            chrom_rpos=otu.chrom_rpos_ranges)

        otu.S_BR, otu.sitewise_eBR, otu.eBR = modify_scores(
            cov=otu.coverageBR,
            scores=otu.S_BR,
            ent=otu.sitewise_eBR,
            chrom_rpos=otu.chrom_rpos_ranges)

        # add new coverage counts to previously observed coverage
        otu.add_coverage(covADD_naive=covADD_naive,
                         covADD_RU=covADD_RU,
                         covADD_BR=covADD_BR)




    def update_times(self, basesNV, basesRU, basesBR, unmapped_RU, unmapped_BR, reject_RU, reject_BR):
        '''
        SIMULATIONS: after each batch of sequencing reads
        add the time taken for sequencing by each condition

        Parameters
        ----------
        basesNV: int
            total number of observed bases for control
        basesRU: int
            total number of observed bases for readfish
        basesBR:
            total number of observed bases for BOSS-RUNS
        unmapped_RU: int
            number of unmapped reads for readfish
        unmapped_BR: int
            number of unmapped reads for bossruns
        reject_RU: int
            number of reads rejected by readfish
        reject_BR: int
            number of reads rejected by BOSS-RUNS

        Returns
        -------

        '''
        # NAIVE: total number of bases + alpha for each read
        timeToNaive = basesNV
        timeToNaive += (self.args.batch_size * self.args.alpha)
        # RU: all accepted bases and ((mu + rho) * number of unmapped/rejected reads)) + (alpha * batch_size)
        timeToRU = basesRU
        timeToRU += (unmapped_RU * (self.args.mu + self.args.rho))
        timeToRU += (reject_RU * (self.args.mu + self.args.rho))
        timeToRU += (self.args.batch_size * self.args.alpha)
        # BR: all accepted bases and ((mu + rho) * number of unmapped/rejected reads)) + (alpha * batch_size)
        timeToBR = basesBR
        timeToBR += (unmapped_BR * (self.args.mu + self.args.rho))
        timeToBR += (reject_BR * (self.args.mu + self.args.rho))
        timeToBR += (self.args.batch_size * self.args.alpha)
        # update attr
        self.timeNaive += timeToNaive
        self.timeRU += timeToRU
        self.timeBR += timeToBR
        # also keep track of yields
        self.yield_NV += basesNV
        self.yield_RU += (basesRU + ((reject_RU + unmapped_RU) * self.args.mu))
        self.yield_BR += (basesBR + ((reject_BR + unmapped_BR) * self.args.mu))



    def relaunch_checkpoint(self, ckp):
        '''
        SIMULATIONS: launch an experiment with already observed data
        Parameters
        ----------
        ckp: file
            compressed numpy pickle that holds coverage and time passed
        Returns
        -------

        '''
        # replace coverage and times for all conditions
        logging.info(f"reloading checkpoint file: {ckp}")
        if not ckp:
            return
        loaded = np.load(ckp)

        self.timeBR = int(loaded['timeBR'])
        self.timeRU = int(loaded['timeRU'])
        self.timeNaive = int(loaded['timeNV'])
        self.batch = int(loaded['batch'])

        self.update_score(covADD_naive=loaded['covNV'],
                          covADD_RU=loaded['covRU'],
                          covADD_BR=loaded['covBR'])
        print(np.allclose(self.otu.sitewise_eBR, loaded['entropy']))




    def make_decisions_modes(self, paf_full, paf_trunc):
        # NAIVE: make decisions with dummy strategy, i.e. cigar_dict_reject will remain empty
        self.cigars_full, self.cigars_reject, _, _ = self.current_batch.make_decision(
            paf_output=paf_full, otu=self.otu, mode='NV')
        logging.info(f'Naive: accept {len(self.cigars_full)}, reject {len(self.cigars_reject)}')
        # RU: decisions with initial strategy
        self.cigars_accept_RU, self.cigars_reject_RU, self.unmapped_RU, _ = self.current_batch.make_decision(
            paf_output=paf_trunc, otu=self.otu, mode='RU')
        logging.info(f'RU: accept {len(self.cigars_accept_RU)}, reject {len(self.cigars_reject_RU)}')
        # BR: decision with truncated read mappings
        self.cigars_accept_BR, self.cigars_reject_BR, self.unmapped_BR, self.unmapped_ids_BR = self.current_batch.make_decision(
            paf_output=paf_trunc, otu=self.otu, mode='BR')
        logging.info(f'BR: accept {len(self.cigars_accept_BR)}, reject {len(self.cigars_reject_BR)}')



    def write_read_ids(self):
        # PLAYBACK: optionally write the read ids of accepted and rejected reads to a file
        with open(f'{self.args.out_dir}/{self.fq_name}.sim.accepted', 'w') as f:
            accepted_ids = [rid for rid in self.cigars_accept_BR.keys()]
            f.write('\n'.join(accepted_ids))
            f.write('\n')

        with open(f'{self.args.out_dir}/{self.fq_name}.sim.rejected', 'w') as f:
            rejected_ids = [rid for rid in self.cigars_reject_BR.keys()]
            f.write('\n'.join(rejected_ids))
            f.write('\n')
            # also write the unmapped ids to file
            unmapped = list(self.unmapped_ids_BR)
            f.write('\n'.join(unmapped))
            f.write('\n')



    def get_new_coverage(self, paf_raw, workers=12):
        # SIMULATIONS
        # NV: here cigar_dict_full is used for full and accepted (because of accepting everything)
        covADD_naive, _ = self.current_batch.convert_coverage(
            cigars_accept=self.cigars_full,
            cigars_reject=self.cigars_reject,
            cigars_full=self.cigars_full,
            otu=self.otu,
            whole_genome=self.args.wgs,
            workers=workers)

        covADD_RU, self.basesRU = self.current_batch.convert_coverage(
            cigars_accept=self.cigars_accept_RU,
            cigars_reject=self.cigars_reject_RU,
            cigars_full=self.cigars_full,
            otu=self.otu,
            whole_genome=self.args.wgs,
            workers=workers)

        covADD_BR, self.basesBR = self.current_batch.convert_coverage(
            cigars_accept=self.cigars_accept_BR,
            cigars_reject=self.cigars_reject_BR,
            cigars_full=self.cigars_full,
            otu=self.otu,
            whole_genome=self.args.wgs,
            workers=workers)

        return covADD_naive, covADD_RU, covADD_BR



    def update_benefit(self, scores):
        '''
        Wrapper to perform update to the read benefit for each position

        Parameters
        ----------
        scores: np.array
            benefit scores of each position

        Returns
        -------
        benefit: np.array
            read benefit of each position
        s_mu: np.array
            benefit score of the first mu bases

        '''
        # calculate s_mu
        s_mu, scores_placed = self.otu.update_S_mu_thread(
            scores=scores, mu=self.args.mu, window=self.args.window)
        # calculate benefit U
        expected_benefit = self.otu.update_U_thread(
            scores_placed=scores_placed, eta=self.args.eta, window=self.args.window)
        # benefit increase by keeping to sequence instead of rejecting
        benefit = expected_benefit - s_mu
        benefit[benefit < 0] = 0
        return benefit, s_mu




    def find_strat_thread(self, benefit, s_mu, fhat):
        '''
        Finding approximate decision strategy from the current read benefits

        Parameters
        ----------
        benefit: np.array
            read benefit of each position
        s_mu: np.array
            benefit score of the first mu positions
        fhat: np.array
            read starting probabilities

        Returns
        -------
        strat: np.array
            boolean array of decisions for each position
        threshold: float
            benefit value at which reads are accepted from a position

        '''
        # take downsampling into consideration
        window = self.args.window
        alpha = self.args.alpha // window
        rho = self.args.rho // window
        mu = self.args.mu // window
        tc = self.merged_genome.timeCost // window
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
        # target array needs to have largest shape of the thread results
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
        ubar0 = np.sum(fhat * s_mu)
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



    @staticmethod
    def check_bucket_switches(otu, threshold, whole_genome):
        '''
        check at which regions of the genome we apply newly calculated strategies
        regions get new strategies when their mean coverage >= threshold

        Parameters
        ----------
        otu: OTU
            object that holds scores and coverage counts
        threshold: int
            cov depth threshold of when to apply new decision strategies
        whole_genome: bool
            flag for algorithm to use for checking buckets

        Returns
        -------
        switches: np.array
            boolean array to indicate strategy activity for each bucket

        '''
        switches = otu.bucket_switches
        # if not all are on, we check whether we should switch buckets
        if not all(switches):
            csum = np.sum(otu.coverageBR, axis=1)

            if whole_genome:
                # simplified version
                csum_buckets = window_sum(csum, otu.bucket_size)
                cmean_buckets = np.divide(csum_buckets, otu.bucket_size)
                cmean_buckets = adjust_length(switches, cmean_buckets)
            else:
                # put the coverage in buckets
                csum_buckets = np.zeros(shape=switches.shape)
                # sum up the coverage of all sites in a bucket
                np.add.at(csum_buckets, otu.roi_indices // otu.bucket_size, csum)
                # divide by num of sites in buckets
                cmean_buckets = np.divide(csum_buckets, otu.num_roi_in_bucket,
                                          out=np.zeros_like(csum_buckets),
                                          where=otu.num_roi_in_bucket != 0)

            # log mean coverage
            meanBR = np.mean(csum)
            otu.meanBR = meanBR
            # logging.info(f"mean cov: {meanBR}")
            # flip strategy switches
            switches[np.where(cmean_buckets >= threshold)] = 1
            switch_count = np.bincount(switches)
            states = len(switch_count)
            if states == 1:
                logging.info(f"switch count: off {switch_count[0]}; on 0")
            elif states == 2:
                logging.info(f"switch count: off {switch_count[0]}; on {switch_count[1]}")
            else:
                logging.info(f'problem with bucket counts')
        # apply new switch states
        otu.bucket_switches = switches
        return switches



    @staticmethod
    def flush_strat(strat, otu, window=1):
        '''
        NON-BUCKET VERSION
        place new decision strategies into the strategy dict
        to make it available to the otu for effecting decisions

        Parameters
        ----------
        strat: np.array
            boolean array representing the decision strategy
        otu: OTU
            object that holds benefit scores and coverage counts
        window: int
            size of downsampling window

        Returns
        -------

        '''
        # for each chromosome, place the strategy into the correct arrays
        ch_cum_sum = otu.ch_cum_sum // window
        for c, c_ind in otu.chromosome_indices_considered.items():
            chrom_start = ch_cum_sum[c_ind]
            chrom_end = ch_cum_sum[c_ind + 1]
            chrom_strat = strat[chrom_start: chrom_end, :]
            # grab previous strat
            prev_strat = otu.strat_dict_BR[c]
            # adjust size difference from expansion
            lendiff = prev_strat.shape[0] - chrom_strat.shape[0]
            if lendiff > 0:
                # original is longer than replacement
                repl = np.append(chrom_strat, chrom_strat[-lendiff:], axis=0)
            elif lendiff < 0:
                # original is shorter than replacement
                repl = chrom_strat[:-abs(lendiff)]
            else:
                repl = chrom_strat

            assert repl.shape == prev_strat.shape
            # replace strat
            otu.strat_dict_BR[c][:, :] = repl
            # log number of accepted sites
            chrom_len = chrom_strat.shape[0]
            f_accept = np.count_nonzero(chrom_strat[:, 0]) * window
            r_accept = np.count_nonzero(chrom_strat[:, 1]) * window
            f_perc = f_accept / chrom_len
            r_perc = r_accept / chrom_len
            logging.info(f'{c} accept: {f_accept}, {r_accept}; {f_perc}, {r_perc}')



    @staticmethod
    def flush_strat_buckets(strat, otu, switches, window=1):
        '''
        BUCKET VERSION
        place new decision strategies into the strategy dict
        to make it available to the otu for effecting decisions

        Parameters
        ----------
        strat: np.array
            boolean array representing the decision strategy
        otu: OTU
            object that holds benefit scores and coverage counts
        switches: np.array
            boolean array where the strategy is active
        window: int
            size of downsampling window

        Returns
        -------

        '''
        # for each chromosome place the strategy into the correct arrays
        ch_cum_sum = otu.ch_cum_sum // window
        ch_cum_sum_buckets = otu.ch_cum_sum // otu.bucket_size
        for c, c_ind in otu.chromosome_indices_considered.items():
            c_start = ch_cum_sum[c_ind]
            c_end = ch_cum_sum[c_ind + 1]
            c_strat = strat[c_start: c_end, :]
            # grab previous strat
            prev_strat = otu.strat_dict_BR[c]
            # adjust size difference from expansion
            c_strat_adj = adjust_length(original=prev_strat, expanded=c_strat)
            # grab the windowing buckets of this chrom
            bucket_start = ch_cum_sum_buckets[c_ind]
            bucket_end = ch_cum_sum_buckets[c_ind + 1]
            c_switches = switches[bucket_start: bucket_end]
            # expand chromosome switches to windowed genome length
            expand_factor = otu.bucket_size // window
            c_switches_rep = np.repeat(c_switches, expand_factor, axis=0)
            c_switches_adj = adjust_length(original=prev_strat, expanded=c_switches_rep)
            # replace strategy where the switch is on
            otu.strat_dict_BR[c][c_switches_adj, :] = c_strat_adj[c_switches_adj, :]
            # log number of accepted sites
            chrom_len = otu.strat_dict_BR[c].shape[0]
            f_accept = np.count_nonzero(otu.strat_dict_BR[c][:, 0])
            r_accept = np.count_nonzero(otu.strat_dict_BR[c][:, 1])
            f_perc = f_accept / chrom_len
            r_perc = r_accept / chrom_len
            logging.info(f'{c}: {f_accept}, {r_accept}; {f_perc}, {r_perc}')



    def update_strategy(self):
        '''
        wrapper for updating read benefits and the decision strategy

        Returns
        -------

        '''
        otu = self.otu
        # flip strategy switches if threshold is reached
        switches = BossRun.check_bucket_switches(
            otu=otu, threshold=self.args.cov_until, whole_genome=self.args.wgs)
        # UPDATE STRATEGY
        if any(switches):
            # update Fhat: unpack and normalise
            fhat = self.merged_genome.updateFpointmass()
            fhat_exp = self.merged_genome.expandF(fhat=fhat, downsample_window=self.args.window)
            # update S_mu and U
            benefit, s_mu = self.update_benefit(scores=otu.S_BR)
            self.benefit = benefit  # tmp for tracking
            # find new strategy
            strat, threshold = self.find_strat_thread(benefit=benefit, s_mu=s_mu, fhat=fhat_exp)
            self.threshold = threshold
            otu.threshold = threshold

            # replace strategies in dictionaries
            if len(switches) > 1:
                BossRun.flush_strat_buckets(strat=strat, otu=otu, switches=switches, window=self.args.window)
            else:
                BossRun.flush_strat(strat=strat, otu=otu, window=self.args.window)
            # save masks to disk
            otu.save_strat()



    def pickle_ckp(self):
        '''
        SIMULATIONS: Save coverage counts and sequencing time spent per condition
        Used to relaunch simulations with some initial data

        Returns
        -------

        '''
        np.savez_compressed(f"{self.args.out_dir}/ckp/ckp_{self.current_batch.num}",
                            covBR=self.otu.coverageBR,
                            covRU=self.otu.coverageRU,
                            covNV=self.otu.coverageNaive,
                            timeBR=self.timeBR,
                            timeRU=self.timeRU,
                            timeNV=self.timeNaive,
                            batch=self.batch,
                            entropy=self.otu.sitewise_eBR,
                            entropyRU=self.otu.sitewise_eRU,
                            entropyNV=self.otu.sitewise_eNV,
                            scoreBR=self.otu.S_BR,
                            threshold=self.threshold)




    def process_batch(self, fastq_file=None, sampler=None):
        # SIMULATIONS
        logging.info(f"\n\n\n Next batch ---------------------------- # {self.batch}")
        tic = time.time()

        if sampler:
            # sample fastq and grab their mappings
            paf_full, paf_trunc = self.read_and_inject(sampler=sampler)
        else:
            # sample from fastq and map reads
            paf_full, paf_trunc = self.read_and_map(fastq_file=fastq_file, workers=16)

        # making decisions
        self.make_decisions_modes(paf_full=paf_full, paf_trunc=paf_trunc)
        # self.write_read_ids()  # PLAYBACK

        # convert coverage
        covADD_naive, covADD_RU, covADD_BR = self.get_new_coverage(paf_raw=paf_full, workers=8)

        # update times
        self.update_times(
            basesNV=self.current_batch.basesTOTAL,
            basesRU=self.basesRU,
            basesBR=self.basesBR,
            unmapped_RU=self.unmapped_RU,
            unmapped_BR=self.unmapped_BR,
            reject_RU=len(self.cigars_reject_RU),
            reject_BR=len(self.cigars_reject_BR))

        # update the observed abundances
        self.abundance_tracker.update(n=len(self.current_batch.read_ids), cigar_dict=self.cigars_full)

        # periodically update read length dist
        self.otu.update_length_dist(batch=self.batch, cigars_full=self.cigars_full)
        # update scores
        self.update_score(covADD_naive=covADD_naive, covADD_RU=covADD_RU, covADD_BR=covADD_BR)

        # count read starting positions
        self.merged_genome.countReadStarts(cigars_full=self.cigars_full)

        logging.info(f'time Naive: {self.timeNaive}')
        logging.info(f'time RU: {self.timeRU}')
        logging.info(f'time BR: {self.timeBR}')
        logging.info(f'yield Naive: {self.yield_NV / 1e9}')
        logging.info(f'yield RU: {self.yield_RU / 1e9}')
        logging.info(f'yield BR: {self.yield_BR / 1e9}')

        # update strategy
        self.update_strategy()

        # save metrics at end
        if (self.batch) % 10 == 0:
            self.pickle_ckp()
            self.otu.save_metrics_chroms(
                batch=self.batch, timeNV=self.timeNaive, timeRU=self.timeRU, timeBR=self.timeBR,
                yield_NV=self.yield_NV, yield_RU=self.yield_RU, yield_BR=self.yield_BR)

        self.batch += 1
        toc = time.time()
        logging.info(f"batch took: {toc - tic}")





class BossRun_live(BossRun):


    def __init__(self, args):
        const = Constants()
        for c, cval in const.__dict__.items():
            setattr(args, c, cval)

        # development option for testing
        if args.testing:
            args.cov_until = 0
            # out_path = "."
            # args.fastq_dir = "."

        self.args = args
        self.processed_files = set()
        self.batch = 0
        self.threshold = 0

        # set output directory for bossruns in cwd
        # make sure the run name does not have any spaces
        assert ' ' not in args.run_name

        args.out_dir = f'./bossruns_{args.run_name}'
        if not os.path.exists(args.out_dir):
            os.mkdir(args.out_dir)
            os.mkdir(f'{args.out_dir}/masks')
            os.mkdir(f'{args.out_dir}/ckp')
            os.mkdir(f'{args.out_dir}/fq')
            os.mkdir(f'{args.out_dir}/fq/RU')
            os.mkdir(f'{args.out_dir}/fq/BR')
            os.mkdir(f'{args.out_dir}/fq/NV')
            os.mkdir(f'{args.out_dir}/metrics')

        # initialise logging file
        init_logger(logfile=f'{args.out_dir}/{args.run_name}.bossruns.log')

        # print args for debugging
        logging.info("BOSS-RUNS")
        for a, aval in self.args.__dict__.items():
            logging.info(f'{a} {aval}')
        logging.info('\n')

        # initialise sequencing device dependent paths
        logging.info("looking for MinKNOW's output path..")
        if not args.fastq_dir:
            try:
                out_path = grab_output_dir(device=args.device, host=args.host, port=args.port)
                logging.info(f"grabbing MinKNOW's output path: \n{out_path}\n")
                args.fastq_dir = f'{out_path}/fastq_pass'
            except:
                logging.info("MinKNOW's output dir could not be inferred from device name. Exiting..")
                logging.info(f'device: {args.device}, host: {args.host}, port: {args.port}')
                # out_path = "/home/lukas/Desktop/BossRuns/playback_target/data/pb01/no_sample/20211021_2209_MS00000_f1_f320fce2"
                # args.fastq_dir = out_path
                sys.exit()

        # grab channels of the condition - only relevant if conditions
        if args.conditions:
            channel_path = f'{out_path}/channels.toml'
            logging.info(f'looking for channels specification at : {channel_path}')
            channels_found = False
            while not channels_found:
                if not os.path.isfile(channel_path):
                    logging.info("channels file does not exist (yet), waiting for 30s")
                    time.sleep(30)
                else:
                    self.channels = grab_br_channels(channels_toml=channel_path, run_name=args.run_name)
                    channels_found = True
            # channels successfully found
            logging.info(f"found channels specification: BR uses {len(self.channels)} channels.")
        else:
            # if we use a whole flowcell, use all channels
            self.channels = set(np.arange(1, 512 + 1))
        # initialise the abundance tracker
        self.abundance_tracker = AbundanceTracker()
        # check which scenario
        args.wgs = 1  # default
        if args.vcf:
            args.wgs = 0




    def scan_dir(self):
        '''
        Periodically scan the sequencers output directory for new files
        gets triggered after some defined time has passed (default 90s)
        creates new batch from all NEW files

        Returns
        -------

        '''
        patterns = ["*.fq.gz", "*.fastq.gz", "*.fastq.gzip", "*.fq.gzip", "*.fastq", "*.fq"]
        all_fq = set()
        for p in patterns:
            all_fq.update(glob.glob(f'{self.args.fastq_dir}/{p}'))

        # use only new files
        new_fq = all_fq.difference(self.processed_files)
        logging.info(f"found {len(new_fq)} new fq files: \n {new_fq}")
        # add new files to processed files
        self.processed_files.update(new_fq)
        return list(new_fq)


    def read_and_map(self, new_fastq, workers=12):
        '''
        Ingest the next batch of sequencing data
        This can be multiple new fastq files that were deposited
        by the sequencer in the last sequencing period

        Parameters
        ----------
        new_fastq: list
            list of paths to new fastq files
        workers: int
            threads to use for mapping reads

        Returns
        -------
        paf_raw: str
            Alignments of reads to the reference in PAF format
        '''
        self.current_batch = CurrentBatch(
            num=self.batch, qt=self.args.base_quality_threshold, window=self.args.window)

        # read all new files
        read_ids = set()
        read_lengths = dict()
        read_sequences = dict()
        read_qualities = dict()
        basesTOTAL = 0

        for fastq_file in new_fastq:
            rids, rlen, rseq, rqual, bTOTAL = self.current_batch.read_batch(fastq_file=fastq_file,
                                                                            channels=self.channels)
            read_ids.update(rids)
            read_lengths.update(rlen)
            read_sequences.update(rseq)
            read_qualities.update(rqual)
            basesTOTAL += bTOTAL

        self.current_batch.read_ids = read_ids
        self.current_batch.read_lengths = read_lengths
        self.current_batch.read_sequences = read_sequences
        self.current_batch.read_qualities = read_qualities
        self.current_batch.basesTOTAL = basesTOTAL
        # map the full reads
        paf_raw = self.mapper.mappy_batch(
            sequences=self.current_batch.read_sequences, truncate=False, workers=workers)
        return paf_raw


    def get_new_coverage(self, paf_raw, workers):
        # convert paf output
        self.cigar_dict = self.current_batch.convert_paf(paf_raw=paf_raw)
        # convert mappings to coverage counts
        # uses same function as when making decisions ourselves
        # pass same dict as full and accept
        coverage_addition, _ = self.current_batch.convert_coverage(
            cigars_accept=self.cigar_dict,
            cigars_reject=dict(),
            cigars_full=self.cigar_dict,
            otu=self.otu,
            whole_genome=self.args.wgs,
            workers=workers)
        return coverage_addition


    def update_score(self, coverage_addition):
        logging.info("updating scores")
        otu = self.otu
        # BR
        otu.S_BR, otu.sitewise_eBR = update_S_scoreArray(
            otu=otu,
            scores=otu.S_BR,
            entropy=otu.sitewise_eBR,
            coverage=otu.coverageBR,
            covADD=coverage_addition)

        otu.S_BR, otu.sitewise_eBR, otu.eBR = modify_scores(
            cov=otu.coverageBR,
            scores=otu.S_BR,
            ent=otu.sitewise_eBR,
            chrom_rpos=otu.chrom_rpos_ranges)

        # add new coverage counts to previous coverage
        otu.add_coverage(covADD_BR=coverage_addition)


    def relaunch_checkpoint(self, ckp):
        '''
        Used to relaunch experiments with some initial data
        E.g. when using multiple flowcells sequenctially

        Parameters
        ----------
        ckp: file
            checkpoint file created by pickle_ckp()

        Returns
        -------

        '''
        # replace empty coverage with saved counts
        logging.info(f"reloading checkpoint file: {ckp}")
        if not ckp:
            return
        loaded = np.load(ckp)
        # load batch number
        self.batch = int(loaded['batch'])
        # use loaded coverage counts to calculate scores
        self.update_score(coverage_addition=loaded['covBR'])
        # check if calculated scores result in same entropy as in the ckp
        print(np.allclose(self.otu.sitewise_eBR, loaded['entropy']))


    def pickle_ckp(self):
        '''
        Save coverage counts so far observed by BOSS-RUNS
        Used to relaunch experiments with some initial data
        E.g. when using multiple flowcells sequenctially

        Returns
        -------

        '''
        np.savez_compressed(f"{self.args.out_dir}/ckp/ckp_{self.current_batch.num}",
                            covBR=self.otu.coverageBR,
                            batch=self.batch,
                            entropy=self.otu.sitewise_eBR,
                            threshold=self.threshold)



    def process_batch(self):
        # LIVE VERSION
        logging.info(f"\n\n\n Next batch ---------------------------- # {self.batch}")
        tic = time.time()

        # find new fastq files
        new_fastq = self.scan_dir()
        if not new_fastq:
            logging.info("no new files, deferring update.. ")
            return self.args.wait

        # read fq of the batch and map
        paf_raw = self.read_and_map(new_fastq=new_fastq)
        coverage_addition = self.get_new_coverage(paf_raw=paf_raw, workers=12)
        self.abundance_tracker.update(n=len(self.current_batch.read_ids), cigar_dict=self.cigar_dict)
        self.otu.update_length_dist(batch=self.batch, cigars_full=self.cigar_dict)
        self.update_score(coverage_addition=coverage_addition)
        self.merged_genome.countReadStarts(cigars_full=self.cigar_dict)
        self.update_strategy()

        if (self.batch) % 10 == 0:
            self.pickle_ckp()
            self.otu.save_metrics_acc(batch=self.batch)

        toc = time.time()
        passed = toc - tic
        next_update = int(self.args.wait - passed)
        logging.info(f"batch took: {passed}")
        logging.info(f"finished updating masks, waiting for {next_update} ... \n")
        self.batch += 1
        return next_update


    def process_batch_analysis(self, fastq):
        # VERSION FOR ANALYSIS
        # RUNS THROUGH FASTQ BATCHES
        logging.info(f"\n\n\n Next batch ---------------------------- # {self.batch}")
        # read fq of the batch and map
        paf_raw = self.read_and_map(new_fastq=fastq)
        coverage_addition = self.get_new_coverage(paf_raw=paf_raw, workers=12)
        self.update_score(coverage_addition=coverage_addition)
        self.batch += 1
        return self.otu.coverageBR, self.otu.S_BR, self.otu.sitewise_eBR






def get_reference(reference_file_path, soft_masking=False):
    '''
    Load the reference sequence and prepare a dictionary of all chromosomes

    Parameters
    ----------
    reference_file_path: str
        Path to the reference file
    soft_masking: bool
        whether to ignore soft masked sites when counting read start positions
        can avoid problems if many TEs map to the same place

    Returns
    -------
    genome: np.array
        integer array of the genome (concatenated chromosomes)
    soft_mask: set
        indices of where to ignore read starts
    chNum: int
        number of chromosomes in the reference
    chArray: np.array
        array of chromosome lengths
    '''
    # Read in chromosome file to get lengths and stats
    chr_seqs = dict()
    chr_seqs_raw = dict()

    with open(reference_file_path, "r") as fh:
        for header, seq in read_fa(fh):
            chr_seqs[header] = seq.upper()
            chr_seqs_raw[header] = seq

    # transform the genome into integers
    genome = integer_genome(chromosome_sequences=chr_seqs)
    unknown_bases = np.where(genome == 99)[0].shape
    if unknown_bases[0] > 0:
        logging.info("unknown characters found in genome")

    if soft_masking:
        # find where genome was soft masked
        # only used when counting read start sites
        genome_sm = integer_genome(chromosome_sequences=chr_seqs_raw)
        soft_mask_ind = np.where(genome_sm == 99)[0]
        # add flanking sites to account for mapping uncertainty
        soft_mask = set()
        if soft_mask_ind.shape != (0,):
            soft_mask.update(soft_mask_ind)
            to_be_added = set()
            for s in soft_mask:
                to_be_added.update(np.arange(s - 2, s + 2))
            soft_mask.update(to_be_added)
    else:
        soft_mask = set()

    chNum = len(chr_seqs)
    chArray = np.asarray([len(seq) for seq in chr_seqs.values()], dtype=np.dtype("uint32"))
    return genome, soft_mask, chNum, chArray



def sitewise_score(scores, pos_posterior, _len_g, _len_b, phi):
    """
    Uses the current posterior probabilities to find the benefit score of a given position.

    Parameters
    ----------
    scores: np.ndarray
        Array of scores to update
    pos_posterior: np.ndarray
        posterior probabilites to calculate the scores
    _len_g: int
        Number of possible genotypes
    _len_b: int
        Number of possible bases
    phi: np.ndarray
        observation probabilities phi

    Returns
    -------
    scores: np.ndarray
        Updated array of scores
    entropy: np.ndarray
        Updated array of entropies

    """
    len_s = len(scores)
    # Shannon entropy of the current posterior distribution.
    # "entropy" is part of the calculation of S(i) in section 0.1.1
    logs = np.log(pos_posterior, where=pos_posterior > 0.0)
    logs[logs == np.nan] = 0

    entropy = np.sum(-pos_posterior * logs, axis=1)
    # Expected entropy after a new batch of reads ("newEntropy").
    # this is also part of the calculation of S(i) in section 0.1.1
    # (part of the final line of the big series of equations)
    new_entropy = np.zeros(len_s)
    # probability of observing a new base of each type at the considered position
    # these observation_probabilities are the P(d_{n+1,i}) currently defined in equation (6) in section0.1.1
    observation_probabilities = np.zeros(len_s)
    # posterior probabilities after reading a certain base i
    # these are the f_i(b|D') calculated as in the last line in equation (4) of section 0.1.1
    new_post = np.zeros((len_s, _len_g))
    for i in range(_len_b):  # new read base
        np.multiply(pos_posterior, phi[i], out=new_post)
        np.sum(new_post, axis=1, out=observation_probabilities)
        observation_probabilities[observation_probabilities == 0] = 1e-300
        new_post /= observation_probabilities[:, np.newaxis]
        np.log(new_post, where=new_post > 0.0, out=logs)
        logs[logs == np.nan] = 0
        for j in range(_len_g):  # genotype
            new_entropy -= observation_probabilities * new_post[:, j] * logs[:, j]
    # update scores
    scores[:] = entropy - new_entropy
    return scores, entropy



def calcCoveragePost(prior, phi_stored, targetCoverage):
    """
    Calculate posterior for site patterns, i.e. not updating posterior
    probabilities with counts but rather from prior straight to a
    certain coverage. Used as private function to create multidimensional
    array for performing lookups.

    Parameters
    ----------
    prior: np.ndarray
        Array that contains the priors for each possible base
    phi_stored: np.ndarray
        Array that stores phi for different counts of each base.
    targetCoverage: np.ndarray
        Array with different patterns of coverage for which to calc the posterior

    Returns
    -------
    coveragePost: np.ndarray
        3D array with all posteriors for the coverage patterns that were input
        shape is 4 x 5 x n

    """
    # prevent overflow when indexing into phi_stored - hacky but needed
    targetCoverage[targetCoverage > 990] = 990

    # len_b = prior.shape[0]  # this is wrong if gap is considered as base
    # len_g = prior.shape[1]

    len_b = phi_stored.shape[0]
    len_g = phi_stored.shape[1]

    # how many patterns of coverage are input
    n = targetCoverage.shape[0]
    # repeat the priors for all input patterns
    coveragePost = np.repeat(prior[:, np.newaxis], repeats=n, axis=1)

    # init empty phi, since we will use the stored version
    phis = np.full(n, 1.0)

    # calculate the posteriors for all coverage patterns
    for j in range(len_g):
        if j > 0:
            phis.fill(1.0)
        for i in range(len_b):
            phis *= phi_stored[i, j, targetCoverage[:, i]]
        for h in range(4):
            coveragePost[h, :, j] *= phis

    # sumB is Z_i(D) in manuscript section 0.1.1. Normalizing factor
    # for the posterior (the likelihood of the new data at the position).
    for _h in range(4):
        sum_b = np.sum(coveragePost[_h, :, :], axis=1)
        sum_b[sum_b < 1e-300] = 1e-300  # prevent the next line to divide by zero
        coveragePost[_h, :, :] /= sum_b[:, np.newaxis]

    return coveragePost



def update_S_scoreArray(otu, scores, entropy, coverage, covADD, array_limit=30, model=0):
    '''
    Main function to update benefit scores using a pre-computed array

    Parameters
    ----------
    otu: OTU
        OTU that holds coverage and score info about the genomes of interest
    scores: np.array
        benefit scores of each site
    entropy: np.array
        entropy values of each site
    coverage: np.array
        5xN array that holds the observed coverage values used to calc benefit scores
    covADD:
        newly observed coverage used for updates
    array_limit: int
        threshold at which benefit scores are not updated anymore
    model: int
        experimental flag to use alternative model

    Returns
    -------
    scores: np.array
        updated benefit scores of each site
    entropy: np.array
        updated entropy values of each site

    '''
    # add newly observed coverage to previous
    targetCoverage = np.add(coverage, covADD)
    # find updated indices
    change_mask = np.where(np.sum(covADD, axis=1) >= 1, True, False)
    # if deletions as genotype are ignored, set to 0
    if otu.len_b == 4:
        targetCoverage[:, 4] = 0

    # indices where array_limit is reached
    maxed_indices = np.where(targetCoverage.sum(axis=1) >= array_limit)[0]
    change_mask[maxed_indices] = False
    # indices of change
    cc_pos = np.nonzero(change_mask)[0]
    # target coverage patterns
    cC = targetCoverage[cc_pos]
    # grab base of changed positions
    refBases = otu.genome[cc_pos]
    # main operation: indexing in pre-computed values using coverage patterns
    scores[cc_pos] = otu.scoreArray[cC[:, 0], cC[:, 1], cC[:, 2], cC[:, 3], cC[:, 4], refBases]
    # prevent maxed sites from being recalculated
    scores[maxed_indices] = np.finfo(float).tiny
    # check positions where scores were not found in the scoreArray
    # i.e. if they remain 0.0
    missing = np.argwhere(scores == 0.0).flatten()
    nmiss = missing.shape[0]

    # if any positions without scores, calc and add to array
    if nmiss != 0:
        missingPatterns = targetCoverage[missing]
        # calculate new scores and entropies
        missingEntropies, missingScores = otu.calcPostAndScore(coveragePatterns=missingPatterns, model=model)
        # place into pre-computed array
        bases = otu.genome[missing]
        scores[missing] = missingScores[bases, np.arange(nmiss)]
        entropy[missing] = missingEntropies[bases, np.arange(nmiss)]
        ind = np.swapaxes(missingPatterns, 0, 1)
        # assign the calculated scores and entropies to the large array
        for i in range(ind.shape[1]):
            otu.scoreArray[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = missingScores[:, i]
            otu.entropyArray[ind[0, i], ind[1, i], ind[2, i], ind[3, i], ind[4, i]] = missingEntropies[:, i]

    # grab the entropy values of the changed sites
    entropy[cc_pos] = otu.entropyArray[cC[:, 0], cC[:, 1], cC[:, 2], cC[:, 3], cC[:, 4], refBases]
    return scores, entropy



def modify_scores(cov, scores, ent, chrom_rpos):
    '''
    Modify benefit scores after updating, specifically to ignore large deletions
    or missing regions in the sampled genome, e.g. accessory genes
    For each chromosome separately, in case there is abundance bias

    Parameters
    ----------
    cov: np.array
        coverage counts of all positions
    scores: np.array
        benefit scores of each position
    ent: np.array
        entropy values of each position
    chrom_rpos: dict
        positions of ROIs within genomes

    Returns
    -------
    scores: np.array
        updated benefit scores of each position
    ent: np.array
        updated entropy values of each position
    ent_sum: float
        total remaining entropy

    '''
    for chrom, rpos in chrom_rpos.items():
        croi_start, croi_end = rpos

        chrom_cov = cov[croi_start: croi_end]
        chrom_covsum = np.sum(chrom_cov, axis=1)

        if np.mean(chrom_covsum) > 5:
            dropout_idx = _find_dropout(chrom_covsum)
            logging.info(f'detected {dropout_idx.shape[0]} dropouts')
            # copies to avoid double fancy-indexing
            chrom_scores = np.copy(scores[croi_start: croi_end])
            chrom_scores[dropout_idx] = 0  # np.finfo(float).tiny
            # optionally also set entropy to 0
            # chrom_ent = np.copy(ent[croi_start: croi_end])
            # chrom_ent[dropout_idx] = 0  # np.finfo(float).tiny
            # assign back to original
            scores[croi_start: croi_end] = chrom_scores
            # ent[croi_start: croi_end] = chrom_ent

    # calculate sum of remaining entropy
    ent_sum = np.nansum(ent)
    return scores, ent, ent_sum



def _find_dropout(covsum, mod=8):
    '''
    If there are sites that have not had any coverage after some time
    we don't expect them to gain any at all and ignore their benefit

    Parameters
    ----------
    covsum: np.array
        coverage depth at each position
    mod: int
        threshold modifier

    Returns
    -------
    dropout: np.array
        array of positions that will most likely not get any coverage

    '''
    cov_mean = np.mean(covsum)
    # ignore threshold is dependent on mean coverage
    threshold = int(cov_mean / mod)
    dropout = np.where(covsum <= threshold)[0]
    return dropout



def update_U_forward(scores, scores_div, approx_CCL, eta, linear):
    ch_len = scores.shape[0]
    # FORWARD STARTING VALUE: sum of scores from first partition
    start_val = np.sum(scores[: approx_CCL[0]])
    # rest of the segments of piece-wise constant f
    for bd in range(len(approx_CCL) - 1):
        start_val += np.sum(scores_div[approx_CCL[bd]: approx_CCL[bd + 1]]) * (eta - 2 - bd)
    # forward benefits
    to_be_added = np.zeros(shape=(eta + 1, ch_len))
    to_be_added[0, :] = -scores

    for bd_index in range(eta - 1):
        bd = approx_CCL[bd_index]
        to_be_added[1 + bd_index, : ch_len - bd] = scores_div[bd:]
        if not linear:
            to_be_added[1 + bd_index, ch_len - bd :] = scores_div[: bd]

    site_sums = np.sum(to_be_added, axis=0)
    to_sum = np.concatenate((np.array([start_val]), site_sums[:-1]))
    forward_u = np.cumsum(to_sum)
    return forward_u



def update_U_reverse(scores, scores_div, approx_CCL, eta, linear):
    ch_len = scores.shape[0]
    # expected benefit for reverse reads
    start_val = scores[0]
    if not linear:
        start_val += np.sum(scores[- approx_CCL[0] + 1:])
        for bd in range(len(approx_CCL) - 1):
            start_val += np.sum(scores_div[- approx_CCL[bd + 1] + 1: - approx_CCL[bd] + 1]) * (eta - 2 - bd)

    # reverse benefits
    to_be_added = np.zeros(shape=(eta + 1, ch_len))
    # to_be_added[0] = start_val
    to_be_added[0, : -1] = scores[1:]

    for bd_index in range(eta - 1):
        bd = approx_CCL[bd_index]
        to_be_added[1 + bd_index, bd - 1:] = -scores_div[: - bd + 1]
        if not linear:
            to_be_added[1 + bd_index, : bd - 1] = -scores_div[- bd + 1:]

    site_sums = np.sum(to_be_added, axis=0)
    to_sum = np.concatenate((np.array([start_val]), site_sums[:-1]))
    reverse_u = np.cumsum(to_sum)
    return reverse_u



def update_Smu_chr(scores, mu, linear):
    # create a new array of 2 * chromosome length
    # instead of removing and adding one value at the beginning and at the
    # end of each position, we have an array of twice the length with
    # negative and positive values, which can be used with cumsum
    chr_len = scores.shape[0]
    to_be_summed = np.zeros(2 * chr_len)
    # score at position 0 (0:mu)
    tbs0 = np.nansum(scores[0: mu])
    to_be_summed[0] = tbs0
    # from pos 2 to 2N - 2mu + 2 with step 2
    # transfer the score of bases mu until the end (only to every other element)
    to_be_summed[2: 2 * chr_len - 2 * mu + 2: 2] = scores[mu:]
    # from pos 3 until the end, with step 2
    # transfer the negative score of pos 0 until end-1
    to_be_summed[3::2] = -scores[0 : -1]  # or end -1 or -2?
    # now to_be_summed has a high value at 0 and then fluctuates around 0
    # with the pos and neg of scores
    if not linear:
        # from 2N - 2mu+2 until the end, with step 2
        # transfer scores from pos 0 until pos mu-1
        # smooth out edge effects at the end of chromosome
        to_be_summed[2 * chr_len - 2 * mu + 2:: 2] = scores[0: mu - 1]

    # array that fluctuates around the previous value at pos 0
    # with +- score value
    partial_sums = np.cumsum(to_be_summed)
    return partial_sums, tbs0



def calc_u(arguments):
    # helper function for multithreading
    func, scores_placed, ccl, eta = arguments
    sdiv = scores_placed / (eta - 1)
    u = func(scores=scores_placed,
             scores_div=sdiv,
             approx_CCL=ccl, eta=eta, linear=True)
    return u



def calc_smu(arguments):
    # helper function for multithreading
    s, mu = arguments
    smu, ch_start = update_Smu_chr(scores=s, mu=mu, linear=True)
    return smu, ch_start



def grab_output_dir(device, host='localhost', port=None):
    '''
    Capture the output directory of MinKNOW,
    i.e. where fastq files are deposited during sequencing
    host and port should be optional if BOSS-RUNS is run on the sequencing machine

    Parameters
    ----------
    device: str
        device name of the 'position' in the sequencing machine
    host: str
        hostname to connect to for MinKNOW
    port: int
        override default port to connect to

    Returns
    -------

    '''
    logging.info(f"minknow API Version {minknow_api_version}")
    # minknow_api.manager supplies Manager (wrapper around MinKNOW's Manager gRPC)
    if minknow_api_version.startswith("5"):
        if not port:
            port = 9502
        manager = Manager(host=host, port=port)
    elif minknow_api_version.startswith("4"):
        if not port:
            port = 9501
        manager = Manager(host=host, port=port, use_tls=False)
    else:
        logging.info("unsupported version of minknow_api")
        sys.exit()
    # Find a list of currently available sequencing positions.
    positions = list(manager.flow_cell_positions())
    pos_dict = {pos.name: pos for pos in positions}
    # index into the dict of available devices
    try:
        target_device = pos_dict[device]
    except KeyError:
        logging.info(f"Error: target device {device} not available")
        logging.info("Error: Please make sure to supply correct name of sequencing position in MinKNOW.")
        sys.exit()
    # connect to the device and navigate api to get output path
    device_connection = target_device.connect()
    current_run = device_connection.protocol.get_current_protocol_run()
    run_id = current_run.run_id
    logging.info(f"connected to run_id: {run_id}")
    out_path = current_run.output_path
    return out_path







