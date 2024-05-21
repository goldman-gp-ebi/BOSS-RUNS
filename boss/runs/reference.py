import logging
from collections import defaultdict
from pathlib import Path
from string import ascii_letters

import numpy as np
from numpy.typing import NDArray
import bottleneck as bn

from boss.utils import read_fa, window_sum, adjust_length
from boss.mapper import Indexer
from boss.runs.sequences import Scoring




class Contig:

    def __init__(self, name: str, seq: str, ploidy: int = 1, rej: bool = False):
        """
        Initialise a contig object

        :param name: Name as in fasta header
        :param seq: DNA string of the contig.
        :param ploidy: Contig from haploid or diploid organism
        :param rej: Always reject from this contig?
        """
        self.name = name.strip().split(" ")[0]
        self.seq = seq.upper()
        self.length = len(self.seq)
        self.rej = rej   # flag whether to reject all reads from this contig
        self.seq_int = self._seq2int()
        self._init_coverage()
        self._init_buckets()
        self.scoring = Scoring(ploidy=ploidy)
        self.len_b = self.scoring.priors.len_b
        self.score0, self.ent0 = self.scoring.score0, self.scoring.ent0
        self._init_scores()
        self._init_strat()


    def _seq2int(self) -> NDArray:
        """
        Transform nucleotide strings to integer arrays.
        Unknown bases are set to integer 0.

        :return: Array of integers representing nucleotides
        """
        # define the dict for translating the bases
        transDict = defaultdict(str)
        transDict.update({'A': '0', 'C': '1', 'G': '2', 'T': '3'})
        for letter in ascii_letters:
            if letter not in transDict.keys():
                transDict[letter] = '0'
        # create translation table from dictionary
        base2int = str.maketrans(dict(transDict))
        # translate and convert from byte string
        read_integer = self.seq.translate(base2int)
        seq_int = np.frombuffer(read_integer.encode(), 'u1') - ord('0')
        assert set(seq_int) == {0, 1, 2, 3}
        return seq_int


    def _init_coverage(self) -> None:
        """
        Initialise containers to hold the observed coverage values

        :return:
        """
        self.coverage = np.zeros(shape=(self.length, 5), dtype="uint16")
        # to indicate where changes happened for score calculations
        self.change_mask = np.zeros(shape=self.length, dtype="bool")



    def _init_buckets(self, bucket_size: int = 20_000) -> None:
        """
        Each contig has multiple 'buckets' in which the strategy activates independently

        :param bucket_size: Size of the bucket windows
        :return:
        """
        self.bucket_size = bucket_size
        self.bucket_switches = np.ones(shape=int(self.length // bucket_size) + 1, dtype="bool")
        self.switched_on = False



    def _init_scores(self) -> None:
        """
        Initialise the benefit scores at each site using the initial priors

        :return:
        """
        # initialise entropy
        initial_ent = np.zeros(shape=self.length)
        initial_ent.fill(self.ent0[0])
        # entropy_sum = np.sum(initial_ent)
        # score array
        self.initial_scores = np.zeros(shape=self.length)
        self.initial_scores.fill(self.score0[0])
        # create copies of the scores, entropy & sum of entropy
        self.scores = np.copy(self.initial_scores)
        self.entropy = np.copy(initial_ent)
        # self.eBR = entropy_sum



    def _init_strat(self, window: int = 100) -> None:
        """
        Initialise dummy strategies for the contig
        :param window: Downsampling window size
        :return:
        """
        if self.rej:
            self.strat = np.zeros(dtype="bool", shape=(self.length // window, 2))
        else:
            self.strat = np.ones(dtype="bool", shape=(self.length // window, 2))



    def increment_coverage(self, increment_list: list[int, int, NDArray, NDArray]) -> None:
        """
        Increment the coverage counts after parsing mappings

        :param increment_list: List of increment instructions for this contig
        :return:
        """
        # reset the change mask
        self.change_mask.fill(0)
        # temporary container for coverage
        tmp_cov = np.zeros(shape=self.coverage.shape, dtype="uint16")
        for (start, end, query_arr, addition) in increment_list:
            # range to index into the bit we want to increment
            indices = np.arange(query_arr.shape[0])
            # add the "addition" to the corresponding bases
            np.add.at(tmp_cov[start: end], (indices, query_arr), addition)
        # set the change mask
        self.change_mask[np.where(tmp_cov)[0]] = 1
        # add the coverage to the previous one
        self.coverage += tmp_cov



    def modify_scores(self) -> None:
        """
        Modify benefit scores after updating, specifically to ignore large deletions
        or missing regions in the sampled genome, e.g. accessory genes
        For each chromosome separately, in case there is abundance bias

        :return:
        """
        # coverage depth
        covsum = np.sum(self.coverage, axis=1)
        if np.mean(covsum) > 5:
            dropout_idx = self._find_dropout(covsum)
            logging.info(f'detected {dropout_idx.shape[0]} dropouts')
            self.scores[dropout_idx] = 0



    @staticmethod
    def _find_dropout(covsum: NDArray, mod: int = 8) -> NDArray:
        """
        If there are sites that have not had any coverage after some time
        we don't expect them to gain any at all and ignore their benefit

        :param covsum: coverage depth at each position
        :param mod: Threshold modifier
        :return: Array of dropout indices
        """
        cov_mean = np.mean(covsum)
        # ignore threshold is dependent on mean coverage
        threshold = int(cov_mean / mod)
        dropout = np.where(covsum <= threshold)[0]
        return dropout



    def check_buckets(self) -> None:
        """
        Check if strategy needs to be activated in some of the buckets of a contig
        This happens if the mean coverage in a bucket is greater than the threshold
        :return:
        """
        # coverage depth
        csum = np.sum(self.coverage, axis=1)
        # coverage in buckets
        csum_buckets = window_sum(csum, self.bucket_size)
        cmean_buckets = np.divide(csum_buckets, self.bucket_size)
        cmean_buckets = adjust_length(original_size=self.bucket_switches.shape[0], expanded=cmean_buckets)
        # flip strategy switches
        threshold = 5
        self.bucket_switches[np.where(cmean_buckets >= threshold)] = 1
        switch_count = np.bincount(self.bucket_switches)
        states = len(switch_count)
        # log the first time a contig's strategy is switched on
        if states == 2 and not self.switched_on:
            self.switched_on = True
            logging.info(f"Activated strategy for: {self.name}")



    def calc_smu(self, window: int = 100, mu: int = 400) -> None:
        """
        Calculate S_mu, i.e. scores of mu-sized fragment at each position

        :param window: Size of downsampling window
        :param mu: Length of mu in model
        :return:
        """
        # downsample the scores
        self.scores_ds = np.zeros(shape=int(self.length // window) + 1)
        site_indices = np.arange(0, self.length) // window
        # avoid buffering
        np.add.at(self.scores_ds, site_indices, self.scores)
        # calculate smu - fwd needs double reversal due to how bn.move_sum() operates
        smu_fwd = bn.move_sum(self.scores_ds[::-1], window=mu // window, min_count=1)[::-1]
        smu_rev = bn.move_sum(self.scores_ds, window=mu // window, min_count=1)
        # assign smu as attribute
        self.smu = np.zeros(shape=(self.scores_ds.shape[0], 2))
        self.smu[:, 0] = smu_fwd
        self.smu[:, 1] = smu_rev



    def calc_u(self, approx_ccl: NDArray, window: int = 100) -> None:
        """
        Calculate the expected benefit of new fragments at each position
        New implementation using bn.move_sum

        :param approx_ccl: Approximation of read length distribution
        :param window: Downsampling window size
        :return:
        """
        # downsample read length dist
        approx_ccl_ds = approx_ccl // window
        mult = np.arange(0.05, 1, 0.1)[::-1]
        # temporary container
        tmp_benefit = np.zeros(shape=(self.scores_ds.shape[0], 2))
        for i in range(10):
            b_part_fwd = bn.move_sum(self.scores_ds[::-1], window=int(approx_ccl_ds[i]), min_count=1)[::-1]
            b_part_rev = bn.move_sum(self.scores_ds, window=int(approx_ccl_ds[i]), min_count=1)
            # apply weighting by length
            wgt = mult[i]
            tmp_benefit[:, 0] += (b_part_fwd * wgt)
            tmp_benefit[:, 1] += (b_part_rev * wgt)
        # assign as attribute
        self.expected_benefit = tmp_benefit
        self.additional_benefit = self.expected_benefit - self.smu
        # correct numerical issues
        self.additional_benefit[self.additional_benefit < 0] = 0




class Reference:

    def __init__(self, ref: str, mmi: str = None, reject_refs: str = None):
        """
        Initialise a reference object. Loads contigs, load or create index
        and set contigs from which to always reject

        :param ref: Path to the reference fasta file
        :param mmi: Optional path to minimap index file
        :param reject_refs: Optional comma-sep list of headers in fasta
        """
        self.ref = ref
        self.mmi = mmi
        if not Path(ref).is_file():
            raise FileNotFoundError("Reference file not found")
        ref_suff = Path(ref).suffixes
        if not any([r in {".fa", ".fasta"} for r in ref_suff]):
            raise ValueError("Reference needs to be either fa, fasta")
        # create index if not present
        if not self.mmi:
            logging.info(f"Indexing reference: {self.ref}")
            mmi = f"{self.ref}.mmi"
            Indexer(fasta=self.ref, mmi=mmi)
            self.mmi = mmi
        else:
            if not Path(self.mmi).is_file():
                raise FileNotFoundError("Given mmi file not found")


        # headers of reference sequences from which to always reject
        if reject_refs:
            self.reject_refs = set(reject_refs.split(','))
        else:
            self.reject_refs = set()

        # load the sequences and their lengths
        logging.info("Reading reference file")
        self.contigs = self._load_contigs()
        self.n_sites = self._total_sites()



    def _load_contigs(self, min_len: int = 1e5, ploidy: int = 1) -> dict[str, Contig]:
        """
        Load contigs from fasta file

        :param min_len: Contigs shorter than this are ignored
        :param ploidy: Indicate haploid or diploid organism
        :return: dictionary of Contig objects
        """
        contigs = {}
        with open(self.ref, 'r') as fasta:
            for cname, cseq in read_fa(fasta):
                # drop >
                cname = cname[1:]
                # filter short reference sequences that can cause issues
                if len(cseq) < min_len:
                    continue
                # load reference sequences
                if cname not in self.reject_refs:
                    cseq_upper = cseq.upper()
                    contigs[cname] = Contig(name=cname, seq=cseq_upper, ploidy=ploidy)
                # for ref seqs that we always reject, set sequence empty
                else:
                    contigs[cname] = Contig(name=cname, seq="ACGT", ploidy=ploidy, rej=True)
        return contigs




    def _total_sites(self) -> int:
        """
        :return: Number of total sites in all contigs
        """
        return np.sum(list(self.contig_lengths().values()))


    def contig_sequences(self) -> dict[str, str]:
        """
        :return: Dictionary of contig sequences
        """
        return {c.name: c.seq for c in self.contigs.values()}


    def contig_lengths(self) -> dict[str, int]:
        """
        :return: Dictionary of contig names and their lengths
        """
        return {c.name: c.length for c in self.contigs.values()}



    def get_strategy_dict(self) -> dict[str, NDArray]:
        """
        Compile a dictionary of strategies from all contigs
        :return:
        """
        strat_dict = {}
        for cname, cont in self.contigs.items():
            strat_dict[cname] = cont.strat
        return strat_dict






