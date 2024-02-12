import logging
from collections import defaultdict

import numpy as np
from numpy.typing import NDArray

from boss.core import Boss
from boss.mapper import Mapper
from boss.utils import adjust_length

from boss.runs.sequences import CoverageConverter, Scoring
from boss.runs.reference import Reference
from boss.runs.abundance_tracker import AbundanceTracker
from boss.runs.readstartdist import ReadStartDist




class BossRuns(Boss):


    def init(self) -> None:
        """
        Initialise necessary components for a BOSS-RUNS experiment
        These components are shared across the live and simulation experiments,
        i.e. this init function is also run in the init of the simulation

        :return:
        """
        # initialise reference
        self.ref = Reference(ref=self.args.ref, reject_refs=self.args.reject_refs)
        self.contigs = self.ref.contigs
        self.contigs_filt = {n: c for n, c in self.contigs.items() if not c.rej}
        # initialise a mapper using the reference
        self.mapper = Mapper(ref=self.ref.mmi)
        # initialise a translator for coverage conversions
        self.cc = CoverageConverter()
        # initialise the abundance tracker
        self.tracker = AbundanceTracker(contigs=self.contigs)
        # initialise the tracker for read start position distribution
        self.read_starts = ReadStartDist(contigs=self.contigs_filt)
        # initialise scoring array
        self.scoring = Scoring(ploidy=self.args.ploidy)
        self.scoring.init_score_array()
        # write initial strategies to file
        strat_dict = self.ref.get_strategy_dict()
        self._write_contig_strategies(contig_strats=strat_dict)



    def _write_contig_strategies(self, contig_strats: dict[str, NDArray]) -> None:
        """
        Write the strategies for all contigs to a single file.

        :param contig_strats: A dictionary containing the strategies for contigs.
        """
        cpath = f'{self.out_dir}/masks/boss'
        np.savez(cpath, **contig_strats)
        # Example how to load these:
        # container = np.load(f'{cpath}.npz')
        # data = {key: container[key] for key in container}



    def _effect_increments(self, increments: defaultdict) -> None:
        """
        Loop through the increments of each contig and add the new coverage
        :param increments: defaultdict of lists of coverage counts per contig
        :return:
        """
        for cname, cont in self.contigs_filt.items():
            # grab the list of increments for this contig
            inc_list = increments[cname]
            cont.increment_coverage(increment_list=inc_list)



    def _update_scores_contigs(self) -> None:
        """
        Update the scores for each contig
        :return:
        """
        for cname, cont in self.contigs_filt.items():
            # main score updating function
            cont.scores, cont.entropy = self.scoring.update_scores(contig=cont)
            # modify scores after updating
            cont.modify_scores()


    def _check_buckets_contigs(self) -> bool:
        """
        Check if the buckets within the contigs have enough coverage to be switched on
        :return: Boolean if any strategy has been switched on
        """
        for cname, cont in self.contigs_filt.items():
            cont.check_buckets()
        # check if any switches are on
        switched_on = [c.switched_on for c in self.contigs.values()]
        return any(switched_on)


    def _update_benefits(self) -> None:
        """
        Update the expected benefit for all contigs
        :return:
        """
        for cname, cont in self.contigs_filt.items():
            cont.calc_smu()
            cont.calc_u(approx_ccl=self.rl_dist.approx_ccl)



    def _distribute_strategy(self, strat: NDArray, window: int = 100) -> None:
        """
        Place new decision strategies into the contigs strategies

        :param strat: Merged array of the updated strategy
        :param window: Downsampling window size
        :return:
        """
        i = 0
        for cname, cont in self.contigs_filt.items():
            # get the buckets of this contig and expand
            expand_fac = cont.bucket_size // window
            buckets_exp = np.repeat(cont.bucket_switches, expand_fac, axis=0)
            buckets = adjust_length(original_size=cont.strat.shape[0], expanded=buckets_exp)
            assert buckets.shape[0] == cont.strat.shape[0]
            # grab the new strategy
            cstrat = strat[i: i + cont.length // window, :]
            assert cstrat.shape == cont.strat.shape
            # assign new strat
            cont.strat[buckets, :] = cstrat[buckets, :]
            # log number of accepted sites
            f_perc = np.count_nonzero(cont.strat[:, 0]) / cont.strat.shape[0]
            r_perc = np.count_nonzero(cont.strat[:, 1]) / cont.strat.shape[0]
            logging.info(f'{cname}: {f_perc}, {r_perc}')
            i += cont.length // window


    def cleanup(self) -> None:
        """
        To be implemented

        :return:
        """
        pass


    def update_wrapper(self) -> None:
        """
        Second part of updates after a new batch of data.
        This is run after counting the coverage from new data

        :return:
        """
        # update the scores of all contigs
        self._update_scores_contigs()
        # flip strategy switches if threshold is reached
        switched_on = self._check_buckets_contigs()
        # UPDATE STRATEGY
        if switched_on:
            # update Fhat: unpack and normalise
            fhat_exp = self.read_starts.update_f_pointmass()
            self._update_benefits()
            # merge the benefits into one array for combined calculation
            benefit, smu = self.scoring.merge_benefit(self.contigs_filt)
            target_size = self.ref.n_sites // 100
            benefit_adj = adjust_length(original_size=target_size,
                                        expanded=benefit)
            smu_adj = adjust_length(original_size=target_size,
                                        expanded=benefit)
            assert fhat_exp.shape == benefit_adj.shape == smu_adj.shape
            # find the current decision strategy
            strat, threshold = self.scoring.find_strat_thread(
                benefit=benefit_adj,
                smu=smu_adj,
                fhat=fhat_exp,
                time_cost=self.rl_dist.time_cost
            )
            # distribute the strategy to the contigs
            self._distribute_strategy(strat=strat)
            # write strategies to file
            strat_dict = self.ref.get_strategy_dict()
            self._write_contig_strategies(contig_strats=strat_dict)



    def process_batch_runs(self, new_reads: dict[str, str], new_quals: dict[str, str]) -> None:
        """
        Process a batch of new data for the BOSS-RUNS mode
        This function is to be passed into the process_batch of the superclass

        :param new_reads: Dictionary of new sequences
        :param new_quals: Dictionary of qualities of new data
        :return:
        """
        # map the new reads to the reference
        paf_dict = self.mapper.map_sequences(sequences=new_reads)
        # convert coverage counts to increment arrays
        increments = self.cc.convert_records(paf_dict=paf_dict, seqs=new_reads, quals=new_quals)
        # effect the coverage increments for each contig
        self._effect_increments(increments=increments)
        # update the abundance tracker with new data
        self.tracker.update(n=len(new_reads), paf_dict=paf_dict)
        # update read starting position distribution
        self.read_starts.count_read_starts(paf_dict=paf_dict)
        # note: read length dist is updated in _get_new_data() of superclass
        self.update_wrapper()

