import logging
from collections import defaultdict
from io import StringIO
from copy import deepcopy

from boss.runs.core import BossRuns
from boss.sampler import Sampler
from boss.paf import Paf, PafLine, paf_dict_type
from boss.batch import ReadCache


class BossRunsSim(BossRuns):

    def init_sim(self):
        # run the init from super
        self.init()

        # initialise the sampler
        self.sampler = Sampler(
            source=self.args.fq,
            paf_full=self.args.paf_full,
            paf_trunc=self.args.paf_trunc,
            maxbatch=self.args.maxb,
            batchsize=self.args.batchsize
        )
        # initialise pseudotiming object
        self.read_cache = ReadCache(batchsize=self.args.batchsize, dumptime=self.args.dumptime)
        self.mu = 400



    def make_decisions(
        self,
        seqs: dict[str, str],
        paf_full: str,
        paf_trunc: str,
        window: int = 100
    ) -> tuple[paf_dict_type, dict[str, str], int, int, int, int]:
        """
        Make decisions about sampled data and build the paf_dict

        :param seqs: Dict of raw sequences
        :param paf_full: Raw output of mapping full-length reads
        :param paf_trunc: Raw output of mapping truncated reads
        :param window: downsampling size
        :return: pad_dict, and numbers of unmapped and rejected reads
        """
        # build a paf dict and either accept or reject reads
        paf_dict = defaultdict(list)
        mapped_reads = set()
        n_rejected = 0
        n_accepted = 0
        # return a dict of sequences after decision process
        reads_decision = deepcopy(seqs)

        # parse the raw paf data into dictionaries
        paf_dict_full = Paf.parse_PAF(StringIO(paf_full))
        paf_dict_trunc = Paf.parse_PAF(StringIO(paf_trunc))

        # loop over the mu-sized mappings to make decisions
        for rid, rlist in paf_dict_trunc.items():
            rec = Paf.choose_best_mapper(rlist)[0]
            mapped_reads.add(rid)
            # deal with strandedness
            if rec.rev:
                start_pos = rec.tend - 1
            else:
                start_pos = rec.tstart
            # actual decision look-up
            try:
                strat = self.contigs_filt[str(rec.tname)].strat
                decision = strat[start_pos // window, rec.rev]

            except (KeyError, IndexError):
                # in case the read maps to a chromosome that we don't have a strategy for
                # reject by default (can happen in testing or for references with short scaffolds)
                decision = 0

            if decision:
                # ACCEPT READ
                # grab the full-length version
                rec_full_list = paf_dict_full[str(rec.qname)]
                rec_full = Paf.choose_best_mapper(rec_full_list)[0]
                paf_dict[str(rec.qname)].append(rec_full)
                n_accepted += 1
            else:
                # REJECT READ
                paf_dict[str(rec.qname)].append(rec)
                n_rejected += 1
                reads_decision[rid] = reads_decision[rid][self.mu:]


        n_mapped = len(mapped_reads)
        n_unmapped = len(self.sampler.fq_stream.read_ids - mapped_reads)
        return paf_dict, reads_decision, n_mapped, n_unmapped, n_accepted, n_rejected



    def filter_paf_dict(self, paf_dict: dict[str, list[PafLine]]) -> dict[str, list[PafLine]]:
        """
        Filter the paf_dict down to only the reads that were accepted
        :param paf_dict:
        :return:
        """
        paf_dict_acc = {}
        for rid, reclist in paf_dict.items():
            rec = reclist[0]
            if rec.qlen != self.read_cache.mu:
                paf_dict_acc[rid] = reclist
        return paf_dict_acc



    def process_batch_runs_sim(self) -> None:
        """
        Process a new batch of simulations, only the first part is different to live runs
        The second part is run inside the update_wrapper() of the superclass
        :return:
        """
        # trigger the sampling of reads and their mappings
        read_seqs, read_quals, paf_f, paf_t = self.sampler.sample()
        # make decisions and generate the paf_dict
        paf_dict, reads_decision, n_mapped, n_unmapped, n_accepted, n_rejected = (
            self.make_decisions(seqs=read_seqs,
                                paf_full=paf_f,
                                paf_trunc=paf_t
                                )
        )
        logging.info(f"mapped {n_mapped}, not mapped {n_unmapped}")
        logging.info(f"accepted {n_accepted}, rejected {n_rejected}")
        # filter paf_dict to accepted reads only
        paf_dict_acc = self.filter_paf_dict(paf_dict=paf_dict)
        # update read length distribution
        self.rl_dist.update(read_lengths={n: r[0].qlen for n, r in paf_dict_acc.items()})
        # NOTE from here similar but not identical to live version
        # convert coverage counts to increment arrays
        increments = self.cc.convert_records(paf_dict=paf_dict, seqs=read_seqs, quals=read_quals)
        # effect the coverage increments for each contig
        self._effect_increments(increments=increments)
        # update the abundance tracker with new data
        self.tracker.update(n=n_accepted, paf_dict=paf_dict_acc)
        # update read starting position distribution
        self.read_starts.count_read_starts(paf_dict=paf_dict_acc)
        # update pseudotimes and check if it's time to write cache to file
        self.read_cache.update_times_runs(
            total_bases=self.sampler.fq_stream.total_bases,
            paf_dict=paf_dict,
            n_unmapped=n_unmapped,
            n_reject=n_rejected
        )
        self.read_cache.fill_cache(read_sequences=self.sampler.fq_stream.read_sequences, reads_decision=reads_decision)
        # update wrapper of superclass
        self.update_wrapper()










