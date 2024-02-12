import logging

import numpy as np

from boss.paf import Paf, paf_dict_type



class AbundanceTracker:

    def __init__(self, contigs: dict):
        """
        Initialise the abundance tracker to count observed molecules

        :param contigs: Dictionary of contigs
        """
        self.total_reads = 0
        self.read_counts = {cname: 0 for cname in contigs.keys()}


    def _count_read_targets(self, paf_dict: paf_dict_type) -> None:
        """
        Count up the target sequences of the mappings in the paf_dict

        :param paf_dict: Dict of mappings
        :return:
        """
        # in case there was no mapped read in the batch
        if len(paf_dict) == 0:
            return

        for rid, rec in paf_dict.items():
            # select best mapper if there are multiple
            if len(rec) > 1:
                rec = Paf.choose_best_mapper(rec)
            # grab target of the mapping
            t = rec[0].tname
            self.read_counts[t] += 1



    def _report_proportions(self) -> None:
        """
        Calculate and report the proportions of reads per target

        :return:
        """
        props = {tname: (n / self.total_reads) for tname, n in self.read_counts.items()}
        logging.info("Counts and rel. proportions of observed reads:")
        for t in list(props.keys()):
            logging.info(f"{t}: {self.read_counts[t]} {np.round(props[t], 3)}")



    def update(self, n: int, paf_dict: paf_dict_type) -> None:
        """
        At each update we increment the observed read count,
        and count the targets of observed reads

        :param n: Total number of new reads
        :param paf_dict: Dict of mappings
        :return:
        """
        self.total_reads += n
        self._count_read_targets(paf_dict)
        self._report_proportions()




