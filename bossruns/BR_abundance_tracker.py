from collections import defaultdict
import logging

# non std-library
import numpy as np
from natsort import natsorted

# custom imports
from .BR_batch import CurrentBatch



class AbundanceTracker:

    def __init__(self):
        self.total_reads = 0
        self.read_counts = defaultdict(int)


    def increment_observed_reads(self, n):
        self.total_reads += n


    def count_read_targets(self, cigar_dict):
        # in case there was no mapped read in the batch
        if len(cigar_dict) == 0:
            return

        # loop through reads to record lengths
        for read_id, records in cigar_dict.items():
            # select best mapper if there are multiple
            if len(records) > 1:
                records = CurrentBatch._choose_best_mapper(records)
            # grab the target of the mapping
            target = records[0].target_name
            self.read_counts[target] += 1


    def report_proportions(self):
        # calculate proportions of reads per target
        props = {tname: (n / self.total_reads) for tname, n in self.read_counts.items()}
        tnames_sorted = natsorted(list(props.keys()))
        # print to logfile
        logging.info("Counts and rel. proportions of observed reads:")
        for t in tnames_sorted:
            logging.info(f"{t}: {self.read_counts[t]} {np.round(props[t], 3)}")


    def update(self, n, cigar_dict):
        # at each update:
        # - increment the observed read count
        # - count the targets of observed reads
        # - report current numbers
        self.increment_observed_reads(n)
        self.count_read_targets(cigar_dict)
        self.report_proportions()




