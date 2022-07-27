import logging

# non-std library
import mappy
from natsort import natsorted
import numpy as np
# custom import
from BR_utils import integer_genome, read_fa


class ReferenceWGS:

    def __init__(self, ref, mmi=None, min_len=None):
        self.ref = ref
        self.mmi = mmi

        if not mmi:
            logging.info(f"indexing reference: {self.ref}")
            self.mmi = self.index_reference()

        # load the sequences and their lengths
        min_len = min_len if min_len else 1e5
        logging.info("reading reference file")
        seq_dict, len_dict, ind_dict = self.chromosome_dict(min_len)
        self.chromosome_lengths = len_dict
        self.chromosome_indices = ind_dict

        # transform chromosomes into integer array
        self.genome = integer_genome(chromosome_sequences=seq_dict)
        # number of sites
        self.roi_length = self.genome.shape[0]

        # create necessary arrays
        logging.info("generating array maps for reference")
        genome2roi_arr, chrom_rpos, chrom_rpos_ranges, roi_indices = self.gen_array_maps(len_dict)
        self.genome2roi_arr = genome2roi_arr
        self.chrom_rpos = chrom_rpos
        self.chrom_rpos_ranges = chrom_rpos_ranges
        self.roi_indices = roi_indices

        chNum, chArray, ch_cum_sum = self.gen_chrom_numbers(len_dict)
        self.chNum = chNum
        self.chArray = chArray
        self.ch_cum_sum = ch_cum_sum

    def index_reference(self):
        # name for index file
        idx = f"{self.ref}.mmi"
        # calling mappy with fn_idx_out writes index to file
        _ = mappy.Aligner(self.ref, preset="map-ont", n_threads=24, fn_idx_out=idx)
        return idx

    def chromosome_dict(self, min_len):
        chrom_dict = {}
        with open(self.ref, 'r') as fasta:
            for cname, cseq in read_fa(fasta):
                if len(cseq) > min_len:
                    chrom_dict[cname] = cseq.upper()
        # apply natsorting
        cnames_sorted = natsorted(list(chrom_dict.keys()))
        chrom_dict_sort = {cname: chrom_dict[cname] for cname in cnames_sorted}
        chrom_lengths_sort = {cname: len(chrom_dict[cname]) for cname in cnames_sorted}
        chrom_indices = {cname: i for i, cname in enumerate(cnames_sorted)}
        return chrom_dict_sort, chrom_lengths_sort, chrom_indices

    def gen_array_maps(self, len_dict):
        # dict that holds an array of chromosome length
        genome2roi_arr = {cname: np.arange(clen) for cname, clen in len_dict.items()}
        n = 0
        # dict that maps chromosomes to the concatenated array
        chrom_rpos = {}
        chrom_rpos_ranges = {}
        for cname, carr in genome2roi_arr.items():
            chrom_rpos[cname] = np.arange(n, n + carr[-1] + 1)
            chrom_rpos_ranges[cname] = (n, n + carr[-1])
            n += carr[-1] + 1
        roi_indices = np.arange(0, n, dtype="uint64")
        return genome2roi_arr, chrom_rpos, chrom_rpos_ranges, roi_indices

    def gen_chrom_numbers(self, len_dict):
        # number of chromosomes
        chNum = len(len_dict.keys())
        # array of chrom lengths
        chArray = np.array(list(len_dict.values()))
        # cumulative sum of chromosome lengths, starting with 0
        ch_cum_sum = np.append(np.zeros(1, dtype='uint64'),
                               np.cumsum(chArray, dtype='uint64'))
        return chNum, chArray, ch_cum_sum



