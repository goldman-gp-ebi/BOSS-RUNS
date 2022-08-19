import logging
from collections import defaultdict

# non-std library
import mappy
from natsort import natsorted
import numpy as np
# custom import
from .BR_utils import integer_genome, read_fa



class Reference:

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
        self.chromosome_sequences = seq_dict
        self.chromosome_lengths = len_dict
        self.chromosome_indices = ind_dict

        chNum, chArray, ch_cum_sum = self.gen_chrom_numbers()
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

    def gen_chrom_numbers(self):
        # number of chromosomes
        chNum = len(self.chromosome_lengths.keys())
        # array of chrom lengths
        chArray = np.array(list(self.chromosome_lengths.values()))
        # cumulative sum of chromosome lengths, starting with 0
        ch_cum_sum = np.append(np.zeros(1, dtype='uint64'),
                               np.cumsum(chArray, dtype='uint64'))
        return chNum, chArray, ch_cum_sum


class ReferenceWGS(Reference):


    def load_reference(self):
        # transform chromosomes into integer array
        self.genome = integer_genome(chromosome_sequences=self.chromosome_sequences)
        # number of sites
        self.roi_length = self.genome.shape[0]

        # create necessary arrays
        logging.info("generating array maps for reference")
        genome2roi_arr, chrom_rpos, chrom_rpos_ranges, roi_indices = self.gen_array_maps()
        self.genome2roi_arr = genome2roi_arr
        self.chrom_rpos = chrom_rpos
        self.chrom_rpos_ranges = chrom_rpos_ranges
        self.roi_indices = roi_indices
        self.genome2roi = 0  # dummy init


    def gen_array_maps(self):
        # dict that holds an array of chromosome length
        genome2roi_arr = {cname: np.arange(clen) for cname, clen in self.chromosome_lengths.items()}
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




class ReferenceVCF(Reference):

    def load_reference(self, vcf):
        # get all sites from the vcf file
        logging.info("reading vcf file")
        vcf_sites = self.read_vcf(vcf)
        # check chromosome oder
        logging.info("checking chromosome order in reference and vcf")
        self.check_vcf_order(vcf_sites=vcf_sites)
        # generate all necessary maps between chrom : gpos : rpos
        logging.info("generating array maps for reference")
        genome2roi = self.gen_genome2roi(vcf_sites)
        roi_length, genome, genome2roi_arr, chrom_rpos, chrom_rpos_ranges = \
            self.gen_array_maps(vcf_sites, genome2roi)
        self.genome = genome
        self.roi_length = roi_length
        self.genome2roi = genome2roi
        self.genome2roi_arr = genome2roi_arr
        self.chrom_rpos = chrom_rpos
        self.chrom_rpos_ranges = chrom_rpos_ranges

        # generate indices with offsets of chromosomes
        roi_indices = self.gen_roi_indices_offset(self.chromosome_indices, genome2roi, self.ch_cum_sum)
        self.roi_indices = roi_indices



    def read_vcf(self, vcf):
        # collect sites in dict
        sites = defaultdict(set)
        # i = 0  # limit of sites
        # lim = 100000
        with open(vcf, 'r') as vc:
            for line in vc:
                # skip empty
                if len(line) < 1:
                    continue
                # skip headers
                if line.startswith('#'):
                    continue
                # parse line
                ll = line.split('\t')

                chrom = ll[0].split(':')[-1]
                if chrom not in self.chromosome_indices.keys():
                    continue
                # this is sometimes necessary depending on vcf formatting
                # chrom = chrom.replace('chr', '')
                # position of site
                position = ll[1]
                # base at site
                # ref = ll[3]
                sites[chrom].add(int(position))
                # testing
                # if lim:
                #     i += 1
                #     if i > lim:
                #         break
        # convert sets of sites to arrays and sort
        pos = {k: np.sort(np.array(list(v))) for k, v in sites.items()}
        return pos


    def check_vcf_order(self, vcf_sites):
        # check that the order of variants in the file is the same as the chromosome dict
        chrom_names_ref = self.chromosome_sequences.keys()
        chrom_names_vcf = list(vcf_sites.keys())
        chrom_names_ref_filt = [c for c in chrom_names_vcf if c in chrom_names_ref]
        assert chrom_names_ref_filt == chrom_names_vcf


    def gen_genome2roi(self, vcf_sites):
        # coordinate translation dict
        r = 0
        genome2roi = defaultdict(dict)
        for cname, gpos_arr in vcf_sites.items():
            chrom_n_sites = gpos_arr.shape[0]
            # range of the positions within all of the ROIs
            crpos = np.arange(r, r + chrom_n_sites)
            site_map =  dict(zip(gpos_arr, crpos))
            r += chrom_n_sites
            genome2roi[cname] = site_map
        return genome2roi


    def gen_array_maps(self, vcf_sites, genome2roi):
        # total number of rois
        roi_length = np.sum([site_arr.shape[0] for site_arr in vcf_sites.values()])

        # save as array for each chrom (used in coverage conv)
        genome2roi_arr = dict()
        # coordinate container holds rpos for each chrom
        chrom_rpos = dict()
        # more efficient version: rois from chromosomes are always consecutive
        chrom_rpos_ranges = dict()

        for c, roi_dict in genome2roi.items():
            genome2roi_arr[c] = np.array(list(roi_dict.keys()))
            chrom_rois = list(roi_dict.values())
            chrom_rpos[c] = np.array(chrom_rois)
            chrom_rpos_ranges[c] = (chrom_rois[0], chrom_rois[-1])  # hsap

        # no reference bases needed for uniform priors
        genome = np.zeros(shape=roi_length, dtype="uint8")
        return roi_length, genome, genome2roi_arr, chrom_rpos, chrom_rpos_ranges


    def gen_roi_indices_offset(self, chromosome_indices, genome2roi, ch_cum_sum):
        # indices to locate ROIs (used for scores_placed)
        roi_indices = []
        for c, c_ind in chromosome_indices.items():
            chrom_pos = np.array(list(genome2roi[c].keys()))
            chrom_offset = ch_cum_sum[c_ind]
            roi_indices.extend(chrom_pos + chrom_offset)
        roi_indices = np.array(roi_indices, dtype="uint64")
        return roi_indices


