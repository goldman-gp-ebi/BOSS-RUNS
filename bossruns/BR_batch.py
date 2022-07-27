import re
import logging
from pathlib import Path
import gzip
from collections import defaultdict
from io import StringIO
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor as TPexe

# custom imports
from .BR_utils import readfq, range_intersection
from .BR_cov import convert_records_helper
from .BR_sample_fq import parallel_batches

# non-std lib
import numpy as np





class CurrentBatch:


    def __init__(self, num, qt, window):
        self.num = num
        self.qt = qt  # base quality threshold, set as constant
        self.window = window
        # set up translation dict for bases
        transDict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.base2int = str.maketrans(transDict)
        # translation dict for qualities
        transDict_qual = {chr(min(126, q + 33)): q for q in range(94)}  # translator for phred encoding
        self.qual2int = str.maketrans(transDict_qual)
        # and a translation dict for the cigar ops
        transDict_cig = {'M': 0, 'D': 1, 'I': 2, 'S': 3}
        self.cig2int = str.maketrans(transDict_cig)
        # compile regex for cigar ops
        self.cigar_regex = re.compile("(\d+)([MIDNSHP=XB])")


    def read_batch(self, fastq_file, channels=None):
        '''
        method to ingest one new fastq file. Fills containers of read ids, lengths, sequences, etc.

        Parameters
        ----------
        fastq_file: str
            path to the next fastq file to read
        channels: set
            channel numbers assigned to the bossruns condition on the flowcell
            only necessary if there are multiple conditions in the experiment

        Returns
        -------

        '''
        read_lengths = {}
        read_sequences = {}
        read_qualities = {}
        basesTOTAL = 0
        n_reads = 0

        logging.info(f"reading file: {fastq_file}")

        # to make sure its a path object, not a string
        if type(fastq_file) is str:
            fastq_file = Path(fastq_file)

        # check whether fastq is gzipped
        if fastq_file.name.endswith(('.gz', '.gzip')):
            fh = gzip.open(fastq_file, 'rt')
        else:
            fh = open(fastq_file, 'rt')

        # loop over all reads in the fastq file
        if not channels:
            # if we don't do any filtering by channel numbers
            for desc, name, seq, qual in readfq(fh):
                bases_in_read = len(seq)
                read_lengths[str(name)] = bases_in_read
                read_sequences[str(name)] = seq
                read_qualities[str(name)] = qual
                basesTOTAL += bases_in_read
                n_reads += 1
        else:
            # if we filter the incoming batch by the channel that the read comes from
            for desc, name, seq, qual in readfq(fh):
                # find the source channel
                try:
                    # regex to get the channel number from the header
                    # \s=whitespace followed by 'ch=' and then any amount of numeric characters
                    curr_channel = re.search("\sch=[0-9]*", desc).group()
                    ch_num = int(curr_channel.split('=')[1])
                except AttributeError:
                    # if the pattern is not in the header, skip the read
                    logging.info("ch= not found in header of fastq read")
                    continue

                if ch_num in channels:
                    # check if the read comes from a channel that is in the set of bossruns channels
                    bases_in_read = len(seq)
                    read_lengths[str(name)] = bases_in_read
                    read_sequences[str(name)] = seq
                    read_qualities[str(name)] = qual
                    basesTOTAL += bases_in_read
                    n_reads += 1
        fh.close()

        logging.info(f"processing {n_reads} reads in this batch")
        # assign attributes used in mapping and makeDecision
        read_ids = set(read_sequences.keys())
        self.read_ids = read_ids
        self.read_lengths = read_lengths
        self.read_sequences = read_sequences
        self.read_qualities = read_qualities
        self.basesTOTAL = basesTOTAL
        # also return them for the case of processing multiple fastq files
        return read_ids, read_lengths, read_sequences, read_qualities, basesTOTAL


    def read_batch_mmap(self, fq_stream):
        '''
        method to get new sequencing reads from a memory mapped fastq file
        used for simulations onlt

        Parameters
        ----------
        fq_stream: FastqStream
            class representing very large memory-mapped fastqfile

        Returns
        -------

        '''
        # fill dicts read_lengths and read_sequences
        read_lengths, read_sequences, basesTOTAL = fq_stream.get_batch()
        self.read_ids = set(read_sequences.keys())
        self.read_lengths = read_lengths
        self.read_sequences = read_sequences
        self.basesTOTAL = basesTOTAL


    def read_batch_parallel(self, fq, bsize):
        '''
        method to get new sequencing reads from a very large fastq file
        uses multiprocessing and random access. Only useful if files in flash storage

        Parameters
        ----------
        fq: str
            path to fastq file
        bsize: int
            number of reads to get

        Returns
        -------

        '''
        # fill the dicts read_lengths and read_sequences
        read_lengths, read_sequences, basesTOTAL = parallel_batches(fq, bsize, workers=10)
        self.read_ids = set(read_sequences.keys())
        self.read_lengths = read_lengths
        self.read_sequences = read_sequences
        self.basesTOTAL = basesTOTAL


    def make_decision(self, paf_output, otu, mode):
        '''
        SIMULATION ONLY:
        Makes decisions about which reads to accept and reject in a PAF file

        Parameters
        ----------
        paf_output: str
            alignments of reads to reference from minimap2
        otu: OTU
            object of BR
        mode: str
            decisions are made for this condition

        Returns
        -------

        '''
        # transform from PAF to dictionary of paf records
        # requires finding mapping pos in respect to full genome
        # init cigar dicts (read_ids as keys)
        # They have a list per read_id for multimappers
        # The list elements are PAF instances from parsing the PAF
        # dict -> list -> PAF record
        cigar_dict_accept = defaultdict(list)
        cigar_dict_reject = defaultdict(list)
        mapped_reads = set()
        count_mappers = 0

        # grab strategy
        strat = getattr(otu, f'strat_dict_{mode}')

        # yield each mapping as a named tuple that we append to the cigar dict directly
        for record in PAF.parse_PAF(StringIO(paf_output)):
            # record = list(parse_PAF(StringIO(paf_output)))[n] # DEBUGGING

            # increment mapping counter and add to set
            count_mappers += 1
            mapped_reads.add(record.query_name)

            # find start and end position relative to whole genome
            if record.strand == '+':
                record.genome_start = record.target_start
                record.genome_end = record.target_end - 1

            elif record.strand == '-':
                record.genome_start = record.target_end - 1
                record.genome_end = record.target_start
            else:
                continue

            # decision process and creation of cigar dict
            # index dictionary of strategies with otu,
            # then index target chromosome
            # then index starting position and strand
            try:
                decision = strat[str(record.target_name)][record.genome_start // self.window][record.strand_index]

            except (KeyError, IndexError):
                # in case the read maps to a chromosome that we don't have a strategy for
                # can happen in testing or for references with scaffolds
                # these are rejected by default (except in control condition of sims)
                if mode == 'NV':
                    decision = 1
                else:
                    decision = 0

            # ACCEPT READ
            if decision:
                cigar_dict_accept[str(record.query_name)].append(record)

            # REJECT READ
            else:
                cigar_dict_reject[str(record.query_name)].append(record)

        unmapped_reads_ID = self.read_ids - mapped_reads
        unmapped_reads = len(unmapped_reads_ID)
        return cigar_dict_accept, cigar_dict_reject, unmapped_reads, unmapped_reads_ID


    def convert_paf(self, paf_raw):
        '''
        convert the output of mapping the full reads from the sequencer
        inorder to parse the coverage counts at each site
        shortened version of make_decision used in simulations

        Parameters
        ----------
        paf_raw: str
            alignments of reads to reference from minimap2 in PAF format

        Returns
        -------
        cigar_dict: dict
            dict of read_ids: list of parsed alignments

        '''
        # init cigar dict
        # list per read_id for multimappers
        # list elements are PAF records
        # dict -> list -> PAF record
        cigar_dict = defaultdict(list)

        # yield each mapping as a named tuple that is appended to the cigar dict
        for record in PAF.parse_PAF(StringIO(paf_raw)):
            # record = list(parse_PAF(StringIO(paf_output)))[n] # DEBUGGING

            # find start and end position
            if record.strand == '+':
                record.genome_start = record.target_start
                record.genome_end = record.target_end - 1

            elif record.strand == '-':
                record.genome_start = record.target_end - 1
                record.genome_end = record.target_start
            else:
                continue

            # transfer the read to the cigar dict
            cigar_dict[str(record.query_name)].append(record)
        return cigar_dict


    def convert_coverage(self, cigars_accept, cigars_reject, cigars_full, otu, workers, whole_genome):
        '''
        Takes the dictionary of PAF entries from mapped reads and converts to coverage counts
        cigars_reject is only used for simluations, in live application all mappings are in cigars_full

        Parameters
        ----------
        cigars_accept: dict
            read ids: list of PAF records
        cigars_reject: dict
            read ids: list of PAF records
        cigars_full: dict
            read ids: list of PAF records
        otu: OTU
            br specific object
        workers: int
            number of threads for parallelisation
        whole_genome: bool
            flag to specify parsing algorithm

        Returns
        -------
        coverage_counts: np.ndarray
            Coverage counts in 5 by N array.
        bases: int
            the number of total sequenced bases, used during simulations to increment time
        '''
        # empty array to hold new counts at each position
        coverage_counts = np.zeros((otu.roi_length, 5), dtype=np.uint16)
        bases = 0

        accepted_reads = set(cigars_accept.keys())
        rejected_reads = set(cigars_reject.keys())

        # count the reads that were mapped in full but not in trunc mode
        count_trunc_non_mappers = 0
        # add up the coverage from all cigar strings
        read_ids = list(cigars_full.keys())
        # collect either the full or trunc mappings for each read
        record_list = []

        for read_id in read_ids:
            # read_id = list(cigars_full.keys())[19]  # DEBUGGING
            # if the read was accepted we use the full version
            # this is the case in live application
            if read_id in accepted_reads:
                records = cigars_full[read_id]
                # SIM: add the length of the sequence to bases, i.e. time count
                # this is for RU and BR - depending on which strategy is used
                bases += len(self.read_sequences[read_id])

            elif read_id in rejected_reads:
                records = cigars_reject[read_id]
            else:
                count_trunc_non_mappers += 1
                continue

            # select best mapper or multiple if justified
            if len(records) > 1:
                records = CurrentBatch._choose_best_mapper(records)

            record_list.append(records)

        # setup for multi-threading
        # conv is a container for everything necessary for cigar parsing
        conv = SimpleNamespace(base2int=self.base2int, qual2int=self.qual2int,
                               cig2int=self.cig2int, cigar_regex=self.cigar_regex,
                               genome2roi=otu.genome2roi, genome2roi_arr=otu.genome2roi_arr,
                               read_seqs=self.read_sequences, read_quals=self.read_qualities,
                               chrom_rpos=otu.chrom_rpos)

        # separate records into chunks for each worker
        record_arrays = np.array_split(record_list, workers)

        # make argument list
        arg_list = []
        for r in range(len(record_arrays)):
            arg_list.append((record_arrays[r], conv, self.qt, whole_genome))

        # ri = convert_records_helper(arguments=arg_list[0]) # DEBUG

        # multithread magic
        with TPexe(max_workers=workers) as executor:
            roi_increments = executor.map(convert_records_helper, arg_list)

        # combine the results from all workers
        inc_list = [x for x in list(roi_increments) if x is not None]

        if inc_list:
            increment_stack = np.concatenate(inc_list)
            # add to the coverage counts - add.at takes 2d indices as tuple
            np.add.at(coverage_counts, (increment_stack[:, 0], increment_stack[:, 1]), 1)
        return coverage_counts, bases


    @staticmethod
    def _choose_best_mapper(records):
        # structured array to decide between ties, by using the score of the DP algorithm
        mapq = [(record.mapq, record.align_score) for record in records]
        custom_dtypes = [('q', int), ('dp', int)]
        mapping_qualities = np.array(mapq, dtype=custom_dtypes)
        sorted_qual = np.argsort(mapping_qualities, order=["q", "dp"])
        record = [records[sorted_qual[-1]]]
        return record


    @staticmethod
    def choose_multimapper(records):
        # safety check that no large amount of different records is present
        if len(records) < 20:
            # get starts and ends of each mapping
            t_ranges = [range(r.target_start, r.target_end) for r in records]
            q_ranges = [range(r.query_start, r.query_end) for r in records]

            # check for overlaps between the different mappings
            t_overlaps = [range_intersection(x, y) for i, x in enumerate(t_ranges)
                          for j, y in enumerate(t_ranges) if i > j]

            q_overlaps = [range_intersection(x, y) for i, x in enumerate(q_ranges)
                          for j, y in enumerate(q_ranges) if i > j]

            # if there are no mappings that overlap more than 10bp on the query or target: accept them all
            if not any([q > 10 for q in q_overlaps]) and not any([t > 10 for t in t_overlaps]):
                return records

            else:
                # if there are too many overlaps, just use the one with the best mapping
                record = CurrentBatch._choose_best_mapper(records)
        else:
            # if there are many little mappings, use the best one instead
            record = CurrentBatch._choose_best_mapper(records)
        return record



    # def write_seqs(self, out_dir, mode):
    #     # write all sequences to fastq file
    #     fname = f'{out_dir}/fq/{mode}/reads.fq'
    #     with open(fname, 'w') as reads:
    #         for rid, seq in self.read_sequences.items():
    #             reads.write(f'>{rid}\n')
    #             reads.write(f'{seq}\n')
    #     return fname
    #
    # def write_seqs_trunc(self, accepted, mu, out_dir, mode):
    #     # write sequences, taking rejections into account
    #     acc_BR = set(accepted.keys())
    #     fname = f'{out_dir}/fq/{mode}/reads.fq'
    #     with open(fname, 'w') as reads:
    #         for rid, seq in self.read_sequences.items():
    #             reads.write(f'>{rid}\n')
    #
    #             if rid in acc_BR:
    #                 reads.write(f'{seq}\n')
    #             else:
    #                 reads.write(f'{seq[: mu]}\n')
    #     return fname





class PAF:

    def __init__(self, query_name, query_len, query_start, query_end, strand, target_name, target_len, target_start,
                 target_end, num_matches, alignment_block_length, mapq):
        self.query_name = query_name
        self.query_len = query_len
        self.query_start = query_start
        self.query_end = query_end
        self.strand = strand
        # to index into the forward/reverse strategy
        self.strand_index = 0 if self.strand == "+" else 1
        self.target_name = target_name
        self.target_len = target_len
        self.target_start = target_start
        self.target_end = target_end
        self.num_matches = num_matches
        self.alignment_block_length = alignment_block_length
        self.mapq = mapq
        self.map_len = self.target_end - self.target_start


    @staticmethod
    def parse_PAF(file_like, fields=None):
        """
        Read a minimap2 PAF file as an iterator
        Parameters
        ----------
        file_like : file-like object
            Object with a read() method, such as a file handler or io.StringIO.
        fields : list
            List of field names to use for records, must have 13 entries.
            Default: ["read_name", "query_length", "query_start", "query_end", "strand",
             "target_name", "target_length", "target_start", "target_end", "residue_matches",
              "alignment_block_length", "mapping_quality", "tags"]
        Returns
        -------
        iterator
        """
        if fields is None:
            fields = [
                "query_name",
                "query_len",
                "query_start",
                "query_end",
                "strand",
                "target_name",
                "target_len",
                "target_start",
                "target_end",
                "num_matches",
                "alignment_block_length",
                "mapq",
                "tags",
            ]

        return PAF._paf_generator(file_like, fields=fields)


    @staticmethod
    def _paf_generator(file_like, fields=None):
        """Generator that returns records from a PAF file

        Parameters
        ----------
        file_like : file-like object
            File-like object
        fields : list
            List of field names to use for records, must have 13 entries.

        Yields
        ------
        record
            Correctly formatted PAF record and a dict of extra tags

        Raises
        ------
        ValueError
        """
        if len(fields) != 13:
            raise ValueError(f"{len(fields)} fields provided, expected 13")
        # PAF = namedtuple("PAF", fields)
        for record in file_like:
            record = record.strip().split("\t")
            # new version that uses mutable object rather than namedtuple
            f = PAF._format_records(record[:12])
            paf = PAF(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11])
            tags = PAF._parse_tags(record[12:])
            # paf.tags = tags  # we don't need to keep all the tags, just AS and cigar for now
            paf.align_score = int(tags.get("AS", 0))
            paf.cigar = tags.get("cg", None)
            yield paf


    @staticmethod
    def _format_records(record):
        # Helper function to convert fields to the correct type
        return [PAF._conv_type(x, int) for x in record]


    @staticmethod
    def _conv_type(s, func):
        """
        Generic converter to convert strings to other types

        Parameters
        ----------
        s : str
            string that represents another type
        func : type
            type to apply to s

        Returns
        -------
        The type of func, otherwise str
        """
        try:
            return func(s)
        except ValueError:
            return s


    @staticmethod
    def _parse_tags(tags):
        """
        Convert list of SAM style tags from a PAF file to a dict

        Parameters
        ----------
        tags : list
            list of SAM style tags

        Returns
        -------
        dict
            dict of SAM style tags
        """
        c = {"i": int, "A": str, "f": float, "Z": str}
        return {
            key: PAF._conv_type(val, c[tag])
            for key, tag, val in (x.split(":") for x in tags)
        }








