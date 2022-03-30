from bisect import bisect_left, bisect_right
import logging

# non-std lib
import numpy as np
from numba import njit

# custom imports
from .BR_utils import reverse_complement


def convert_records_helper(arguments):
    '''
    helper function for multi-threading, where only a single arg can be passed
    this is the workload of a single worker, i.e. looping over records

    Parameters
    ----------
    arguments: list
        list of arguments passed to single worker

    Returns
    -------
    increments: np.array
        array of counts that indicate which base at which position should be incremented

    '''
    # unpack args
    record_list, conv, qt, whole_genome = arguments
    # record_list, conv, qt = arg_list[0]
    increment_list = []
    for records in record_list:
        # loop through all valid mappings
        increments = convert_records(
            records=records,
            conv=conv,
            qt=qt,
            whole_genome=whole_genome)
        increment_list.append(increments)

    increments_filt = [x for x in increment_list if x is not None]
    if len(increments_filt) > 0:
        return np.concatenate(increments_filt)
    else:
        return


def convert_records(records, conv, qt, whole_genome):
    '''
    wrapper for full parsing of paf records to coverage counts

    Parameters
    ----------
    records: list
        list of valid mappings from a single read
    conv: obj
        container with translation dicts
    qt: int
        minimum base quality threshold
    whole_genome: bool
        flag to indicate which parsing algorithm to use

    Returns
    -------
    increments: np.array
        which bases at which position should be incremented

    '''
    for record in records:
        # check if the read covers any sites of interest
        start = min(record.genome_start, record.genome_end)
        end = max(record.genome_start, record.genome_end)
        # intersect area that read covers and ROIs of chromosome
        target_chromosome_rois = conv.genome2roi[str(record.target_name)]
        try:
            chrom_roi_arr = conv.genome2roi_arr[str(record.target_name)]
        except KeyError:
            # we might have reads mapping to some chrom, but no rois
            return
        left_idx = bisect_left(chrom_roi_arr, start)
        right_idx = bisect_right(chrom_roi_arr, end)

        # indices within current chromosome
        covered_rois = chrom_roi_arr[left_idx: right_idx]

        # if read does not cover any sites of interest, return
        if len(covered_rois) == 0:
            return

        # it might happen that left is 0 and right is -1
        # then we would suddenly cover all rois
        if left_idx == 0 and right_idx in {0, -1}:
            return

        # otherwise parse cigar
        if record.strand == '-':
            seq = reverse_complement(conv.read_seqs[record.query_name])
            qual = conv.read_quals[record.query_name][::-1]
            offcut = len(seq) - record.query_end
        elif record.strand == '+':
            seq = conv.read_seqs[record.query_name]
            qual = conv.read_quals[record.query_name]
            offcut = record.query_start
        else:
            logging.info("non-standard strand symbol")
            break

        # parse cigar takes the cigar string, the sequence and the mapping length
        # returns an array of counts for each base and indels
        query_arr, qual_arr = _parse_cigar(cigar_string=record.cigar, read_string=seq, qual_string=qual,
                                          chunk_length=record.map_len, offcut=offcut,
                                          base2int=conv.base2int, qual2int=conv.qual2int,
                                          cig2int=conv.cig2int, cigar_regex=conv.cigar_regex)

        # transform into increments
        if whole_genome:
            increments = _collect_increments_wg(covered_rois=covered_rois,
                                                         target_chromosome_rois=target_chromosome_rois,
                                                         query_arr=query_arr)

        else:
            increments = _collect_increments(covered_rois=covered_rois,
                                            start=start,
                                            target_chromosome_rois=target_chromosome_rois,
                                            query_arr=query_arr,
                                            qual_arr=qual_arr,
                                            qt=qt)

        # if none of the covered rois had coverage of sufficient quality, just return
        if len(increments) == 0:
            return
        else:
            return increments


def _parse_cigar(cigar_string, read_string, qual_string, chunk_length, offcut, base2int, qual2int, cig2int, cigar_regex):
    """
    Return array of integer bases

    Parameters
    ----------
    cigar_string : str
        A CIGAR string
    read_string : str
        Read that the CIGAR string describes
    qual_string: str
        string of base qualities (phred)
    chunk_length : int
        Length of the reference covered by the alignment.
    offcut: int
        piece of the read at the start, which did not map. Offset for the query_pos
    base2int: str.maketrans
        translation of bases to integer
    qual2int: str.maketrans
        translation of qualities to integers
    cig2int: str.maketrans
        translation of cigar ops to integers
    cigar_regex: re.pattern
        regex patterns to parse cigar ops

    Returns
    -------
    query_array: list
        integerised bases of query at each pos
    qual_array: list
        quality scores of query at each pos

    """
    # translate read seq into integers
    # base2int is a translations table from str.maketrans
    read_integer = read_string.translate(base2int)
    intSeq = np.fromstring(read_integer, dtype='int8', sep='')
    # translate qual to integers
    qual_integer = qual_string.translate(qual2int)
    int_qual = np.fromstring(qual_integer, dtype='int8', sep='')

    # split up the cigar
    lengths, ops = _prep_cigar(cigar_string, cigar_regex, cig2int)

    # loop through all cigar operators
    query_array, qual_array = _loop_cigar(cigar_lengths=lengths, cigar_ops=ops,
                                         int_seq=intSeq, int_qual=int_qual,
                                         chunk_length=chunk_length, qpos=offcut)
    return query_array, qual_array



def _prep_cigar(cigar_string, cigar_regex, cig2int):
    '''
    takes a cigar string and returns arrays for lengths and codes of cigar ops

    Parameters
    ----------
    cigar_string: str
        cigar string of alignment
    cigar_regex: re.pattern
        regex for parsing cigar ops
    cig2int: str.maketrans
        translates cigar opcodes into integers

    Returns
    -------
    lengths_arr: np.array
        array of lengths of individual cigar operations
    opsSeq: np.array
        array of cigar ops as integers

    '''
    # takes a cigar string and returns arrays for lengths and opcodes
    # first split cigar into tuples of length, opcode
    parts = cigar_regex.findall(cigar_string)
    # cast them into lists
    lengths, ops = zip(*parts)
    # convert to array
    lengths_arr = np.array(lengths, dtype=np.uint32)
    # for the opcodes, first put them into a string and translate to integers
    ops_seq = ''.join(ops)
    ops_tr = ops_seq.translate(cig2int)
    opsSeq = np.fromstring(ops_tr, dtype='int8', sep='')
    return lengths_arr, opsSeq


@njit
def _loop_cigar(cigar_lengths, cigar_ops, int_seq, int_qual, chunk_length, qpos):
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
            print('operation ' + str(operation) + ' not accounted for in parsing the cigar \n')

    return query, qual


def _collect_increments(covered_rois, start, target_chromosome_rois, query_arr, qual_arr, qt):
    # for each of the ROIs that was overlapped, grab the counts
    increments = []
    for roi in covered_rois:
        # calculate the offset of the ROI from the beginning of the read
        roi_offset = roi - start
        # get the base
        if qual_arr[roi_offset] >= qt:
            roi_coordinate = target_chromosome_rois[roi]
            increments.append((roi_coordinate, query_arr[roi_offset]))
    # transform to numpy array
    increments = np.array(increments, dtype=np.uint32)
    return increments


def _collect_increments_wg(covered_rois, target_chromosome_rois, query_arr):
    # faster version if we consider whole genomes where ROIs are always consecutive
    roi_coords = np.arange(target_chromosome_rois[covered_rois[0]],
                           target_chromosome_rois[covered_rois[-1]] + 1)
    increments = np.swapaxes(np.vstack((roi_coords, query_arr)), 0, 1)
    return increments








