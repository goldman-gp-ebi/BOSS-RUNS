import logging
from shutil import which
import sys
import subprocess
import argparse
from itertools import groupby
from collections import defaultdict
from string import ascii_letters

# non-std lib
import numpy as np
import toml


def init_logger(logfile):
    logging.basicConfig(format='%(asctime)s %(message)s',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(f"{logfile}"), logging.StreamHandler()])


def find_exe(exe):
    exe_path = which(exe, path='/'.join(sys.executable.split('/')[0:-1]))
    return exe_path.strip()


def execute(command):
    running = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               encoding='utf-8', shell=True)
    # run on a shell and wait until it finishes
    stdout, stderr = running.communicate()
    return stdout, stderr


def spawn(comm):
    running = subprocess.Popen(comm, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               encoding='utf-8', shell=True)
    return running


def empty_file(path):
    # create file and make sure its empty
    with open(path, 'w'):
        pass
    return


def binc(arguments):
    # helper function for multi-threaded bincounting
    a, b = arguments
    res = np.bincount(a, weights=b)
    return res


def append_row(frame_dict, new_row):
    # this is used by save metrics
    # takes 2 dicts as input: one is a metric frame,
    # the other is a row to append to that
    for colname, rowval in new_row.items():
        frame_dict[colname].append(rowval)
    return frame_dict


def range_intersection(r1, r2):
    return len(range(max(r1.start, r2.start), min(r1.stop, r2.stop)))


def window_sum(arr, w):
    # sums of non-overlapping windows
    sumw = np.sum(arr[: (len(arr) // w) * w].reshape(-1, w), axis=1)
    return sumw


def adjust_length(original, expanded):
    # after expanding arrays from binning
    lendiff = original.shape[0] - expanded.shape[0]
    if lendiff > 0:
        # original is longer than replacement
        repl = np.append(expanded, expanded[-lendiff:], axis=0)
    elif lendiff < 0:
        # original is shorter than replacement
        repl = expanded[: -abs(lendiff)]
    else:
        repl = expanded

    assert repl.shape[0] == original.shape[0]
    return repl


def reverse_complement(dna):
    '''
    Return the reverse complement of a dna string. Used when parsing the cigar
    of a read that mapped on the reverse strand.

    Parameters
    ----------
    dna: str
        string of characters of the usual alphabet

    Returns
    -------
    rev_comp: str
        the reverse complemented input dna

    '''
    trans = str.maketrans('ATGC', 'TACG')
    rev_comp = dna.translate(trans)[::-1]
    return rev_comp


def read_fa(fh):
    # iterator for all headers in the file
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip().split(' ')[0]
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield headerStr, seq


def readfq(fp):
    """
    GENERATOR FUNCTION
    Read a fastq file and return the sequence
    Parameters
    ----------
    fp: _io.IO
        File handle for the fastq file.

    Yields
    -------
    desc: str
        The fastq read header
    name: str
        The read ID
    seq: str
        The sequence

    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for ll in fp:  # search for the start of the next record
                if ll[0] in ">@":  # fasta/q header line
                    last = ll[:-1]  # save this line
                    break
        if not last:
            break
        desc, name, seqs, last = last[1:], last[1:].partition(" ")[0], [], None
        for ll in fp:  # read the sequence
            if ll[0] in "@+>":
                last = ll[:-1]
                break
            seqs.append(ll[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield desc, name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for ll in fp:  # read the quality
                seqs.append(ll[:-1])
                leng += len(ll) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield desc, name, seq, "".join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield desc, name, seq, None  # yield a fasta record instead
                break


def integer_genome(chromosome_sequences):
    '''
    transform nucleotide genome to integer values.
    Unknown bases are set to integer 99.
    Parameters
    ----------
    chromosome_sequences: dict
        dictionary with strings of nucleotides as values

    Returns
    -------
    genome: np.array
        1d array of integers corresponding to nucleotides. Concatenated chromosomes

    '''
    # define the dict for translating the bases
    transDict = defaultdict(str)
    transDict.update({'A': 0, 'C': 1, 'G': 2, 'T': 3})
    for letter in ascii_letters:
        if letter not in transDict.keys():
            transDict[letter] = 0
    # create translation table from dictionary
    base2int = str.maketrans(dict(transDict))

    genome = []
    for chr_name, seq in chromosome_sequences.items():
        # translate and convert from byte string
        read_integer = seq.translate(base2int)
        intSeq = np.fromstring(read_integer, dtype='int8', sep='')
        genome.append(intSeq)
    # transform to array at the end
    genome = np.concatenate(genome)
    return genome


def find_blocks_generic(arr, x, min_len):
    # find run starts
    x_pos = np.where(arr == x)[0]

    if x_pos.shape[0] == 0:
        return []

    # diff between neighboring loc
    x_diff = np.diff(x_pos)
    # if diff > 1: new block
    big_dist = np.where(x_diff > 1)[0]
    # the first entry is a block start and then all other where a big change happens
    # also each change is a block end, and the last entry of course as well
    block_ranges = np.concatenate((np.array([x_pos[0]]), x_pos[big_dist + 1],
                                   x_pos[big_dist] + 1, np.array([x_pos[-1] + 1])))
    blocks = block_ranges.reshape(big_dist.shape[0] + 1, 2, order='F')
    # only report blocks longer than min_len
    blocks_filt = blocks[np.where(blocks[:, 1] - blocks[:, 0] > min_len)[0], :]
    return blocks_filt


def grab_br_channels(channels_toml, run_name):
    '''
    Get the channel numbers assigned to a condition on the flowcell

    Parameters
    ----------
    channels_toml: str
        path to a channels.toml file from readfish
    run_name: str
        header of the condition in the toml file

    Returns
    -------
    br_channels: set
        channel numbers of the BR condition
    '''
    # parse the channels TOML file
    toml_dict = toml.load(channels_toml)
    # find the condition that corresponds to BR
    br_key = 0
    for key in toml_dict["conditions"].keys():
        name = toml_dict["conditions"][key]["name"]
        if name == run_name:
            br_key = key
            break
    try:
        br_channels = set(toml_dict["conditions"][br_key]["channels"])
        logging.info("grabbing channel numbers for bossruns condition")
        return br_channels
    except UnboundLocalError:
        logging.info("run_name in parameter file not found in channel-specification toml. Exiting")
        exit()


class MyArgumentParser(argparse.ArgumentParser):
    '''
    subclassing argument parser from argparse
    overwriting how the lines of argument files are read
    by default ArgumentParser treats each line as a contiguous word
    i.e. --test 0 1 becomes '--test 0 1'
    this way, they are split into name, val
    '''
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()



def read_args_fromfile(parser, file):
    with open(file, 'r') as f:
        arguments = [ll for line in f for ll in line.rstrip().split()]
    args = parser.parse_args(args=arguments)
    return args




def setup_parser_sim():
    # SIMULATION VERSION
    parser = MyArgumentParser(fromfile_prefix_chars='@')
    #
    parser.add_argument('--ref', type=str, required=True)
    parser.add_argument('--ref_idx', type=str, default=None)
    parser.add_argument('--vcf', type=str, default=None)
    parser.add_argument('--run_name', default="br", type=str)
    #
    parser.add_argument('--ckp', default=None, type=str, help="path to checkpoint file for relaunching")
    parser.add_argument('--channels', dest='channels', type=str, default=None, help='path to channels.toml file')
    parser.add_argument('--fastq_dir', dest='fastq_dir', type=str, default=None, help='dir with playback data')
    #
    parser.add_argument('--fq', type=str, default=None, help='path to fastq for streaming')
    parser.add_argument('--paf_full', dest='paf_full', type=str, default=None, help='path to paf full for streaming')
    parser.add_argument('--paf_trunc', dest='paf_trunc', type=str, default=None, help='path to paf trunc for streaming')
    parser.add_argument('--batch_size', dest='batch_size', default=100, type=int, help='Number of reads in a batch')
    #
    parser.add_argument('--ploidy', default=1, help='1 == haploid, 2==diploid')
    return parser



def setup_parser_live():
    # LIVE VERSION
    parser = MyArgumentParser(fromfile_prefix_chars='@')
    #
    parser.add_argument('--ref', type=str, required=True, help='Path to reference')
    parser.add_argument('--ref_idx', type=str, default=None, help='Optional minimap index of reference')
    parser.add_argument('--vcf', type=str, default=None, help='Path to vcf file for ROIs')
    parser.add_argument('--run_name', type=str, default="br",
                        help='Experiment identifier. If multiple conditions: must match name in channels.toml')
    parser.add_argument('--ploidy', default=1, help='1 == haploid, 2 == diploid')
    parser.add_argument('--conditions', action='store_true',
                        help="Multiple conditions on a single flowcell, used to assign channels")
    #
    parser.add_argument('--device', required=True, type=str, help="Name of device/sequencing position in MinKNOW")
    parser.add_argument('--host', default='localhost', type=str, help="hostname of sequencing device")
    parser.add_argument('--port', default=None, type=str, help="port of sequencing device")
    #
    parser.add_argument('--wait', default=90, type=int, help='Period between strategy updates (sec.)')
    parser.add_argument('--ckp', default=None, type=str, help="path to checkpoint file for relaunching. Experimental")
    return parser


