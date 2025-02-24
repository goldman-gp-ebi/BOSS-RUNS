import subprocess
from types import SimpleNamespace
from importlib.metadata import version

import numpy as np
from numpy.typing import NDArray


def empty_file(path: str) -> None:
    """
    Create an empty file at the specified path.

    :param path: The path to the file.
    """
    with open(path, 'w'):
        pass
    return


def init_logger(logfile: str, args: SimpleNamespace) -> None:
    """
    Initialize the logger with the given logfile and log the arguments.

    :param logfile: The path to the logfile.
    :param args: The arguments to log.
    """
    empty_file(logfile)
    import logging
    logging.basicConfig(format='%(asctime)s %(message)s',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(f"{logfile}"), logging.StreamHandler()])

    logging.info("BOSS RUNS/AEONS")
    logging.info(f"{version('boss_runs')}")
    logging.info('\n')
    for a, aval in args.__dict__.items():
        logging.info(f'{a} {aval}')
    logging.info('\n')


def execute(command: str) -> tuple[str, str]:
    """
    Execute a command in a shell and return the stdout and stderr.

    :param command: The command to execute.
    :return: The stdout and stderr as a tuple.
    """
    # create the unix process
    running = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               encoding='utf-8', shell=True)
    # run on a shell and wait until it finishes
    stdout, stderr = running.communicate()
    return stdout, stderr



def write_logs(stdout: str, stderr: str, basename: str) -> None:
    """
    Write the stdout and stderr from a subprocess to files.

    :param stdout: The stdout output.
    :param stderr: The stderr output.
    :param basename: The base name for the log files.
    """
    # write stdout and stderr from subprocess to file
    with open(f'{basename}.out', 'a') as outf:
        outf.write(stdout)
        outf.write('\n')
    with open(f'{basename}.err', 'a') as errf:
        errf.write(stderr)
        errf.write('\n')


def binc(arguments: tuple[NDArray, NDArray]) -> NDArray:
    """
    Helper function for multithreaded bincounting

    :param arguments: Container of array and weights to count
    :return:
    """
    a, b = arguments
    res = np.bincount(a, weights=b)
    return res


def reverse_complement(dna: str) -> str:
    '''
    Return the reverse complement of a dna string. Used when parsing the cigar
    of a read that mapped on the reverse strand.

    :param dna: string of characters of the usual alphabet
    :return: reverse complemented input dna
    '''
    trans = str.maketrans('ATGC', 'TACG')
    rev_comp = dna.translate(trans)[::-1]
    return rev_comp


def random_id(k: int = 20) -> str:
    """
    Generate a random alphanumeric ID of length k.

    :param k: The length of the ID.
    :return: The generated random ID.
    """
    import random
    import string
    x = ''.join(random.choices(string.ascii_letters + string.digits, k=k))
    return x


def load_gfa(gfa_path: str) -> dict[str, str]:
    """
    Load a GFA file and return a dictionary of sequences.

    :param gfa_path: The path to the GFA file.
    :return: A dictionary mapping header strings to sequence strings.
    """
    def _load_gfa(infile: str):
        with open(infile, 'r') as gfa_file:
            for line in gfa_file:
                if line.startswith('S'):
                    ll = line.split('\t')
                    header = ll[1]
                    seq = ll[2]
                    yield header, seq

    sequences = {}
    for hd, sq in _load_gfa(gfa_path):
        sequences[hd] = sq
    return sequences


# def find_blocks_generic(arr: NDArray, x: int, min_len: int) -> NDArray:
#     """
#     Find blocks in the array that match x.
#
#     :param arr: The input array.
#     :param x: The value to find blocks of.
#     :param min_len: The minimum length of blocks to report.
#     :return: An array containing the start and end positions of the blocks.
#     """
#     # find run starts
#     x_pos = np.where(arr == x)[0]
#
#     if x_pos.shape[0] == 0:
#         return np.array([])
#
#     # diff between neighboring loc
#     x_diff = np.diff(x_pos)
#     # if diff > 1: new block
#     big_dist = np.where(x_diff > 1)[0]
#     # the first entry is a block start and then all other where a big change happens
#     # also each change is a block end, and the last entry of course as well
#     block_ranges = np.concatenate((np.array([x_pos[0]]), x_pos[big_dist + 1],
#                                    x_pos[big_dist] + 1, np.array([x_pos[-1] + 1])))
#     blocks = block_ranges.reshape(big_dist.shape[0] + 1, 2, order='F')
#     # only report blocks longer than min_len
#     blocks_filt = blocks[np.where(blocks[:, 1] - blocks[:, 0] > min_len)[0], :]
#     return blocks_filt


def find_blocks_ge(arr: NDArray, x: float, min_len: int) -> NDArray:
    """
    Find blocks in the array where the value is greater than or equal to x.

    :param arr: The input array.
    :param x: The value to compare against.
    :param min_len: The minimum length of blocks to report.
    :return: An array containing the start and end positions of the blocks.
    """
    # find run starts
    x_pos = np.where(arr >= x)[0]

    if x_pos.shape[0] == 0:
        return np.array([])

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



def window_sum(arr: NDArray, w: int) -> NDArray:
    """
    Calculate sums of non-overlapping windows

    :param arr: Input array
    :param w: Window size to calculate sums in
    :return: Reduced array of window sums
    """
    sumw = np.sum(arr[: (len(arr) // w) * w].reshape(-1, w), axis=1)
    return sumw



def adjust_length(original_size: int, expanded: NDArray) -> NDArray:
    """
    Adjust the size of an array after expanding during downsampling

    :param original_size: Target size of expanded array
    :param expanded: Input array after expansion
    :return: Adjusted input array
    """
    # after expanding arrays from binning
    lendiff = original_size - expanded.shape[0]
    if lendiff > 0:
        # original is longer than replacement
        repl = np.append(expanded, expanded[-lendiff:], axis=0)
    elif lendiff < 0:
        # original is shorter than replacement
        repl = expanded[: -abs(lendiff)]
    else:
        repl = expanded

    assert repl.shape[0] == original_size
    return repl




