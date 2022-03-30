from itertools import groupby
from collections import defaultdict
import argparse

# non-std lib
import numpy as np
import pandas as pd
from natsort import natsorted



def read_fa(fh):
    # iterator for all headers in the file
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip().split(' ')[0]
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield headerStr, seq


def chromosome_dict(fa):
    chrom_dict = {}
    with open(fa, 'r') as fasta:
        for cname, cseq in read_fa(fasta):
            chrom_dict[cname] = cseq
    # apply natsorting
    cnames_sorted = natsorted(list(chrom_dict.keys()))
    chrom_dict_sort = {cname: chrom_dict[cname] for cname in cnames_sorted}
    chrom_lengths_sort = {cname: len(chrom_dict[cname]) for cname in cnames_sorted}
    return chrom_dict_sort, chrom_lengths_sort


def write_bases(seq_dict, out):
    base_raw = open(f'{out}.ref.fa', 'w')
    base_raw.write('>chr\n')

    for cname, cseq in seq_dict.items():
        # replace ambiguous bases
        iupac = "NRYKMSWBDHV"
        for i in iupac:
            cseq = cseq.replace(i, "A")
        base_raw.write(cseq)
    base_raw.close()


def position_dict(length_dict):
    posdict = dict()
    for name, size in length_dict.items():
        posdict[name] = np.arange(size)
    return posdict


def genome2roi_mapping(pos):
    # array mask
    genome2roi = []
    # array of tuples with (chrom, genome_coord, roi_coord)
    for chrom, pos_arr in pos.items():
        for p in pos_arr:
            genome2roi.append((chrom, int(p)))

    genome2roi = pd.DataFrame(genome2roi)
    # add col with roi indices
    nrows = genome2roi.shape[0]
    roi_indices = np.arange(nrows)
    genome2roi.loc[:, 'roi_ind'] = roi_indices
    # back to array
    genome2roi_arr = np.array(genome2roi)
    # final shape: coordinate dictionary
    genome2roi = defaultdict(dict)
    for chrom, genome_pos, roi_pos in genome2roi_arr:
        genome2roi[chrom][int(genome_pos)] = int(roi_pos)
    return genome2roi


def write_masks(genome2roi, chrom_lengths, out):
    # prep the chromosome lengths
    chrom_length_arr = np.array([(c, clen) for c, clen in chrom_lengths.items()])

    # save the masks to a file
    np.savez_compressed(f'{out}.mask',
                        chrom_lengths=chrom_length_arr,
                        genome2roi=genome2roi)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fa', type=str, required=True)
    parser.add_argument('--out', type=str, required=True)
    args = parser.parse_args()
    return args

#%%

# -------
# SimpleNamespace ...
# fa = ""
# out = 'test'
# from types import SimpleNamespace
# args = SimpleNamespace(fa=fa, out=out)
# -------


args = get_args()
# get chromosome names and lengths
seq_dict, length_dict = chromosome_dict(fa=args.fa)
# write bases to single ref
write_bases(seq_dict=seq_dict, out=args.out)
# position dict contains ranges for whole genomes
pos = position_dict(length_dict=length_dict)
# mapping of gpos to rpos
genome2roi = genome2roi_mapping(pos=pos)
# write outfile
write_masks(genome2roi=genome2roi, chrom_lengths=length_dict, out=args.out)
print("DONE")



