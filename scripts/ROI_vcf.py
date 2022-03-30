import argparse
import requests
import time
import sys
from collections import defaultdict

# non-std lib
import numpy as np
import pandas as pd




def fetch_endpoint(req, content_type):
    '''
    semi-generic call for ensembl API
    retrieve either info about a sequence or the sequence itself
    Parameters
    ----------
    req: str
        the request
    content_type: str
        either get info or sequence

    Returns
    -------
    r: str
        results are either returned as json (info) or text (sequence)
    '''
    r = requests.get(f"http://rest.ensembl.org/{req}", headers={"Content-Type": content_type})

    # check return value
    if r.status_code in (500, 502, 503, 504):
        time.sleep(2)
        return fetch_endpoint(req, content_type)
    elif r.status_code == 429:
        # get retry-after
        if int(r.headers['X-RateLimit-Remaining']) < 50:
            print("waiting for rate limit reset...")
            time.sleep(float(r.headers['X-RateLimit-Reset']))
        # relaunch after waiting time
        print("sleeping just a little...")
        time.sleep(3)
        return fetch_endpoint(req, content_type)

    elif r.status_code == 200:
        # info is returned as json, sequence as text
        if content_type == 'application/json':
            return r.json()
        else:
            return r.text

    else:
        r.raise_for_status()
        sys.exit()



def load_chroms(chroms):
    cset = set()
    with open(chroms, 'r') as cn:
        for line in cn:
            ll = line.strip()
            cset.add(ll)
    return cset


def grab_sites(vcf, chroms=None, lim=0):
    # collect sites in a list
    sites = []
    i = 0  # limit of sites
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
            # if no chromosome names are given to subset
            if not chroms:
                chrom = ll[0].split(':')[-1]
                chrom = chrom.replace('chr', '')
            else:
                chrom_vcf = ll[0]
                if not chrom_vcf in chroms:
                    # skip this site
                    continue
            # position of site
            position = ll[1]
            # base at site
            ref = ll[3]
            sites.append((chrom, int(position), ref))
            # testing
            if lim:
                i += 1
                if i > lim:
                    break
    # put sites into a dictionary
    pos = defaultdict(list)
    for chrom, p, ref in sites:
        pos[chrom].append(p)
    pos = {k: np.array(v) for k, v in pos.items()}
    return pos



def chrom_length_dict(species, main_only=True, chroms=None):
    genome_request = f"info/assembly/{species}"
    info = fetch_endpoint(req=genome_request, content_type='application/json')

    chromdict = dict()
    # grab top level regions
    top_level_regions = info['top_level_region']
    for t in top_level_regions:
        chromdict[t['name']] = t['length']

    if main_only:
        chromdict = {k: v for k, v in chromdict.items() if len(k) < 5}
    # optionally subset to specific chroms
    if chroms:
        chromdict = {k: v for k, v in chromdict.items() if k in chroms}
    print("finished fetching chromosomes of species")
    return chromdict



def genome2roi_mapping(pos, n=0):
    # array that can be pickled and packed
    genome2roi = []
    # array of tuples with (chrom, genome_coord, roi_coord)
    for chrom, pos_arr in pos.items():
        for p in pos_arr:
            neighbors = np.arange(p - n, p + n + 1)
            for neighbor_pos in neighbors:
                genome2roi.append((chrom, int(neighbor_pos)))

    genome2roi = pd.DataFrame(genome2roi)
    # deduplicate in case we have the same position multiple times
    genome2roi_dedup = genome2roi.drop_duplicates()
    # add col with roi indices
    nrows = genome2roi_dedup.shape[0]
    roi_indices = np.arange(nrows)
    genome2roi_dedup.loc[:, 'roi_ind'] = roi_indices
    # back to array
    genome2roi_arr = np.array(genome2roi_dedup)
    # make coordinate translation dict
    genome2roi = defaultdict(dict)
    for chrom, genome_pos, roi_pos in genome2roi_arr:
        genome2roi[chrom][int(genome_pos)] = int(roi_pos)
    return genome2roi



def write_masks(genome2roi, chrom_lengths, out):
    # prep chromosome lengths
    chrom_length_arr = np.array([(c, clen) for c, clen in chrom_lengths.items()])
    # save ROI mask
    np.savez_compressed(f'{out}.mask',
                        chrom_lengths=chrom_length_arr,
                        genome2roi=genome2roi)

    for key, val in genome2roi.items():
        print(key)
        print(len(val))


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', type=str, required=True)
    parser.add_argument('--out', type=str, required=True)
    parser.add_argument('--species', type=str, default='hsap', required=True)
    parser.add_argument('--chroms', type=str, default=None)
    args = parser.parse_args()
    return args

#%%

# -------
# SimpleNamespace ...
# species = "hsap"
# vcf = ""
# out = 'test'
# from types import SimpleNamespace
# args = SimpleNamespace(species=species, vcf=vcf, out=out)
# -------

args = get_args()

if args.chroms:
    chroms = load_chroms(args.chroms)
else:
    chroms = None
# grab sites from vcf file
sites = grab_sites(vcf=args.vcf, chroms=chroms)
# get chromosome lengths
chrom_lengths = chrom_length_dict(species=args.species, chroms=chroms)
# create site mapping
genome2roi = genome2roi_mapping(sites)
# write output file
write_masks(genome2roi=genome2roi, chrom_lengths=chrom_lengths, out=args.out)
print("DONE")







