import logging
import pytest
from pathlib import Path
from io import StringIO

import numpy as np

import boss.batch
from boss.sampler import Sampler
from boss.paf import Paf


np.random.seed(1)


sources = ["../data/ERR3152366_10k.fq", "../data/ERR3152366_10k.fq.gz"]
paf_files = ["../data/ERR3152366_10k.paf", "../data/ERR3152366_10k_trunc.paf"]


@pytest.fixture
def sampler():
    # initialise the wrapper class for the fastq and paf streaming
    s = boss.sampler.Sampler(
            source=sources[0],
            paf_full=paf_files[0],
            paf_trunc=paf_files[1],
            maxbatch=10,
            batchsize=5,
    )
    return s


@pytest.fixture
def fq_files():
    f = Path("../data/fastq_pass/").glob("F*")
    f = list(f)
    f.sort()
    # make one filepath str to test input types
    f[1] = str(f[1])
    # add file without ch=
    f.append("../data/fastq_pass_no_ch/ERR3152366_1k_0.fq")
    return f



@pytest.mark.parametrize("channels, rlen, total", [
    ({1, 2, 3, 4}, 741, 69275),
    ({464, 456, 32, 24, 16, 8, 351, 343}, 518, 135997),
    ({464, 456, 32, 24, 16, 8, 351, 343, 0, 555}, 518, 135997),
    ({}, 319, 8261497),
])
def test_batch(fq_files, channels, rlen, total):
    logging.info(fq_files)
    batch = boss.batch.FastqBatch(fq_files, channels=channels)
    # grab a random element for testing
    rid, seq = list(batch.read_sequences.items())[0]
    assert rid in batch.read_ids
    assert batch.read_lengths[rid] == rlen
    assert batch.total_bases == total



# let fail intentionally
@pytest.mark.xfail(raises=TypeError)
def test_batch_type():
    batch = boss.batch.FastqBatch(fq_files=[1])
    assert batch



@pytest.fixture
def read_cache():
    rc = boss.batch.ReadCache(batchsize=10, dumptime=10_000)
    return rc


def test_init(read_cache):
    assert read_cache
    assert Path("00_reads/control_0.fa").is_file()


def test_update_times_aeons(read_cache, fq_files):
    np.random.seed(1)
    batch = boss.batch.FastqBatch(fq_files)
    # create some arbitrary decisions
    reads_decision = {}
    for rid, seq in batch.read_sequences.items():
        if np.random.choice((0, 1)):
            reads_decision[rid] = seq
        else:
            reads_decision[rid] = seq[: 400]
    read_cache.update_times_aeons(read_sequences=batch.read_sequences, reads_decision=reads_decision)
    assert read_cache.time_control == 8_264_497
    assert read_cache.time_boss == 5_819_980


def test_update_times_runs(read_cache):
    s = boss.sampler.Sampler(
        source=sources[0],
        paf_full=paf_files[0],
        paf_trunc=paf_files[1],
        maxbatch=10,
        batchsize=50,
    )
    r_seqs, r_quals, paf_f, paf_t = s.sample()
    total_bases = np.sum([len(s) for s in r_seqs.values()])
    total_n = len(r_seqs)
    # sample some unmapped and rejected reads
    pafd_f = Paf.parse_PAF(StringIO(paf_f))
    pafd_t = Paf.parse_PAF(StringIO(paf_t))
    s0 = np.random.choice(np.array([rid for rid in r_seqs.keys()]), size=total_n // 2)
    unmapped = np.random.choice(s0, size=s0.shape[0] // 2)
    rejected = set(s0) - set(unmapped)
    n_unmapped = len(unmapped)
    n_rejected = len(rejected)
    logging.info(f"unmapped {n_unmapped}, reject {n_rejected}")
    paf_dict = {}
    for rid in r_seqs.keys():
        if rid in pafd_t.keys() and rid in pafd_f.keys():
            if rid in rejected:
                paf_dict[rid] = pafd_t[rid]
            else:
                paf_dict[rid] = pafd_f[rid]

    read_cache.update_times_runs(
        total_bases=total_bases,
        paf_dict=paf_dict,
        n_unmapped=n_unmapped,
        n_reject=n_rejected)

    assert read_cache.time_control == 249464
    assert read_cache.time_boss == 190181



def test_fill_cache(read_cache):
    s = boss.sampler.Sampler(
        source=sources[0],
        maxbatch=10,
        batchsize=50,
    )
    r_seqs, r_quals, paf_f, paf_t = s.sample()
    # create some arbitrary decisions
    reads_decision = {}
    for rid, seq in r_seqs.items():
        if np.random.choice((0, 1)):
            reads_decision[rid] = seq
        else:
            reads_decision[rid] = seq[: 400]

    read_cache.fill_cache(read_sequences=r_seqs,
                          reads_decision=reads_decision)

    assert len(read_cache.cache_boss) == 50
    assert len(read_cache.cache_control) == 50

    # update time, then fill cache again, which triggers dump
    read_cache.update_times_aeons(read_sequences=r_seqs,
                                  reads_decision=reads_decision)

    assert read_cache.time_control == 249_464
    assert read_cache.time_boss == 151_124

    read_cache.fill_cache(read_sequences=r_seqs,
                          reads_decision=reads_decision)

    assert len(read_cache.cache_boss) == 0
    assert len(read_cache.cache_control) == 0