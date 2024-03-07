import pytest
from pathlib import Path

import numpy

import boss.sampler


sources = ["../data/ERR3152366_10k.fq", "../data/ERR3152366_10k.fq.gz"]
paf_files = ["../data/ERR3152366_10k.paf", "../data/ERR3152366_10k_trunc.paf"]


@pytest.fixture
def fqs(request):
    # Initialise a FastqStream but remove offset files first to make sure they are recreated
    Path(f'{request.param}.offsets.npy').unlink(missing_ok=True)
    fqs = boss.sampler.FastqStream_mmap(source=request.param, batchsize=5, maxbatch=10)
    return fqs


@pytest.mark.parametrize('fqs', sources, indirect=True)
def test_init(fqs):
    # initialise the fastqstream with two different fastq files
    # one of which is gzipped
    assert fqs.offsets.shape == (11, 5)
    assert type(fqs.offsets[5, 2]) is numpy.uint64
    arb_val = fqs.offsets[5, 2]
    assert arb_val == 268869
    # reload with shuffling
    fqs._load_offsets(shuffle=True, batchsize=5, maxbatch=10)
    assert fqs.offsets.shape == (11, 5)
    assert type(fqs.offsets[5, 2]) is numpy.uint64


@pytest.mark.parametrize('fqs', sources, indirect=True)
def test_read_batch(fqs):
    # load a batch of reads from a stream, again with fastq and gzipped fastq
    fqs.read_batch()
    assert len(fqs.read_ids) == 5
    assert len(fqs.read_lengths) == 5
    assert len(fqs.read_sequences) == 5
    assert fqs.total_bases == 23070


@pytest.fixture
def paf_stream():
    # initialise a pafstream and remove offsets first, as above
    Path(f'{paf_files[0]}.offsets').unlink(missing_ok=True)
    Path(f'{paf_files[1]}.offsets').unlink(missing_ok=True)
    p = boss.sampler.PafStream(paf_full=paf_files[0], paf_trunc=paf_files[1])
    return p


def test_paf_init(paf_stream):
    assert len(paf_stream.offsets_full) == 9120
    assert paf_stream.offsets_full['ERR3152366.998'][0] == 1213253



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
def sampler_aeons():
    # initialise the wrapper class for the fastq and paf streaming
    s = boss.sampler.Sampler(
            source=sources[0],
            maxbatch=10,
            batchsize=5,
    )
    return s


@pytest.mark.xfail(raises=ValueError)
def test_request_too_much():
    _ = boss.sampler.Sampler(
        source=sources[0],
        maxbatch=300,
        batchsize=50,
    )


def test_sample(sampler):
    # this is the main function to retrieve new data
    # samples reads first and then grabs their mappings
    r_seqs, r_quals, paf_f, paf_t = sampler.sample()
    assert len(sampler.fq_stream.read_ids) == 5
    assert len(sampler.fq_stream.read_lengths) == 5
    assert len(r_seqs) == 5
    assert sampler.fq_stream.total_bases == 23070
    assert paf_f.startswith("ERR3152366")
    assert paf_t.startswith("ERR3152366")



def test_sample_nopaf(sampler_aeons):
    # this is the main function to retrieve new data
    # samples reads first and then grabs their mappings
    r_seqs, r_quals, paf_f, paf_t = sampler_aeons.sample()
    assert len(sampler_aeons.fq_stream.read_ids) == 5
    assert len(sampler_aeons.fq_stream.read_lengths) == 5
    assert len(r_seqs) == 5
    assert sampler_aeons.fq_stream.total_bases == 23070
    assert not paf_f
    assert not paf_t


@pytest.mark.xfail(raises=ValueError)
def test_exhaust_sampler(sampler):
    while True:
        _ = sampler.sample()








