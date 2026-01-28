import pytest
from pathlib import Path

import numpy

import boss.sampler

from ..constants import PATHS



@pytest.fixture
def fq_mmap():
    # Initialise a FastqStream but remove offset files first to make sure they are recreated
    Path(f'{PATHS.fastq}.offsets.npy').unlink(missing_ok=True)
    fq_mmap = boss.sampler.FastqStream_mmap(source=PATHS.fastq, batchsize=5, maxbatch=10)
    return fq_mmap


@pytest.fixture
def fq_mmap_gz():
    # Initialise a FastqStream but remove offset files first to make sure they are recreated
    Path(f'{PATHS.fastq_gz}.offsets.npy').unlink(missing_ok=True)
    fq_mmap_gz = boss.sampler.FastqStream_mmap(source=PATHS.fastq_gz, batchsize=5, maxbatch=10)
    return fq_mmap_gz


@pytest.mark.parametrize('fq', ['fq_mmap', 'fq_mmap_gz'])
def test_init(fq, request):
    fq = request.getfixturevalue(fq)
    # initialise the fastqstream with two different fastq files
    # one of which is gzipped
    assert fq.offsets.shape == (11, 5)
    assert type(fq.offsets[5, 2]) is numpy.uint64
    arb_val = fq.offsets[5, 2]
    assert arb_val == 268869
    # reload with shuffling
    fq._load_offsets(shuffle=True, batchsize=5, maxbatch=10)
    assert fq.offsets.shape == (11, 5)
    assert type(fq.offsets[5, 2]) is numpy.uint64


@pytest.mark.parametrize('fq', ['fq_mmap', 'fq_mmap_gz'])
def test_read_batch(fq, request):
    fq = request.getfixturevalue(fq)
    # load a batch of reads from a stream, again with fastq and gzipped fastq
    fq.read_batch()
    assert len(fq.read_ids) == 5
    assert len(fq.read_lengths) == 5
    assert len(fq.read_sequences) == 5
    assert fq.total_bases == 23070


def test_paf_init():
    # initialise a pafstream and remove offsets first, as above
    Path(f'{PATHS.paf}.offsets').unlink(missing_ok=True)
    Path(f'{PATHS.paf_trunc}.offsets').unlink(missing_ok=True)
    p = boss.sampler.PafStream(paf_full=PATHS.paf, paf_trunc=PATHS.paf_trunc)
    assert len(p.offsets_full) == 9120
    assert p.offsets_full['ERR3152366.998'][0] == 1213253  # type: ignore



@pytest.mark.xfail(raises=ValueError)
def test_request_too_much():
    _ = boss.sampler.Sampler(
        source=PATHS.fastq,
        maxbatch=300,
        batchsize=50,
    )


def test_sample(sampler):
    # this is the main function to retrieve new data
    # samples reads first and then grabs their mappings
    r_seqs, r_quals, r_barcodes, paf_f, paf_t = sampler.sample()
    assert len(sampler.fq_stream.read_ids) == 50
    assert len(sampler.fq_stream.read_lengths) == 50
    assert len(r_seqs) == 50
    assert sampler.fq_stream.total_bases == 246464
    assert paf_f.startswith("ERR3152366")
    assert paf_t.startswith("ERR3152366")



def test_sample_nopaf(sampler_nopaf):
    # this is the main function to retrieve new data
    # samples reads first and then grabs their mappings
    r_seqs, r_quals, r_barcodes, paf_f, paf_t = sampler_nopaf.sample()
    assert len(sampler_nopaf.fq_stream.read_ids) == 50
    assert len(sampler_nopaf.fq_stream.read_lengths) == 50
    assert len(r_seqs) == 50
    assert sampler_nopaf.fq_stream.total_bases == 246464
    assert not paf_f
    assert not paf_t


@pytest.mark.xfail(raises=ValueError)
def test_exhaust_sampler(sampler):
    while True:
        _ = sampler.sample()








