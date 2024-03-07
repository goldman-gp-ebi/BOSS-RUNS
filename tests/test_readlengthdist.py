import pytest

import numpy as np

import boss.readlengthdist
from boss.batch import FastqBatch



@pytest.fixture
def readlengthdist():
    return boss.readlengthdist.ReadlengthDist()



@pytest.fixture
def zymo_reads():
    fq = "../data/ERR3152366_10k.fq"
    batch = FastqBatch(fq_files=[fq])
    return batch



def test_update(readlengthdist, zymo_reads):
    readlengthdist.update(read_lengths=zymo_reads.read_lengths)
    assert int(readlengthdist.lam) == 4700
    assert readlengthdist.longest_read == 32099
    assert np.allclose(
        readlengthdist.approx_ccl,
        np.array([1647, 2280, 2810, 3305, 3840, 4379, 5045, 5867, 7015, 9768])
    )



def test_update_whale(readlengthdist, zymo_reads):
    # add a whale read
    rid = list(zymo_reads.read_lengths.keys())[3]
    zymo_reads.read_lengths[rid] = 2_222_222
    readlengthdist.update(read_lengths=zymo_reads.read_lengths)
    assert int(readlengthdist.lam) == 4802
    assert readlengthdist.longest_read == 999_999
    # these are also slightly different due to the whale
    assert np.allclose(
        readlengthdist.approx_ccl,
        np.array([1647, 2280, 2810, 3305, 3840, 4381, 5045, 5867, 7016, 9769])
    )







