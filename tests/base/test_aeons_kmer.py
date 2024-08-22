import pytest
import numpy as np

import boss.aeons.kmer as kmer
from boss.aeons.sequences import Sequence


@pytest.fixture
def kmer_counter():
    return kmer.KmerCounter()


@pytest.fixture
def tetramer_dist():
    return kmer.TetramerDist()


def test_tetra_zscores(sampler, kmer_counter):
    r_seqs, _, _, _ = sampler.sample()
    seq = list(r_seqs.values())[0]
    tm_zscores = kmer_counter.tetra_zscores(seq=seq)
    assert np.isclose(tm_zscores['ACGT'], 1.4174142762462025)
    assert len(tm_zscores) == 256


def test_euclidean_dist(tetramer_dist, sampler):
    r_seqs, _, _, _ = sampler.sample()
    seq1_rid, seq1_s = list(r_seqs.items())[0]
    seq2_rid, seq2_s = list(r_seqs.items())[1]
    seqo1 = Sequence(seq1_rid, seq1_s)
    seqo2 = Sequence(seq2_rid, seq2_s)
    euc = tetramer_dist.euclidean_dist(seqo1=seqo1, seqo2=seqo1)
    assert euc == 0.0
    euc = tetramer_dist.euclidean_dist(seqo1=seqo1, seqo2=seqo2)
    assert np.isclose(euc, 0.015520185523573026)
    euc = tetramer_dist.euclidean_dist(seqo1=seqo2, seqo2=seqo1)
    assert np.isclose(euc, 0.015520185523573026)





