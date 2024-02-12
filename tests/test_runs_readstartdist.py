import pytest
import numpy as np

import boss.runs.readstartdist as br_rsd
import boss.runs.reference as refc
import boss.paf



@pytest.fixture
def reference():
    r = refc.Reference(ref="../data/zymo.fa", mmi="../data/zymo.fa.mmi")
    return r

@pytest.fixture
def contigs():
    r = refc.Reference(ref="../data/zymo.fa", mmi="../data/zymo.fa.mmi")
    return r.contigs


@pytest.fixture
def read_start_dist(contigs):
    rsd = br_rsd.ReadStartDist(contigs=contigs)
    return rsd


@pytest.fixture
def paf_file():
    p = "../data/ERR3152366_10k.paf"
    return p


@pytest.fixture
def paf_dict(paf_file):
    return boss.paf.Paf.parse_PAF(paf_file, min_len=1)


def test_init(read_start_dist, reference):
    assert len(read_start_dist.read_starts) == 9
    # number of total windows
    assert read_start_dist.total_len == 15502
    # target size for expanding fhat
    # this should be: total genome size / 100
    downsampled_genome_size = reference.n_sites // 100
    assert read_start_dist.target_size == downsampled_genome_size


def test_count_read_starts(read_start_dist, paf_dict, reference):
    read_starts_merged = read_start_dist.merge()
    assert read_starts_merged.shape == (15502, 2)
    assert np.sum(read_starts_merged) == 0
    read_start_dist.count_read_starts(paf_dict=paf_dict)
    read_starts_merged = read_start_dist.merge()
    assert read_starts_merged.shape == (15502, 2)
    len_paf_dict = len(paf_dict)
    assert np.sum(read_starts_merged) == len_paf_dict - 1
    # update fhat after counting read starts
    fhat = read_start_dist.update_f_pointmass()
    # assert fhat.shape == (15502, 2)
    assert np.isclose(fhat.min(), 1.21414333e-06)
    assert np.isclose(fhat.max(), 1.578982876202093e-05)
    downsampled_genome_size = reference.n_sites // 100
    assert fhat.shape == (downsampled_genome_size, 2)

