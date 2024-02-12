import pytest

import numpy as np

import boss.runs.sequences as brs
from boss.batch import FastqBatch
from boss.paf import Paf


# with deletions (default) there are 5 states
# which means 15 genotypes for diploids
@pytest.mark.parametrize('ploidy, b, g', [
    (1, 4, 5),
    (2, 4, 15)
])
def test_priors(ploidy, b, g):
    p = brs.Priors(ploidy=ploidy)
    assert p.len_b == b + 1
    assert p.len_g == g
    assert p.phi_stored.shape == (b + 1, g, 1000)
    assert p.priors.shape == (b, g)


@pytest.mark.xfail(raises=ValueError)
def test_priors_ploidy():
    brs.Priors(ploidy=3)


def test_uniform():
    p = brs.Priors(ploidy=1)
    p.uniform_priors()
    assert np.all(np.isclose(p.priors, p.priors[0]))


# without deletions there are 4 states
# and 10 genotypes for diploids
@pytest.mark.parametrize('diploid, del_err, b, g', [
    (False, 0, 4, 4),
    (True, 0, 4, 10)
])
def test_generate_phi(diploid, del_err, b, g):
    lenb, leng, phi = brs.Priors._generate_phi(diploid=diploid, deletion_error=del_err)
    assert lenb == b
    assert leng == g
    assert phi.shape == (b, g)



@pytest.mark.parametrize('del_err, b, g', [
    (0, 4, 4),
    (0.5, 4, 5)
])
def test_haploid_priors(del_err, b, g):
    priors = brs.Priors._haploid_priors(deletion_error=del_err)
    assert priors.shape == (b, g)


@pytest.mark.parametrize('del_err, b, g', [
    (0, 4, 10),
    (0.5, 4, 15)
])
def test_diploid_priors(del_err, b, g):
    priors = brs.Priors._diploid_priors(deletion_error=del_err)
    assert priors.shape == (b, g)





# CoverageConverter
@pytest.fixture
def cov_conv():
    return brs.CoverageConverter()


@pytest.fixture
def zymo_reads():
    fq = "../data/ERR3152366_10k.fq"
    batch = FastqBatch(fq_files=[fq])
    return batch


@pytest.fixture
def paf_dict():
    return Paf.parse_PAF(paf_file="../data/ERR3152366_10k.paf", min_len=1)


def test_convert_records(cov_conv, zymo_reads, paf_dict):
    incr = cov_conv.convert_records(
        paf_dict=paf_dict,
        seqs=zymo_reads.read_sequences,
        quals=zymo_reads.read_qualities
    )
    assert incr
    assert len(incr) == 8


@pytest.fixture
def scoring():
    return brs.Scoring()


def test_init_scoring(scoring):
    assert np.isclose(scoring.score0, 0.04969294)
    assert np.isclose(scoring.ent0, 0.09302521)


def test_score_array(scoring):
    scoring.init_score_array()
    assert scoring.score_arr.shape == (40, 40, 40, 40, 40, 4)
    assert scoring.entropy_arr.shape == (40, 40, 40, 40, 40, 4)
    assert np.isclose(scoring.score_arr[28, 0, 0, 0, 0, 3], 3.834200141940696e-44)
    assert np.isclose(scoring.entropy_arr[28, 0, 0, 0, 0, 3], 3.834200141940696e-44)
    assert np.isclose(scoring.score_arr[2, 0, 0, 0, 0, 3], 0.17253973305650225)
    assert np.isclose(scoring.entropy_arr[2, 0, 0, 0, 0, 3], 0.22957118271635163)




