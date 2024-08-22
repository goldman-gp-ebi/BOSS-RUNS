import numpy as np

from boss.readlengthdist import ReadlengthDist




def test_update(zymo_read_batch):
    r = ReadlengthDist()
    r.update(read_lengths=zymo_read_batch.read_lengths)
    assert int(r.lam) == 4700
    assert r.longest_read == 32099
    assert np.allclose(
        r.approx_ccl,
        np.array([1647, 2280, 2810, 3305, 3840, 4379, 5045, 5867, 7015, 9768])
    )




def test_update_noreads():
    r = ReadlengthDist()
    r.update(read_lengths={'a': 1, 'b': 2, 'c': 3})
    # remains at default
    assert int(r.lam) == 6000
    # not set if no valid reads observed
    assert hasattr(r, 'longest_read') is False
    # also stays at default
    assert np.allclose(
        r.approx_ccl,
        np.array([1167,  2729,  3903,  4918,  5866,  6808,  7797,  8912, 10321, 12713])
    )



def test_update_whale(zymo_read_batch):
    # add a whale read
    r = ReadlengthDist()
    rid = list(zymo_read_batch.read_lengths.keys())[3]
    zymo_read_batch.read_lengths[rid] = 2_222_222
    r.update(read_lengths=zymo_read_batch.read_lengths)
    assert int(r.lam) == 4802
    assert r.longest_read == 999_999
    # these are also slightly different due to the whale
    assert np.allclose(
        r.approx_ccl,
        np.array([1647, 2280, 2810, 3305, 3840, 4381, 5045, 5867, 7016, 9769])
    )







