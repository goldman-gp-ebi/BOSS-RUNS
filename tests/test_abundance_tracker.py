import pytest

import boss.runs.abundance_tracker as abt
from boss.runs.reference import Reference
from boss.paf import Paf



@pytest.fixture
def contigs():
    r = Reference(ref="../data/zymo.fa", mmi="../data/zymo.fa.mmi")
    return r.contigs


@pytest.fixture
def paf_file():
    p = "../data/ERR3152366_10k.paf"
    return p

@pytest.fixture
def paf_dict(paf_file):
    return Paf.parse_PAF(paf_file, min_len=1)


def test_abt(contigs):
    a = abt.AbundanceTracker(contigs=contigs)
    assert len(a.read_counts) == 9
    assert sum(a.read_counts.values()) == 0


def test_update(contigs, paf_dict):
    a = abt.AbundanceTracker(contigs=contigs)
    a.update(n=100, paf_dict=paf_dict)
    assert sum(a.read_counts.values()) == 9120



def test_update_empty(contigs):
    a = abt.AbundanceTracker(contigs=contigs)
    a.update(n=100, paf_dict={})
    assert sum(a.read_counts.values()) == 0




