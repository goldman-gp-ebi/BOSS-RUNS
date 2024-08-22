import boss.runs.abundance_tracker as abt


def test_abt(zymo_ref):
    a = abt.AbundanceTracker(contigs=zymo_ref.contigs)
    assert len(a.read_counts) == 9
    assert sum(a.read_counts.values()) == 0


def test_update(zymo_ref, paf_dict):
    a = abt.AbundanceTracker(contigs=zymo_ref.contigs)
    a.update(n=100, paf_dict=paf_dict)
    assert sum(a.read_counts.values()) == 9120


def test_update_empty(zymo_ref):
    a = abt.AbundanceTracker(contigs=zymo_ref.contigs)
    a.update(n=100, paf_dict={})
    assert sum(a.read_counts.values()) == 0


