import pytest
from pathlib import Path

import boss.mapper
from boss.batch import FastqBatch
from boss.paf import PafLine


@pytest.fixture
def dummy_fasta():
    fa = '../data/test.fasta'
    with open(fa, 'w') as f:
        f.write('>chr1\nAAAAAAAAGCTGACTATCGATCGTAGCTAGCTGACTGATCGTG\n'
                '>chr2\nGTGATCGTGATCGATCGTAGTCGTGTGCTAGTGTGTGTGTGATCGT\n')
    return fa


@pytest.fixture
def zymo_fasta():
    fa = "../data/zymo.fa"
    return fa


@pytest.fixture
def mapper(zymo_fasta):
    mapper = boss.mapper.Mapper(ref=zymo_fasta)
    # test init of non-default parameters
    _ = boss.mapper.Mapper(ref=zymo_fasta, default=False)
    return mapper


@pytest.fixture
def zymo_reads():
    fq = "../data/ERR3152366_10k.fq"
    batch = FastqBatch(fq_files=[fq])
    return batch.read_sequences


@pytest.mark.xfail(raises=FileNotFoundError)
def test_non_existing_ref():
    _ = boss.mapper.Mapper(ref="./dummy.fa")


def test_indexer(dummy_fasta):
    _ = boss.mapper.Indexer(fasta=dummy_fasta, mmi=f'{dummy_fasta}.mmi')
    assert Path(f'{dummy_fasta}.mmi').is_file()


def test_indexer_zymo(zymo_fasta):
    _ = boss.mapper.Indexer(fasta=zymo_fasta, mmi=f'{zymo_fasta}.mmi')
    assert Path(f'{zymo_fasta}.mmi').is_file()


def test_map_sequences(mapper, zymo_reads):
    paf_dict = mapper.map_sequences(sequences=zymo_reads)
    assert len(paf_dict) == 9108
    s = paf_dict['ERR3152366.998'][0]
    assert type(s) is PafLine
    assert s.qname == "ERR3152366.998"
    assert s.tname == "NZ_CP041014.1"


def test_map_sequences_trunc(mapper, zymo_reads):
    paf_dict = mapper.map_sequences(sequences=zymo_reads, trunc=True)
    assert len(paf_dict) == 8393
    s = paf_dict['ERR3152366.998'][0]
    assert type(s) is PafLine
    assert s.qname == "ERR3152366.998"
    assert s.tname == "NZ_CP041014.1"


