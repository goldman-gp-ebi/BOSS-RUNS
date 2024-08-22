import pytest
from pathlib import Path

import boss.mapper
from boss.paf import PafLine


@pytest.fixture
def dummy_fasta():
    fa = '../data/test.fasta'
    with open(fa, 'w') as f:
        f.write('>chr1\nAAAAAAAAGCTGACTATCGATCGTAGCTAGCTGACTGATCGTG\n'
                '>chr2\nGTGATCGTGATCGATCGTAGTCGTGTGCTAGTGTGTGTGTGATCGT\n')
    return fa


@pytest.fixture
def mapper(fasta_file):
    mapper = boss.mapper.Mapper(ref=fasta_file)
    # test init of non-default parameters
    _ = boss.mapper.Mapper(ref=fasta_file, default=False)
    return mapper


@pytest.mark.xfail(raises=FileNotFoundError)
def test_non_existing_ref():
    _ = boss.mapper.Mapper(ref="./dummy.fa")


def test_indexer(dummy_fasta):
    _ = boss.mapper.Indexer(fasta=dummy_fasta, mmi=f'{dummy_fasta}.mmi')
    assert Path(f'{dummy_fasta}.mmi').is_file()


def test_indexer_zymo(fasta_file):
    _ = boss.mapper.Indexer(fasta=fasta_file, mmi=f'{fasta_file}.mmi')
    assert Path(f'{fasta_file}.mmi').is_file()


def test_map_sequences(mapper, zymo_read_batch):
    paf_dict = mapper.map_sequences(sequences=zymo_read_batch.read_sequences)
    assert len(paf_dict) == 9108
    s = paf_dict['ERR3152366.998'][0]
    assert type(s) is PafLine
    assert s.qname == "ERR3152366.998"
    assert s.tname == "NZ_CP041014.1"


def test_map_sequences_trunc(mapper, zymo_read_batch):
    paf_dict = mapper.map_sequences(sequences=zymo_read_batch.read_sequences, trunc=True)
    assert len(paf_dict) == 8393
    s = paf_dict['ERR3152366.998'][0]
    assert type(s) is PafLine
    assert s.qname == "ERR3152366.998"
    assert s.tname == "NZ_CP041014.1"


