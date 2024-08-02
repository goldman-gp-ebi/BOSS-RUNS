from pathlib import Path
import pytest

import boss.aeons.repeats
import boss.batch
from boss.aeons.sequences import SequencePool


@pytest.fixture
def batch():
    f = list(Path("../data/fastq_pass/").glob("F*"))
    b = boss.batch.FastqBatch(f)
    return b


@pytest.fixture
def repeat_filter(batch):
    seqpool = SequencePool(batch.read_sequences)
    rf = boss.aeons.repeats.RepeatFilter(name="boss", seqpool=seqpool)
    return rf


@pytest.fixture
def repeat():
    return boss.aeons.repeats.Repeat(rid="dummy")


def test_init(repeat_filter):
    assert Path(f'{repeat_filter.name}.seqs.fa').is_file()
    assert len(repeat_filter.covs) == 726
    assert repeat_filter.lim == 6.0
    assert len(repeat_filter.repeats) > 20


def test_filter_batch(repeat_filter, batch):
    filtered_seqs = repeat_filter.filter_batch(seq_dict=batch.read_sequences)
    assert len(batch.read_sequences) == 5000
    assert len(filtered_seqs) == 4988


def test_get_sequence(repeat, batch):
    seqpool = SequencePool(batch.read_sequences)
    repeat.get_sequence(seqpool=seqpool.sequences)


def test_fasta(repeat, batch):
    # test without grabbing sequence
    fa = repeat.fasta()
    assert fa == ''
    # assign new id and grab sequence
    repeat.rid = "b5ffbdab-29be-4afa-8346-2d5c4eba25ce"
    # logging.info(batch.read_sequences.keys())
    seqpool = SequencePool(batch.read_sequences)
    repeat.get_sequence(seqpool=seqpool.sequences)
    fa = repeat.fasta()
    assert fa.startswith('>b5ffbdab')

