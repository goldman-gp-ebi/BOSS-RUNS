from pathlib import Path
import pytest

import boss.aeons.repeats
from boss.aeons.sequences import SequencePool



@pytest.fixture
def repeat_filter(zymo_read_batch_big):
    seqpool = SequencePool(zymo_read_batch_big.read_sequences)
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


def test_filter_batch(repeat_filter, zymo_read_batch_big):
    filtered_seqs = repeat_filter.filter_batch(seq_dict=zymo_read_batch_big.read_sequences)
    assert len(zymo_read_batch_big.read_sequences) == 5000
    assert len(filtered_seqs) == 4988


def test_get_sequence(repeat, zymo_read_batch_big):
    seqpool = SequencePool(zymo_read_batch_big.read_sequences)
    repeat.get_sequence(seqpool=seqpool.sequences)


def test_fasta(repeat, zymo_read_batch_big):
    # test without grabbing sequence
    fa = repeat.fasta()
    assert fa == ''
    # assign new id and grab sequence
    repeat.rid = "b5ffbdab-29be-4afa-8346-2d5c4eba25ce"
    # logging.info(zymo_read_batch_big.read_sequences.keys())
    seqpool = SequencePool(zymo_read_batch_big.read_sequences)
    repeat.get_sequence(seqpool=seqpool.sequences)
    fa = repeat.fasta()
    assert fa.startswith('>b5ffbdab')

