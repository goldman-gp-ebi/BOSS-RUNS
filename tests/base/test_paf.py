import pytest
import logging
from io import StringIO

import boss.paf
from boss.utils import execute
from boss.dependencies import Dependencies
from boss.config import Config

from ..constants import PATHS


@pytest.fixture
def paf_io():
    with open(PATHS.paf, 'r') as f:
        paf_str = f.read()
    return StringIO(paf_str)


def test_parse_PAF():
    # init does not do anything for Paf
    _ = boss.paf.Paf()
    paf_dict = boss.paf.Paf.parse_PAF(PATHS.paf, min_len=1)
    assert len(paf_dict) == 9120


def test_parse_PAF_io(paf_io):
    paf_dict = boss.paf.Paf.parse_PAF(paf_io, min_len=1000)
    assert len(paf_dict) == 8913


def test_parse_PAF_dummy():
    paf_dict = boss.paf.Paf.parse_PAF("dummy", min_len=1)
    assert len(paf_dict) == 0



@pytest.fixture
def all_vs_all():
    # grab the minimap2 executable
    mm2 = Dependencies().minimap2
    # mm2 command for all-vs-all mapping
    comm = f'{mm2} -x ava-ont {PATHS.fastq} {PATHS.fastq} >{PATHS.fastq}.ava.paf'
    stdout, stderr = execute(comm)
    logging.info(stdout)
    logging.info(stderr)
    return f'{PATHS.fastq}.ava.paf'


def test_choose_best_mapper(paf_dict):
    recs = None
    for rid, recs in paf_dict.items():
        if len(recs) > 1:
            recs = boss.paf.Paf.choose_best_mapper(recs)
            break
    assert recs is not None
    logging.info(recs[0].__dict__)
    r = recs[0]
    assert r.qname == "ERR3152366.12"
    assert r.tname == "NZ_CP041014.1"


def test_parse_filter_classify_records(all_vs_all):
    records, skip = boss.paf.Paf.parse_filter_classify_records(
        all_vs_all,
        filters=Config().args
    )
    assert len(records) == 37014
    assert len(skip) == 6262
    classifications = {rec.c for rec in records}
    assert classifications == {2, 3, 4, 5, 6}  # no IM or TRIM




