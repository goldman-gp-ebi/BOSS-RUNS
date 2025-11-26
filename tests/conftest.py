import glob
import pytest
from pathlib import Path

from boss.runs.reference import Reference
from boss.paf import Paf
from boss.sampler import Sampler
from boss.batch import ReadCache
from boss.batch import FastqBatch
from boss.config import Config
from boss.mapper import Indexer


from .constants import PATHS



@pytest.fixture
def zymo_index():
    mmip = f'{PATHS.fasta}.mmi'
    _ = Indexer(fasta=PATHS.fasta, mmi=mmip)
    assert Path(mmip).is_file()
    return mmip


@pytest.fixture
def zymo_ref(zymo_index):
    r = Reference(ref=PATHS.fasta, mmi=zymo_index)
    return r


@pytest.fixture
def paf_dict():
    return Paf.parse_PAF(PATHS.paf, min_len=1)


@pytest.fixture
def sampler():
    # initialise the wrapper class for the fastq and paf streaming
    s = Sampler(
            source=PATHS.fastq,
            paf_full=PATHS.paf,
            paf_trunc=PATHS.paf_trunc,
            maxbatch=10,
            batchsize=50,
    )
    return s


@pytest.fixture
def sampler_nopaf():
    # initialise the wrapper class for the fastq and paf streaming
    s = Sampler(
            source=PATHS.fastq,
            maxbatch=10,
            batchsize=50,
    )
    return s


@pytest.fixture
def read_cache():
    return ReadCache(batchsize=10, dumptime=10_000)


@pytest.fixture
def fq_file_mix():
    f = glob.glob(f'{PATHS.fastq_dir}/F*')
    f.sort()
    # add file without ch=
    f.append(glob.glob(f'{PATHS.fastq_dir_noch}/*_0.fq')[0])
    return f


@pytest.fixture
def zymo_read_batch():
    return FastqBatch(fq_files=[PATHS.fastq])


@pytest.fixture
def zymo_read_batch_big():
    return FastqBatch(fq_files=[f for f in glob.glob(f'{PATHS.fastq_dir}/F*')])


@pytest.fixture
def args_fake_real():
    conf = Config()
    # assign some args since we don't load the full config
    conf.args.toml_readfish = "TEST"
    conf.args.live_run = True
    # to wait less long
    conf.args.data_wait = 4
    return conf.args


