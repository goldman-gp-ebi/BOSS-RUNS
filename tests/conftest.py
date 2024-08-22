import glob
from types import SimpleNamespace
import pytest
from pathlib import Path

from boss.runs.reference import Reference
from boss.paf import Paf
from boss.sampler import Sampler
from boss.batch import ReadCache
from boss.batch import FastqBatch
from boss.config import Config



paths = SimpleNamespace(
    data='../data',
    fastq="../data/ERR3152366_10k.fq",
    fastq_gz="../data/ERR3152366_10k.fq.gz",
    paf="../data/ERR3152366_10k.paf",
    paf_trunc="../data/ERR3152366_10k_trunc.paf",
    fasta="../data/zymo.fa",
    mmi='../data/zymo.fa.mmi',
    fastq_dir="../data/fastq_pass",
    fastq_dir_noch="../data/fastq_pass_no_ch/",
    readfish_toml='./config/BOSS_RUNS_RF.toml',
)

arg_dicts = {
    "boss-runs": ['--toml', "./config/BOSS_RUNS.toml", '--toml_readfish', "./config/BOSS_RUNS_RF.toml"],
    "boss-runs-sim": ['--toml', "./config/BOSS_RUNS_SIM.toml"],
    "boss-aeons": ['--toml', "./config/BOSS_AEONS.toml", '--toml_readfish', "./config/BOSS_AEONS_RF.toml"],
    "boss-aeons-sim": ['--toml', "./config/BOSS_AEONS_SIM.toml"],
    "boss-aeons-broken": ['--toml', "./config/BOSS_AEONS.toml", '--toml_readfish', "./config/BOSS_AEONS_RF_broken.toml"],
}





@pytest.fixture(autouse=True)
def check_for_new_files(request):
    '''
    this fixture runs with all tests. It checks which files
    and directories are created by the tests, via comparing
    the code before and after the yield.
    :param request:
    :return:
    '''
    p = Path('.')
    before_files = set([x for x in p.glob('**/*') if x.is_file()])
    before_dirs = set([x for x in p.glob('**/*') if x.is_dir()])
    yield
    after_files = set([x for x in p.glob('**/*') if x.is_file()])
    after_dirs = set([x for x in p.glob('**/*') if x.is_dir()])
    new_files = [f for f in after_files if f not in before_files]
    new_dirs = [f for f in after_dirs if f not in before_dirs]
    if new_files:
        print("\n{} created by {}".format(new_files, request.node))
        for f in new_files:
            f.unlink()
    if new_dirs:
        print("\n{} created by {}".format(new_dirs, request.node))
        for d in new_dirs:
            # for root, dirs, _ in d.walk(top_down=False):   # walk requires pathlib 3.12
            #     for name in dirs:
            #         (root / name).rmdir()
            try:
                d.rmdir()
            except OSError:
                pass


@pytest.fixture
def zymo_ref():
    r = Reference(ref=paths.fasta, mmi=paths.mmi)
    return r

@pytest.fixture
def paf_file():
    return paths.paf

@pytest.fixture
def paf_file_trunc():
    return paths.paf_trunc

@pytest.fixture
def paf_dict():
    return Paf.parse_PAF(paths.paf, min_len=1)

@pytest.fixture
def sampler():
    # initialise the wrapper class for the fastq and paf streaming
    s = Sampler(
            source=paths.fastq,
            paf_full=paths.paf,
            paf_trunc=paths.paf_trunc,
            maxbatch=10,
            batchsize=50,
    )
    return s


@pytest.fixture
def sampler_nopaf():
    # initialise the wrapper class for the fastq and paf streaming
    s = Sampler(
            source=paths.fastq,
            maxbatch=10,
            batchsize=50,
    )
    return s


@pytest.fixture
def read_cache():
    return ReadCache(batchsize=10, dumptime=10_000)


@pytest.fixture
def fastq_file():
    return paths.fastq


@pytest.fixture
def fastq_file_gz():
    return paths.fastq_gz

@pytest.fixture
def fasta_file():
    return paths.fasta

@pytest.fixture
def mmi_file():
    return paths.mmi

@pytest.fixture
def fastq_pass_dir():
    return paths.fastq_dir


@pytest.fixture
def fq_file_mix():
    f = glob.glob(f'{paths.fastq_dir}/F*')
    f.sort()
    # add file without ch=
    f.append(glob.glob(f'{paths.fastq_dir_noch}/*_0.fq')[0])
    return f



@pytest.fixture
def zymo_read_batch():
    return FastqBatch(fq_files=[paths.fastq])


@pytest.fixture
def zymo_read_batch_big():
    return FastqBatch(fq_files=[f for f in glob.glob(f'{paths.fastq_dir}/F*')])


@pytest.fixture
def arg_dict():
    return arg_dicts



@pytest.fixture
def args_fake_real():
    conf = Config()
    # assign some args since we don't load the full config
    conf.args.toml_readfish = "TEST"
    conf.args.split_flowcell = False
    conf.args.live_run = True
    # to wait less long
    conf.args.data_wait = 4
    return conf.args


@pytest.fixture
def readfish_toml_loc():
    return paths.readfish_toml


@pytest.fixture
def data_loc():
    return paths.data


