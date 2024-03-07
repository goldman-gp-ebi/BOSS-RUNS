import pytest
import subprocess
import time
from pathlib import Path

import boss.config
import boss.aeons.core


@pytest.fixture
def args():
    conf = boss.config.Config()
    # assign some args since we don't load the full config
    conf.args.toml_readfish = "TEST"
    conf.args.split_flowcell = False
    conf.args.live_run = True
    # to wait less long
    conf.args.data_wait = 4
    return conf.args


@pytest.fixture
def BossAeons(args):
    return boss.aeons.core.BossAeons(args=args)



def test_process_batch(BossAeons):
    tic = time.time()
    subprocess.run("cp ../data/ERR3152366_10k.fq ../data/fastq_pass/ERR3152366_10k.fq", shell=True)
    BossAeons.init()
    subprocess.run("rm ../data/fastq_pass/ERR3152366_10k.fq", shell=True)
    assert BossAeons.batch == 0
    # add some new data
    subprocess.run("rm ../data/fastq_pass/FAT91932_pass_e7bf7751_f43c451e_0.fastq.gz", shell=True)
    subprocess.run("rm ../data/fastq_pass/ERR3152366_10k.fq", shell=True)
    subprocess.run("cp ../data/fastq_pass_ch/FAT91932_pass_e7bf7751_f43c451e_0.fastq.gz ../data/fastq_pass", shell=True)
    subprocess.run("cp ../data/ERR3152366_10k.fq ../data/fastq_pass", shell=True)
    next_update = BossAeons.process_batch(BossAeons.process_batch_aeons)
    subprocess.run("rm ../data/fastq_pass/FAT91932_pass_e7bf7751_f43c451e_0.fastq.gz", shell=True)
    subprocess.run("rm ../data/fastq_pass/ERR3152366_10k.fq", shell=True)
    assert BossAeons.batch == 1
    assert next_update != BossAeons.args.wait
    # check that new contigs were produced
    assert Path("out_boss/contigs/aeons.fa").stat().st_mtime > tic
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic
    assert Path("out_boss/contigs/aeons.fa").read_text().startswith(">utg000001l")



def test_cleanup(BossAeons):
    BossAeons.cleanup()
    assert len(list(Path(".").glob(f"{BossAeons.name}.*"))) == 0





