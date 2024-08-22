import subprocess
import time
import pytest
from pathlib import Path

from boss.aeons.core import BossAeons


@pytest.fixture
def boss_aeons_inst(args_fake_real):
    return BossAeons(args=args_fake_real)


def test_process_batch(boss_aeons_inst):
    tic = time.time()
    # empty dir at first
    subprocess.run("rm -r ../data/fastq_pass && mkdir ../data/fastq_pass", shell=True)
    # put the single 10k file in the dir
    subprocess.run("cp ../data/ERR3152366_10k.fq ../data/fastq_pass/", shell=True)
    boss_aeons_inst.launch_live_components()
    boss_aeons_inst.init()
    # remove that file again
    subprocess.run("rm -r ../data/fastq_pass && mkdir ../data/fastq_pass", shell=True)
    assert boss_aeons_inst.batch == 0
    # add some new data
    subprocess.run("cp ../data/BOSS_test_data/fastq_pass/FAT91932_pass_e7bf7751_f43c451e_0.fastq.gz ../data/fastq_pass/", shell=True)
    subprocess.run("cp ../data/ERR3152366_10k.fq ../data/fastq_pass", shell=True)
    next_update = boss_aeons_inst.process_batch(boss_aeons_inst.process_batch_aeons)
    subprocess.run("rm -r ../data/fastq_pass && mkdir ../data/fastq_pass", shell=True)
    subprocess.run("cp -r ../data/BOSS_test_data/fastq_pass/ ../data/", shell=True)
    assert boss_aeons_inst.batch == 1
    assert next_update != boss_aeons_inst.args.wait
    # check that new contigs were produced
    assert Path("out_boss/contigs/aeons.fa").stat().st_mtime > tic
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic
    assert Path("out_boss/contigs/aeons.fa").read_text().startswith(">utg000001l")



def test_cleanup(boss_aeons_inst):
    boss_aeons_inst.cleanup()
    assert len(list(Path("..").glob(f"{boss_aeons_inst.name}.*"))) == 0





