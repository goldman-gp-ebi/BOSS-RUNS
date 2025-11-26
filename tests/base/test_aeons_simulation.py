import pytest
import time
from pathlib import Path
import subprocess

import boss.config
import boss.aeons.simulation

from ..constants import PATHS


@pytest.fixture
def args():
    conf = boss.config.Config()
    args = conf.args
    # assign some args since we don't load the full config
    args.live_run = False
    args.fq = PATHS.fastq
    args.maxb = 8
    args.batchsize = 100
    args.dumptime = 10000
    return args




@pytest.mark.xfail(raises=ValueError)
def test_init_fail(args):
    b = boss.aeons.simulation.BossAeonsSim(args=args)
    b.init_sim()



def test_init(args):
    args.batchsize = 1000
    b = boss.aeons.simulation.BossAeonsSim(args=args)
    b.init_sim()
    assert Path("00_reads/control_0.fa").is_file()
    assert Path("00_reads/boss_0.fa").is_file()
    assert len(b.strat) == 0
    assert Path(b.pool.contig_fa).is_file()
    subprocess.run("rm -r 00_reads/", shell=True)



def test_process_batch(args):
    tic = time.time()
    args.batchsize = 1000
    args.lowcov = 1
    b = boss.aeons.simulation.BossAeonsSim(args=args)
    b.init_sim()
    assert b.batch == 5
    # add some new data
    next_update = b.process_batch_sim(b.process_batch_aeons_sim)
    assert b.batch == 6
    assert next_update != b.args.wait
    # check that new contigs were produced
    assert Path("out_boss/contigs/aeons.fa").stat().st_mtime > tic
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic
    assert Path("out_boss/contigs/aeons.fa").read_text().startswith(">utg0")
    # run another batch to reject some stuff
    tic = time.time()
    next_update = b.process_batch_sim(b.process_batch_aeons_sim)
    assert b.batch == 7
    assert next_update != b.args.wait
    # check that new contigs were produced
    assert Path("out_boss/contigs/aeons.fa").stat().st_mtime > tic
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic
    assert Path("out_boss/contigs/aeons.fa").read_text().startswith(">utg0")
    assert b.read_cache.time_boss < b.read_cache.time_control
    subprocess.run("rm -r out_boss/", shell=True)




