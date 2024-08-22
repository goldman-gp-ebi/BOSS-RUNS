import pytest
from pathlib import Path
import logging
import time

import boss.runs.simulation
import boss.config



@pytest.fixture
def args(fasta_file, fastq_file, paf_file, paf_file_trunc):
    conf = boss.config.Config()
    args = conf.args
    # assign some args since we don't load the full config
    args.live_run = False
    args.ref = fasta_file
    args.fq = fastq_file
    args.paf_full = paf_file
    args.paf_trunc = paf_file_trunc
    args.maxb = 8
    args.batchsize = 100
    args.dumptime = 10000
    return args



def test_init(args):
    b = boss.runs.simulation.BossRunsSim(args=args)
    b.init_sim()
    assert type(b.ref) is boss.runs.reference.Reference
    assert Path(f"{b.ref.ref}.mmi").is_file()
    assert len(b.contigs) == 9
    logging.info(b.contigs.keys())
    assert b.contigs["NZ_CP041015.1"].length == 4045619
    assert Path("00_reads/control_0.fa").is_file()
    assert Path("00_reads/boss_0.fa").is_file()



def test_process_batch(args):
    args.batchsize = 500
    args.maxb = 9
    b = boss.runs.simulation.BossRunsSim(args=args)
    b.init_sim()
    assert b.batch == 0
    tic = time.time()
    # we need to switch bucket switches manually here
    for cname, cont in b.contigs_filt.items():
        cont.switched_on = True
    next_update = b.process_batch_sim(b.process_batch_runs_sim)
    assert b.batch == 1
    assert next_update != b.args.wait
    # check that new strats were produced
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic
    assert Path("00_reads/control_1.fa").is_file()
    assert Path("00_reads/boss_1.fa").is_file()
    # run another batch to check rejections
    tic = time.time()
    next_update = b.process_batch_sim(b.process_batch_runs_sim)
    assert b.batch == 2
    assert next_update != b.args.wait
    # check that new strats were produced
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic
    assert Path("00_reads/control_2.fa").is_file()
    assert Path("00_reads/boss_2.fa").is_file()
    assert b.read_cache.time_boss < b.read_cache.time_control


