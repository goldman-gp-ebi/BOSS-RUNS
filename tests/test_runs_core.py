import pytest
from pathlib import Path
import logging
import time

import boss.runs.core
import boss.config



@pytest.fixture
def args():
    conf = boss.config.Config()
    args = conf.args
    # assign some args since we don't load the full config
    args.toml_readfish = "TEST"
    args.split_flowcell = False
    args.live_run = True
    args.ref = "../data/zymo.fa"
    return args


@pytest.fixture
def BossRuns(args):
    bossruns = boss.runs.core.BossRuns(args=args)
    return bossruns



def test_init(BossRuns):
    BossRuns.init()
    assert type(BossRuns.ref) is boss.runs.reference.Reference
    assert Path(f"{BossRuns.ref.ref}.mmi").is_file()
    assert len(BossRuns.contigs) == 9
    logging.info(BossRuns.contigs.keys())
    assert BossRuns.contigs["NZ_CP041015.1"].length == 4045619



@pytest.mark.parametrize("modes", [True, False])
def test__init_dummy_strats(BossRuns, modes):
    BossRuns.init()
    strat_dict = BossRuns.ref.get_strategy_dict()
    assert len(strat_dict) == 9
    assert strat_dict["NZ_CP041015.1"].shape == (4045619 // 100, 2)
    # test _write_contig_strategies (run during init)
    assert (Path(BossRuns.out_dir) / "masks" / "boss.npz").is_file()



def test_process_batch(BossRuns):
    BossRuns.init()
    assert BossRuns.batch == 0
    tic = time.time()
    next_update = BossRuns.process_batch(BossRuns.process_batch_runs)
    assert BossRuns.batch == 1
    assert next_update != BossRuns.args.wait
    # check that new strats were produced
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic



