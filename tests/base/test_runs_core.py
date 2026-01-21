import pytest
from pathlib import Path
import logging
import time
import numpy as np

import boss.runs.core
import boss.config

from ..constants import PATHS


@pytest.fixture
def args():
    conf = boss.config.Config()
    args = conf.args
    # assign some args since we don't load the full config
    args.toml_readfish = "TEST"
    args.split_flowcell = False
    args.live_run = True
    args.ref = PATHS.fasta
    args.mmi = PATHS.mmi
    return args


@pytest.fixture
def BossRuns(args):
    bossruns = boss.runs.core.BossRuns(args=args)
    return bossruns



def test_init(BossRuns):
    BossRuns.init()
    assert type(BossRuns.ref) is boss.runs.reference.Reference  # type: ignore
    assert Path(f"{BossRuns.ref.ref}.mmi").is_file()
    assert len(BossRuns.contigs) == 9
    logging.info(BossRuns.contigs.keys())
    assert BossRuns.contigs["NZ_CP041015.1"].length == 4045619



@pytest.mark.parametrize("modes", [True, False])
def test__init_dummy_strats(BossRuns, modes):
    BossRuns.init()
    strat_dict = BossRuns.ref.get_strategy_dict()
    assert len(strat_dict) == 9
    assert strat_dict["NZ_CP041015.1"].shape == (4045619 // 100, 2, 1)
    # test _write_contig_strategies (run during init)
    assert (Path(BossRuns.out_dir) / "masks" / "boss.npz").is_file()



def test_process_batch(BossRuns):
    BossRuns.init()
    BossRuns.launch_live_components()
    assert BossRuns.batch == 0
    tic = time.time()
    # we need to switch bucket switches manually here
    for cname, cont in BossRuns.contigs_filt.items():
        cont.switched_on = np.ones(shape=(len(BossRuns.args.barcodes)), dtype="bool")
    next_update = BossRuns.process_batch(BossRuns.process_batch_runs)
    assert BossRuns.batch == 1
    assert next_update != BossRuns.args.wait
    # check that new strats were produced
    assert Path("out_boss/masks/boss.npz").stat().st_mtime > tic



