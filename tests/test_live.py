import pytest
from pathlib import Path

import boss.live




def test_split_flowcell():
    channels = boss.live.LiveRun.split_flowcell(out_path="../data", run_name="aeons")
    assert isinstance(channels, set)
    assert len(channels) == 256



@pytest.mark.xfail(raises=ValueError)
def test_split_flowcell_dummy():
    channels = boss.live.LiveRun.split_flowcell(out_path="../data", run_name="dummy")
    assert isinstance(channels, set)
    assert len(channels) == 0



def test_connect_sequencer():
    op = boss.live.LiveRun.connect_sequencer(device="TEST")
    assert op == "../data"


# TODO implement a test with simulated device
# def test_connect_sequencer_sim():
#     op = boss.live.LiveRun.connect_sequencer(device="MS00000")
#     assert op == "./data"


@pytest.mark.xfail(raises=ValueError)
def test_connect_sequencer_dummy():
    _ = boss.live.LiveRun.connect_sequencer(device="dummy")



@pytest.mark.parametrize("fastq_pass, processed_files, exp", [
    ("../data/fastq_pass", set(), 11),
    ("../data/fastq_pass", {'../data/fastq_pass/ERR3152366_1k_3.fq', '../data/fastq_pass/ERR3152366_1k_7.fq'}, 9),
    ("../scripts", {}, 0),
])
def test_scan_dir(fastq_pass, processed_files, exp):
    new_fq = boss.live.LiveRun.scan_dir(fastq_pass, processed_files)
    assert len(new_fq) == exp
    if new_fq:
        assert Path(new_fq[0]).is_file()


def test_launch_readfish():
    boss.live.LiveRun.launch_readfish(toml="./config/BOSS_RUNS_RF.toml", device="dummy", name="dummy")

def test_launch_readfish_TEST():
    boss.live.LiveRun.launch_readfish(toml="./config/BOSS_RUNS_RF.toml", device="TEST", name="dummy")

@pytest.mark.xfail(raises=FileNotFoundError)
def test_launch_readfish_DUMMY():
    boss.live.LiveRun.launch_readfish(toml="dummy", device="dummy", name="dummy")


