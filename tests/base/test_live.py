import time
import pytest
from pathlib import Path

import boss.live

from ..constants import PATHS


@pytest.fixture
def dummy_sequencer():
    seq = boss.live.Sequencer()
    return seq



def test_grab_channels(dummy_sequencer):
    dummy_sequencer.grab_channels(run_name="aeons")
    channels = dummy_sequencer.channels
    assert isinstance(channels, set)
    assert len(channels) == 256



@pytest.mark.xfail(raises=ValueError)
def test_grab_channels_wrongname(dummy_sequencer):
    dummy_sequencer.grab_channels(run_name="dummy")
    channels = dummy_sequencer.channels
    assert isinstance(channels, set)
    assert len(channels) == 0




@pytest.mark.xfail(raises=ValueError)
def test_connect_sequencer_dummy():
    _ = boss.live.LiveRun.connect_sequencer(device="dummy")



# here we pass fixtures as parameters
@pytest.mark.parametrize("fastq_pass, processed_files, exp", [
    (PATHS.fastq_dir, set(), 5),
    (PATHS.fastq_dir, {f'{PATHS.fastq_dir}/FAT91932_pass_e7bf7751_f43c451e_1.fastq.gz', f'{PATHS.fastq_dir}/FAT91932_pass_e7bf7751_f43c451e_2.fastq.gz'}, 3),
    ("../scripts", set(), 0),
])
def test_scan_dir(fastq_pass, processed_files, exp):
    new_fq = boss.live.LiveRun.scan_dir(fastq_pass=fastq_pass, processed_files=processed_files)
    assert len(new_fq) == exp
    if new_fq:
        assert Path(new_fq[0]).is_file()


def test_launch_readfish():
    log = boss.live.LiveRun.launch_readfish(toml=PATHS.readfish_toml, device="dummy", name="dummy")
    time.sleep(1)  # required for file to be registered, otherwise assertion fails
    assert Path(log).is_file()
    assert Path(log).stat().st_size > 0


def test_launch_readfish_TEST():
    log = boss.live.LiveRun.launch_readfish(toml=PATHS.readfish_toml, device="TEST", name="dummy")
    assert not log


@pytest.mark.xfail(raises=FileNotFoundError)
def test_launch_readfish_DUMMY():
    boss.live.LiveRun.launch_readfish(toml="dummy", device="dummy", name="dummy")


