import time
import pytest
from pathlib import Path

import boss.live




def test_split_flowcell(data_loc):
    channels = boss.live.LiveRun.split_flowcell(out_path=data_loc, run_name="aeons")
    assert isinstance(channels, set)
    assert len(channels) == 256



@pytest.mark.xfail(raises=ValueError)
def test_split_flowcell_dummy(data_loc):
    channels = boss.live.LiveRun.split_flowcell(out_path=data_loc, run_name="dummy")
    assert isinstance(channels, set)
    assert len(channels) == 0



def test_connect_sequencer(data_loc):
    op = boss.live.LiveRun.connect_sequencer(device="TEST")
    assert op == data_loc



@pytest.mark.xfail(raises=ValueError)
def test_connect_sequencer_dummy():
    _ = boss.live.LiveRun.connect_sequencer(device="dummy")


# some dummy fixtures that allow me to specify all real paths in conftest
@pytest.fixture
def dummy_loc():
    return '../scripts'


@pytest.fixture
def empty_set():
    return set()


@pytest.fixture
def previous_files(fastq_pass_dir):
    return {f'{fastq_pass_dir}/FAT91932_pass_e7bf7751_f43c451e_1.fastq.gz',
            f'{fastq_pass_dir}/FAT91932_pass_e7bf7751_f43c451e_2.fastq.gz'}


# here we pass fixtures as parameters
@pytest.mark.parametrize("fastq_pass, processed_files, exp", [
    ("fastq_pass_dir", 'empty_set', 5),
    ("fastq_pass_dir", 'previous_files', 3),
    ("dummy_loc", 'empty_set', 0),
])
def test_scan_dir(fastq_pass, processed_files, exp, request):
    new_fq = boss.live.LiveRun.scan_dir(
        request.getfixturevalue(fastq_pass),
        request.getfixturevalue(processed_files)
    )
    assert len(new_fq) == exp
    if new_fq:
        assert Path(new_fq[0]).is_file()


def test_launch_readfish(readfish_toml_loc):
    log = boss.live.LiveRun.launch_readfish(toml=readfish_toml_loc, device="dummy", name="dummy")
    time.sleep(1)  # required for file to be registered, otherwise assertion fails
    assert Path(log).is_file()
    assert Path(log).stat().st_size > 0


def test_launch_readfish_TEST(readfish_toml_loc):
    log = boss.live.LiveRun.launch_readfish(toml=readfish_toml_loc, device="TEST", name="dummy")
    assert not log


@pytest.mark.xfail(raises=FileNotFoundError)
def test_launch_readfish_DUMMY():
    boss.live.LiveRun.launch_readfish(toml="dummy", device="dummy", name="dummy")


