import time
import pytest
import os
from pathlib import Path
import subprocess
import logging

import boss.live


@pytest.fixture
def simion():
    return '/home/administrator/proj/simion/'


@pytest.fixture
def playback_dir(simion):
    paths = sorted(Path(f'{simion}/minknow_run/no_group/no_sample_id/').iterdir(), key=os.path.getmtime)
    return str(paths[-1])


@pytest.fixture
def readfish_log(readfish_toml_loc, playback_dir):
    # this is useful to have as first test so that we get a channels toml directly from readfish
    log = boss.live.LiveRun.launch_readfish(toml=readfish_toml_loc, device="MS00000", name="runs")

    # remove all reads sequenced so far
    # otherwise readfish takes forever to start doing things
    assert '*' not in playback_dir
    fqs = (Path(playback_dir) / 'fastq_pass').glob('*.fastq.gz')
    for fq in fqs:
        fq.unlink()
    pods = (Path(playback_dir) / 'pod5').glob('*.pod5')
    for pod in pods:
        pod.unlink()

    logging.info(log)
    time.sleep(40)
    subprocess.run('pkill -f readfish', shell=True)
    return log


def test_launch_readfish(readfish_log, playback_dir):
    log = readfish_log
    assert Path(log).is_file()
    assert Path(log).stat().st_size > 0
    assert Path('live_alignments.paf').is_file()
    assert Path('live_reads.fq').is_file()
    assert Path(f'{playback_dir}/unblocked_read_ids.txt').is_file()
    assert Path(f'{playback_dir}/channels.toml').is_file()

    # check that the logs contain what we are looking for
    # some ugly parsing
    checked = [0] * 3
    with open(log, 'r') as logf:
        ll = logf.readlines()
        for line in ll:
            lstrip = line.strip()
            if lstrip.startswith('strategies:'):
                lsplit = lstrip.split(':')[-1].strip()
                assert lsplit == 'out_runs/masks'
                checked[0] = 1
            if lstrip.startswith('contigs:'):
                lsplit = lstrip.split(':')[-1].strip()
                assert lsplit == 'out_runs/contigs'
                checked[1] = 1
        # check the last line
        last = ll[-1]
        logging.info(last)
        llast = last.split(';')
        seq = int(llast[2].split(':')[-1].replace(',', ''))
        unb = int(llast[3].split(':')[-1].replace(',', ''))
        assert unb > 100
        assert seq > 100
        # ratio of roughly around 1 is expected
        assert seq / unb > 0.8
        assert seq / unb < 1.5
        checked[2] = 1
    assert sum(checked) == 3



def test_connect_sequencer(playback_dir):
    sequencer = boss.live.LiveRun.connect_sequencer(device="MS00000")
    assert Path(sequencer.out_path) == Path(playback_dir)
    assert sequencer.device_type == 'min'



def test_scan_dir(playback_dir):
    sequencer = boss.live.LiveRun.connect_sequencer(device="MS00000")
    assert Path(sequencer.out_path) == Path(playback_dir)
    new_fqs = boss.live.LiveRun.scan_dir(fastq_pass=str(Path(sequencer.out_path) / 'fastq_pass'), processed_files=set())
    assert isinstance(new_fqs, list)
    # not testing the content here on purpose since I'm removing files



def test_grab_channels(playback_dir):
    sequencer = boss.live.LiveRun.connect_sequencer(device="MS00000")
    sequencer.grab_channels(run_name="runs")
    channels = sequencer.channels
    assert isinstance(channels, set)
    assert len(channels) == 256




def test_launch_readfish_singleregion(readfish_toml_loc_singleregion, playback_dir):
    # this is run last so the new channels.toml does not mess with the other tests
    log = boss.live.LiveRun.launch_readfish(toml=readfish_toml_loc_singleregion, device="MS00000", name="runs")
    time.sleep(2)  # required for file to be registered, otherwise assertion fails
    assert Path(log).is_file()
    assert Path(log).stat().st_size > 0
    assert Path(f'{playback_dir}/channels.toml').is_file()
    assert '*' not in playback_dir
    # remove all reads sequenced so far
    # otherwise readfish takes forever to start doing things
    fqs = (Path(playback_dir) / 'fastq_pass').glob('*.fastq.gz')
    for fq in fqs:
        fq.unlink()
    pods = (Path(playback_dir) / 'pod5').glob('*.pod5')
    for pod in pods:
        pod.unlink()

    time.sleep(5)
    subprocess.run('pkill -f readfish', shell=True)

    # check that the channels file contains only one condition
    # and that boss uses all channels accordingly
    sequencer = boss.live.LiveRun.connect_sequencer(device="MS00000")
    sequencer.grab_channels(run_name="runs")
    channels = sequencer.channels
    assert isinstance(channels, set)
    assert len(channels) == 0



