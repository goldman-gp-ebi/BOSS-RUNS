import time
import pytest
import os
from pathlib import Path
import subprocess

import boss.live


@pytest.fixture
def simion():
    return '/home/administrator/proj/simion/'


@pytest.fixture
def playback_dir(simion):
    paths = sorted(Path(f'{simion}/minknow_run/no_group/no_sample_id/').iterdir(), key=os.path.getmtime)
    return str(paths[-1])



def test_launch_readfish(readfish_toml_loc, playback_dir):
    # this is useful to have as first test so that we get a channels toml directly from readfish
    log = boss.live.LiveRun.launch_readfish(toml=readfish_toml_loc, device="MS00000", name="runs")
    time.sleep(1)  # required for file to be registered, otherwise assertion fails
    assert Path(log).is_file()
    assert Path(log).stat().st_size > 0
    assert Path('live_alignments.paf').is_file()
    assert Path('live_reads.fq').is_file()
    assert Path(f'{playback_dir}/unblocked_read_ids.txt').is_file()
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

    # let readfish run for a while
    time.sleep(30)
    # stop readfish
    subprocess.run('pkill -f readfish', shell=True)
    # check that the logs contain what we are looking for
    # some ugly parsing
    checked = [0] * 3
    with open(log, 'r') as logf:
        ll = logf.readlines()
        for line in ll:
            l = line.strip()
            if l.startswith('strategies:'):
                lsplit = l.split(':')[-1].strip()
                assert lsplit == 'out_runs/masks'
                checked[0] = 1
            if l.startswith('contigs:'):
                lsplit = l.split(':')[-1].strip()
                assert lsplit == 'out_runs/contigs'
                checked[1] = 1
        # check the last line
        last = ll[-1]
        llast = last.split(';')
        seq = int(llast[2].split(':')[-1].replace(',', ''))
        unb = int(llast[3].split(':')[-1].replace(',', ''))
        assert unb > 100
        assert seq > 100
        assert seq / unb > 0.8
        assert seq / unb < 1.5
        checked[2] = 1
    assert sum(checked) == 3



def test_connect_sequencer(playback_dir):
    op = boss.live.LiveRun.connect_sequencer(device="MS00000")
    assert Path(op) == Path(playback_dir)



def test_scan_dir(playback_dir):
    op = boss.live.LiveRun.connect_sequencer(device="MS00000")
    assert Path(op) == Path(playback_dir)
    new_fq = boss.live.LiveRun.scan_dir(fastq_pass=str(Path(op) / 'fastq_pass'), processed_files=set())
    assert len(new_fq) > 0
    assert Path(new_fq[0]).is_file()



def test_split_flowcell(playback_dir):
    channels = boss.live.LiveRun.split_flowcell(out_path=playback_dir, run_name="runs")
    assert isinstance(channels, set)
    assert len(channels) == 256




