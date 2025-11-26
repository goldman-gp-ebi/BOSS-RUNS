import os
import pytest
from pathlib import Path
import multiprocessing as mp
import time
import logging

from boss.BOSS import main as boss_main
from boss.readfish_boss import main as readfish_main

from boss.snippet import some_func



@pytest.fixture
def simion():
    return '/home/administrator/proj/simion/'


@pytest.fixture
def playback_dir(simion):
    paths = sorted(Path(f'{simion}/minknow_run/no_group/no_sample_id/').iterdir(), key=os.path.getmtime)
    return str(paths[-1])







def test_readfish_main(readfish_toml_loc, playback_dir):
    # launch the entry-point of readfish
    # this is similar to test_launch_readfish in test_live_playback.py
    # but uses main() instead of subprocess
    # advantage is that we can analyse coverage that way - coverage of dynamic_readfish.py and readfish_boss.py

    # here we expect half of all reads (all of boss_runs reads)
    # to be rejected, since they will not map to the randomly
    # generated reads from trigger_updates.py

    # remove all reads sequenced so far
    # otherwise readfish takes forever to start doing things
    assert '*' not in playback_dir
    fqs = (Path(playback_dir) / 'fastq_pass').glob('*.fastq.gz')
    for fq in fqs:
        fq.unlink()
    pods = (Path(playback_dir) / 'pod5').glob('*.pod5')
    for pod in pods:
        pod.unlink()


    arg_list = [readfish_toml_loc, 'MS00000', 'runs']
    # readfish_main(arg_list=arg_list)  # put this instead of somefunc

    rf = mp.Process(target=some_func, args=(arg_list, ))
    # rf = mp.Process(target=readfish_main, args=(arg_list,))

    rf.start()
    # let readfish run for a while
    time.sleep(10)
    rf.terminate()   # TODO continue here, this lets me run rf as a child process for as long as I like
    # time.sleep(4)
    rf.join()
    # rf.stop()
    #
    #
    # # get the log file of readfish - here the function does not return it
    # # so we have to glob for it instead
    # paths = sorted(Path(f'.').glob('*_readfish.log'), key=os.path.getmtime)
    # log = str(paths[-1])
    #
    # assert Path(log).is_file()
    # assert Path(log).stat().st_size > 0
    # assert Path('live_alignments.paf').is_file()
    # assert Path('live_reads.fq').is_file()
    # assert Path(f'{playback_dir}/unblocked_read_ids.txt').is_file()
    # assert Path(f'{playback_dir}/channels.toml').is_file()
    #
    # # # check that the logs contain what we are looking for
    # # # some ugly parsing
    # checked = [0] * 3
    # with open(log, 'r') as logf:
    #     ll = logf.readlines()
    #     for line in ll:
    #         lstrip = line.strip()
    #         if lstrip.startswith('strategies:'):
    #             lsplit = lstrip.split(':')[-1].strip()
    #             assert lsplit == 'out_runs/masks'
    #             checked[0] = 1
    #         if lstrip.startswith('contigs:'):
    #             lsplit = lstrip.split(':')[-1].strip()
    #             assert lsplit == 'out_runs/contigs'
    #             checked[1] = 1
    #         # grep for reloading occasions
    #         # TODO
    #
    #     # check the last line
    #     last = ll[-1]
    #     llast = last.split(';')
    #     seq = int(llast[2].split(':')[-1].replace(',', ''))
    #     unb = int(llast[3].split(':')[-1].replace(',', ''))
    #     assert unb > 100
    #     assert seq > 100
    #     assert seq / unb > 0.8
    #     assert seq / unb < 1.5
    #     checked[2] = 1
    # assert sum(checked) == 3



# TODO add a trigger updates call to the test above
# check simion readme for deets

# TODO add the same as above but for aeons
# add a call to trigger updates
# I think here we don't expect any rejections
# because reads won't map to anything so they will be sequenced fully


# TODO add tests for the boss_main entry point
# see what we need to grep for in the log files



# @pytest.mark.xfail(raises=ValueError, strict=True)
# def test_main_runs(arg_dict):
#     main(arg_dict["boss-runs"])
#
# @pytest.mark.xfail(raises=ValueError, strict=True)
# def test_main_aeons(arg_dict):
#     main(arg_dict["boss-aeons"])
#
# def test_main_runs_sim(arg_dict):
#     main(arg_dict["boss-runs-sim"])
#
# @pytest.mark.xfail(raises=ValueError, strict=True)
# def test_main_aeons_sim(arg_dict):
#     main(arg_dict["boss-aeons-sim"])



