import pytest
import argparse
from pathlib import Path
import boss.readfish.boss_readfish as rf
import subprocess


"""
This test suite only works if a playback experiment is running
I.e. the code needs to connect to a running experiment to perform 
most of the functions in the readfish module
"""


ARG_LIST_RUNS = ["config/BOSS_RUNS_RF.toml", "MS00000", "bossruns"]
ARG_LIST_AEONS = ["config/BOSS_AEONS_RF.toml", "MS00000", "bossaeons"]


@pytest.fixture
def args_runs():
    _, args = rf.get_args(ARG_LIST_RUNS)
    return args

@pytest.fixture
def args_aeons():
    _, args = rf.get_args(ARG_LIST_AEONS)
    return args


def test_get_args(args_runs):
    assert isinstance(args_runs, argparse.Namespace)
    assert args_runs.host == "127.0.0.1"
    assert args_runs.toml == ARG_LIST_RUNS[0]


@pytest.fixture
def analysis_mod_runs(args_runs):
    """
    Create an AnalysisMod object for testing
    This test is modeled after run() in the boss_readfish.py entry-point
    """
    import logging
    from readfish._read_until_client import RUClient
    from readfish.read_until.read_cache import AccumulatingCache
    from readfish._utils import get_device
    from readfish._config import Conf

    args = args_runs

    # Setup logger used in this entry point, this one should be passed through
    logger = logging.getLogger(f"readfish.{args.command}")

    # Fetch sequencing device
    position = get_device(args.device, host=args.host, port=args.port)

    # Create a read until client
    read_until_client = RUClient(
        mk_host=position.host,
        mk_port=position.description.rpc_ports.secure,
        filter_strands=True,
        cache_type=AccumulatingCache,
        timeout=args.wait_for_ready,
    )

    # Load TOML configuration
    conf = Conf.from_file(args.toml, read_until_client.channel_count, logger=logger)

    worker = rf.AnalysisMod(
        read_until_client,
        conf=conf,
        logger=logger,
        debug_log=args.debug_log,
        unblock_duration=args.unblock_duration,
        throttle=args.throttle,
        dry_run=args.dry_run,
        toml=args.toml,
    )
    return worker


@pytest.fixture
def analysis_mod_aeons(args_aeons):
    """
    Create an AnalysisMod object for testing
    This test is modeled after run() in the boss_readfish.py entry-point
    """
    import logging
    from readfish._read_until_client import RUClient
    from readfish.read_until.read_cache import AccumulatingCache
    from readfish._utils import get_device
    from readfish._config import Conf


    rf.BossBits.gen_dummy_idx()

    args = args_aeons


    # Setup logger used in this entry point, this one should be passed through
    logger = logging.getLogger(f"readfish.{args.command}")

    # Fetch sequencing device
    position = get_device(args.device, host=args.host, port=args.port)

    # Create a read until client
    read_until_client = RUClient(
        mk_host=position.host,
        mk_port=position.description.rpc_ports.secure,
        filter_strands=True,
        cache_type=AccumulatingCache,
        timeout=args.wait_for_ready,
    )

    # Load TOML configuration
    conf = Conf.from_file(args.toml, read_until_client.channel_count, logger=logger)

    worker = rf.AnalysisMod(
        read_until_client,
        conf=conf,
        logger=logger,
        debug_log=args.debug_log,
        unblock_duration=args.unblock_duration,
        throttle=args.throttle,
        dry_run=args.dry_run,
        toml=args.toml,
    )
    return worker



@pytest.fixture
def bossbits_runs(analysis_mod_runs):
    am = analysis_mod_runs
    subprocess.run("rm -r out_runs", shell=True)
    bb = rf.BossBits(conf=am.conf, logger=am.logger, mapper=am.mapper)
    return bb


@pytest.fixture
def bossbits_aeons(analysis_mod_aeons):
    am = analysis_mod_aeons
    subprocess.run("rm -r out_aeons", shell=True)
    bb = rf.BossBits(conf=am.conf, logger=am.logger, mapper=am.mapper)
    return bb



def test_bb_init_runs(bossbits_runs):
    bb = bossbits_runs
    assert len(bb.strand_converter) == 2
    assert isinstance(bb.mask_path, Path)
    assert isinstance(bb.cont_path, Path)
    assert (bb.mask_path / "boss.npz").is_file()


def test_bb_init_aeons(bossbits_aeons):
    bb = bossbits_aeons
    assert len(bb.strand_converter) == 2
    assert isinstance(bb.mask_path, Path)
    assert isinstance(bb.cont_path, Path)
    assert (bb.mask_path / "boss.npz").is_file()
    assert Path("readfish.mmi").is_file()


def test_reload_aeons(analysis_mod_aeons, bossbits_aeons):
    am = analysis_mod_aeons
    bb = bossbits_aeons
    mask_mtime_0 = bb.last_mask_mtime
    cont_mtime_0 = bb.last_contig_mtime
    # update both masks and contigs
    subprocess.run(f"cp ../data/boss.npz {bb.mask_path}", shell=True)
    subprocess.run(f"cp ../data/aeons.fa {bb.cont_path}", shell=True)
    reload_mapper, mapper = bb.reload(conf=am.conf)
    mask_mtime_1 = bb.last_mask_mtime
    cont_mtime_1 = bb.last_contig_mtime
    assert reload_mapper
    assert mapper
    # reload again straight away without new files
    reload_mapper, mapper = bb.reload(conf=am.conf)
    mask_mtime_2 = bb.last_mask_mtime
    cont_mtime_2 = bb.last_contig_mtime
    assert not reload_mapper
    assert mask_mtime_0 < mask_mtime_1
    assert mask_mtime_1 == mask_mtime_2
    assert cont_mtime_0 < cont_mtime_1
    assert cont_mtime_1 == cont_mtime_2



def test_reload_runs(analysis_mod_runs, bossbits_runs):
    am = analysis_mod_runs
    bb = bossbits_runs
    mask_mtime_0 = bb.last_mask_mtime
    cont_mtime_0 = bb.last_contig_mtime
    # update both masks and contigs
    subprocess.run(f"cp ../data/boss.npz {bb.mask_path}", shell=True)
    reload_mapper, mapper = bb.reload(conf=am.conf)
    mask_mtime_1 = bb.last_mask_mtime
    cont_mtime_1 = bb.last_contig_mtime
    assert not reload_mapper
    # reload again straight away without new files
    reload_mapper, mapper = bb.reload(conf=am.conf)
    mask_mtime_2 = bb.last_mask_mtime
    cont_mtime_2 = bb.last_contig_mtime
    assert not reload_mapper
    assert mask_mtime_0 < mask_mtime_1
    assert mask_mtime_1 == mask_mtime_2
    assert cont_mtime_0 == cont_mtime_1 == cont_mtime_2




@pytest.mark.timeout(15)
@pytest.mark.xfail()
def test_main_aeons():
    rf.main(arg_list=ARG_LIST_AEONS)


@pytest.mark.timeout(15)
@pytest.mark.xfail()
def test_main_runs():
    rf.main(arg_list=ARG_LIST_RUNS)

