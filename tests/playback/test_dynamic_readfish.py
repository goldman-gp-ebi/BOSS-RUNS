import pytest
import argparse
from pathlib import Path
import subprocess

from readfish.plugins.abc import AlignerABC

from boss.readfish_boss import run
from boss.dynamic_readfish import BossBits, get_args



# the next two fixtures grab what we need to test the bossbits code,
# i.e. only the conf, logger, and mapper objects that are created
# during the setup of the readfish entry-point
@pytest.fixture
def conf_runs(readfish_toml_loc):
    parser, args = get_args(arg_list=[readfish_toml_loc, 'MS00000', 'runs'])
    args.return_conf = True
    conf, logger = run(parser=parser, args=args, extras=[])
    mapper: AlignerABC = conf.mapper_settings.load_object("Aligner")
    return conf, logger, mapper


@pytest.fixture
def conf_aeons(readfish_toml_aeons_loc):
    parser, args = get_args(arg_list=[readfish_toml_aeons_loc, 'MS00000', 'aeons'])
    args.return_conf = True
    conf, logger = run(parser=parser, args=args, extras=[])
    mapper: AlignerABC = conf.mapper_settings.load_object("Aligner")
    return conf, logger, mapper


# testing the argument generation for the readfish entry-point
def test_get_args(readfish_toml_loc):
    parser, args = get_args(arg_list=[readfish_toml_loc, 'MS00000', 'runs'])
    assert isinstance(args, argparse.Namespace)
    assert args.host == "127.0.0.1"
    assert args.toml == readfish_toml_loc


# testing the idx generation and keeping it around for the other tests
@pytest.mark.keepfiles
def test_gen_dummy_idx():
    BossBits.gen_dummy_idx()
    assert Path("readfish.mmi").is_file()


# here is where the bits from readfish.targets.Analysis are used
# these are required to initialise the bossbits
@pytest.fixture
def bossbits_runs(conf_runs):
    conf, logger, mapper = conf_runs
    subprocess.run("rm -r out_runs", shell=True)
    bb = BossBits(conf=conf, logger=logger, mapper=mapper)
    return bb


@pytest.fixture
def bossbits_aeons(conf_aeons):
    conf, logger, mapper = conf_aeons
    subprocess.run("rm -r out_aeons", shell=True)
    bb = BossBits(conf=conf, logger=logger, mapper=mapper)
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


# these two tests copy around some files and then run the reloading
# then they reload again to tests what happens when the files have
# not been modified since the last update
def test_reload_aeons(conf_aeons, bossbits_aeons):
    conf, _, _ = conf_aeons
    bb = bossbits_aeons
    mask_mtime_0 = bb.last_mask_mtime
    cont_mtime_0 = bb.last_contig_mtime
    # update both masks and contigs
    subprocess.run(f"cp ../data/boss.npz {bb.mask_path}", shell=True)
    subprocess.run(f"cp ../data/aeons.fa {bb.cont_path}", shell=True)
    reload_mapper, mapper = bb.reload(conf=conf)
    mask_mtime_1 = bb.last_mask_mtime
    cont_mtime_1 = bb.last_contig_mtime
    assert reload_mapper
    assert mapper
    # reload again straight away without new files
    reload_mapper, mapper = bb.reload(conf=conf)
    mask_mtime_2 = bb.last_mask_mtime
    cont_mtime_2 = bb.last_contig_mtime
    assert not reload_mapper
    assert mask_mtime_0 < mask_mtime_1
    assert mask_mtime_1 == mask_mtime_2
    assert cont_mtime_0 < cont_mtime_1
    assert cont_mtime_1 == cont_mtime_2


def test_reload_runs(conf_runs, bossbits_runs):
    conf, _, _ = conf_runs
    bb = bossbits_runs
    mask_mtime_0 = bb.last_mask_mtime
    cont_mtime_0 = bb.last_contig_mtime
    # update both masks and contigs
    subprocess.run(f"cp ../data/boss.npz {bb.mask_path}", shell=True)
    reload_mapper, mapper = bb.reload(conf=conf)
    mask_mtime_1 = bb.last_mask_mtime
    cont_mtime_1 = bb.last_contig_mtime
    assert not reload_mapper
    # reload again straight away without new files
    reload_mapper, mapper = bb.reload(conf=conf)
    mask_mtime_2 = bb.last_mask_mtime
    cont_mtime_2 = bb.last_contig_mtime
    assert not reload_mapper
    assert mask_mtime_0 < mask_mtime_1
    assert mask_mtime_1 == mask_mtime_2
    assert cont_mtime_0 == cont_mtime_1 == cont_mtime_2


