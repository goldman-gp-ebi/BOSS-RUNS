from types import SimpleNamespace

import pytest

import boss.config



def test_defaults():
    # just testing some arbitrary values
    conf = boss.config.Config()
    assert type(conf.args) is SimpleNamespace
    assert conf.args.min_seq_len == 2500
    assert conf.args.fq == ""
    assert conf.args.device == "TEST"



def test_config_runs(arg_dict):
    conf = boss.config.Config(parse=True, arg_list=arg_dict['boss-runs'])
    assert conf.args.name == "runs"
    assert conf.args.device == "MS00000"
    assert conf.args.ref == "../data/zymo.fa"
    assert conf.args.live_run is True
    assert conf.args.sim_run is False
    assert conf.args.toml_readfish == "./config/BOSS_RUNS_RF.toml"
    assert conf.args.split_flowcell is True


def test_config_runs_sim(arg_dict):
    conf = boss.config.Config(parse=True, arg_list=arg_dict['boss-runs-sim'])
    assert conf.args.name == "runs"
    assert conf.args.device == "TEST"
    assert conf.args.live_run is False
    assert conf.args.sim_run is True
    assert hasattr(conf.args, "toml_readfish") is False



def test_config_aeons(arg_dict):
    conf = boss.config.Config(parse=True, arg_list=arg_dict['boss-aeons'])
    assert conf.args.name == "aeons"
    assert conf.args.device == "MS00000"
    assert conf.args.ref == ""
    assert conf.args.live_run is True
    assert conf.args.sim_run is False
    assert conf.args.toml_readfish == "./config/BOSS_AEONS_RF.toml"
    assert conf.args.split_flowcell is True
    assert conf.args.lowcov == 1


@pytest.mark.xfail(raises=ValueError)
def test_config_aeons_broken(arg_dict):
    _ = boss.config.Config(parse=True, arg_list=arg_dict['boss-aeons-broken'])


def test_config_aeons_sim(arg_dict):
    conf = boss.config.Config(parse=True, arg_list=arg_dict['boss-aeons-sim'])
    assert conf.args.name == "aeons"
    assert conf.args.device == "TEST"
    assert conf.args.live_run is False
    assert conf.args.sim_run is True
    assert hasattr(conf.args, "toml_readfish") is False




