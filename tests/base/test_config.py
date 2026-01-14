from types import SimpleNamespace

import pytest

import boss.config

from ..constants import DATA_BASE, CONF_BASE, ARG_DICT



def test_defaults():
    # just testing some arbitrary values
    conf = boss.config.Config()
    assert type(conf.args) is SimpleNamespace
    assert conf.args.min_seq_len == 2500
    assert conf.args.fq == ""
    assert conf.args.device == "TEST"



def test_config_runs():
    conf = boss.config.Config(parse=True, arg_list=ARG_DICT['boss-runs'])
    assert conf.args.name == "runs"
    assert conf.args.device == "MS00000"
    assert conf.args.ref == f"{DATA_BASE}/zymo.fa"
    assert conf.args.live_run is True
    assert conf.args.sim_run is False
    assert conf.args.toml_readfish == f"{CONF_BASE}/BOSS_RUNS_RF.toml"


def test_config_runs_sim():
    conf = boss.config.Config(parse=True, arg_list=ARG_DICT['boss-runs-sim'])
    assert conf.args.name == "runs"
    assert conf.args.device == "TEST"
    assert conf.args.live_run is False
    assert conf.args.sim_run is True
    assert hasattr(conf.args, "toml_readfish") is False



def test_config_aeons():
    conf = boss.config.Config(parse=True, arg_list=ARG_DICT['boss-aeons'])
    assert conf.args.name == "aeons"
    assert conf.args.device == "MS00000"
    assert conf.args.ref == ""
    assert conf.args.live_run is True
    assert conf.args.sim_run is False
    assert conf.args.toml_readfish == f"{CONF_BASE}/BOSS_AEONS_RF.toml"
    assert conf.args.lowcov == 1


@pytest.mark.xfail(raises=ValueError)
def test_config_aeons_broken():
    _ = boss.config.Config(parse=True, arg_list=ARG_DICT['boss-aeons-broken'])


def test_config_aeons_sim():
    conf = boss.config.Config(parse=True, arg_list=ARG_DICT['boss-aeons-sim'])
    assert conf.args.name == "aeons"
    assert conf.args.device == "TEST"
    assert conf.args.live_run is False
    assert conf.args.sim_run is True
    assert hasattr(conf.args, "toml_readfish") is False




