import sys

import pytest

import boss.config

from ..constants import DATA_BASE, CONF_BASE



def test_defaults():
    # just testing some arbitrary values
    conf = boss.config.Config()
    assert type(conf.args) is boss.config.BossConfig
    assert conf.args.optional.min_seq_len == 2500
    assert conf.args.simulation.fq is None
    assert conf.args.live.device is None



def test_config_runs(monkeypatch):
    monkeypatch.setattr(sys, "argv", [sys.argv[0], "--toml", f"{CONF_BASE}/BOSS_RUNS.toml"])
    conf = boss.config.Config(parse=True)
    assert conf.args.general.name == "runs"
    assert conf.args.live.device == "MS00000"
    assert conf.args.general.ref == f"{DATA_BASE}/zymo.fa"
    assert conf.args.general.toml_readfish == f"{CONF_BASE}/BOSS_RUNS_RF.toml"
    assert conf.args.simulation.fq is None


def test_config_runs_sim(monkeypatch):
    monkeypatch.setattr(sys, "argv", [sys.argv[0], "--toml", f"{CONF_BASE}/BOSS_RUNS_SIM.toml"])
    conf = boss.config.Config(parse=True)
    assert conf.args.general.name == "runs"
    assert conf.args.live.device is None
    assert conf.args.general.toml_readfish is None




def test_config_aeons(monkeypatch):
    monkeypatch.setattr(sys, "argv", [sys.argv[0], "--toml", f"{CONF_BASE}/BOSS_AEONS.toml"])
    conf = boss.config.Config(parse=True)
    assert conf.args.general.name == "aeons"
    assert conf.args.live.device == "MS00000"
    assert conf.args.general.ref is None
    assert conf.args.general.toml_readfish == f"{CONF_BASE}/BOSS_AEONS_RF.toml"
    assert conf.args.optional.lowcov == 1


@pytest.mark.xfail(raises=ValueError)
def test_config_aeons_broken(monkeypatch):
    monkeypatch.setattr(sys, "argv", [sys.argv[0], "--toml", f"{CONF_BASE}/BOSS_AEONS_broken.toml"])
    _ = boss.config.Config(parse=True)


def test_config_aeons_sim(monkeypatch):
    monkeypatch.setattr(sys, "argv", [sys.argv[0], "--toml", f"{CONF_BASE}/BOSS_AEONS_SIM.toml"])
    conf = boss.config.Config(parse=True)
    assert conf.args.general.name == "aeons"
    assert conf.args.live.device is None
    assert conf.args.general.toml_readfish is None





