import pytest

from BOSS import main


arg_lists = {
    "boss-runs": ['--toml', "./config/BOSS_RUNS.toml", '--toml_readfish', "./config/BOSS_RUNS_RF.toml"],
    "boss-runs-sim": ['--toml', "./config/BOSS_RUNS_SIM.toml"],
    "boss-aeons": ['--toml', "./config/BOSS_AEONS.toml", '--toml_readfish', "./config/BOSS_AEONS_RF.toml"],
    "boss-aeons-sim": ['--toml', "./config/BOSS_AEONS_SIM.toml"],
}


@pytest.mark.xfail(raises=ValueError, strict=True)
def test_main_runs():
    main(arg_lists["boss-runs"])

@pytest.mark.xfail(raises=ValueError, strict=True)
def test_main_aeons():
    main(arg_lists["boss-aeons"])

def test_main_runs_sim():
    main(arg_lists["boss-runs-sim"])

@pytest.mark.xfail(raises=ValueError, strict=True)
def test_main_aeons_sim():
    main(arg_lists["boss-aeons-sim"])

