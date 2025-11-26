from types import SimpleNamespace


DATA_BASE = 'data/BOSS_test_data'
CONF_BASE = 'tests/config'


PATHS = SimpleNamespace(
    data='../data',
    fastq=f"{DATA_BASE}/ERR3152366_10k.fq",
    fastq_gz=f"{DATA_BASE}/ERR3152366_10k.fq.gz",
    paf=f"{DATA_BASE}/ERR3152366_10k.paf",
    paf_trunc=f"{DATA_BASE}/ERR3152366_10k_trunc.paf",
    fasta=f"{DATA_BASE}/zymo.fa",
    mmi=f"{DATA_BASE}/zymo.fa.mmi",
    fastq_dir=f"{DATA_BASE}/fastq_pass",
    fastq_dir_noch=f"{DATA_BASE}/fastq_pass_no_ch/",
    readfish_toml=f"{CONF_BASE}/BOSS_RUNS_RF.toml",
    readfish_toml_singleregion=f"{CONF_BASE}/BOSS_RUNS_RF_singleregion.toml",
    readfish_toml_aeons=f"{CONF_BASE}/BOSS_AEONS_RF.toml",
)

ARG_DICT = {
    "boss-runs": ['--toml', f"{CONF_BASE}/BOSS_RUNS.toml", '--toml_readfish', f"{CONF_BASE}/BOSS_RUNS_RF.toml"],
    "boss-runs-sim": ['--toml', f"{CONF_BASE}/BOSS_RUNS_SIM.toml"],
    "boss-aeons": ['--toml', f"{CONF_BASE}/BOSS_AEONS.toml", '--toml_readfish', f"{CONF_BASE}/BOSS_AEONS_RF.toml"],
    "boss-aeons-sim": ['--toml', f"{CONF_BASE}/BOSS_AEONS_SIM.toml"],
    "boss-aeons-broken": ['--toml', f"{CONF_BASE}/BOSS_AEONS.toml", '--toml_readfish', f"{CONF_BASE}/BOSS_AEONS_RF_broken.toml"],
}

