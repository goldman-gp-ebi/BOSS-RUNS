# BOSS* development notes

## General organisation


BOSS* has two different running modes: either using a given reference or creating its on via real-time contig construction.
The arguments given in the toml file determine which mode is run (providing a reference executes BOSS-RUNS).

The code is organised as a module (boss) and submodules (runs and aeons):

- boss (shared codebase of BOSS*)
  - runs (code specific to *RUNS)
  - aeons (code specific to *AEONS)
  - readfish (modified entry point to readfish for BOSS*)


The main classes that orchestrate an experiment are organised via sub-classing.
The Boss class contains shared functionality, the rest is organised similar to the modules above, 
with separate classes that contain code for simulations.


classes:

|- Boss               |   boss/core.py              |
|-- BossRuns          |   boss/runs/core.py         |
|--- BossRunsSim      |   boss/runs/simulations.py  |
|-- BossAeons         |   boss/aeons/core.py        |
|--- BossAeonsSim     |   boss/aeons/simulations.py |


## Entry points


General entry-point for the package is `main()` in `boss.BOSS`. The hatch build copies this as a script into a callable path.

The superclass Boss initialises and then runs `process_batch()` for each batch of new data. (Or `process_batch_sim()` for simulations)

The `process_batch()` functions take a method of one of the subclasses and execute that as their main body





## Testing

The code is covered by pytests. Test data is available in a submodule in `data/`. 
To get the test data: `git clone --recurse-submodules ...` 


```shell
pip install pytest pytest-cov pytest-timeout 
pip install -e .
```

Most tests can be run stand-alone with

```shell
pytest tests/base/ --cov --cov-report html 
```

Some tests need to have playback sequencing running in the background.

```shell
pytest tests/playback 
```

Upper case TOML config files in `tests/config` are for the test suite. 
The other ones are for the walkthrough example of Chromosome 20 fishing.


## test script 

```
mamba create -n btest python=3.10 pip bioconda::gfatools bioconda::minimap2 bioconda::miniasm && mamba activate btest   
pip install boss_runs
pip install pytest pytest-cov pytest-timeout
git clone --recurse-submodules https://github.com/goldman-gp-ebi/BOSS-RUNS.git
pytest tests/base/ --cov --cov-report html
```

## Endpoints

- boss takes the output path of the sequencing protocol 
- boss takes the channels toml spec from readfish after getting the output path
- boss scans the output dir to get fastq files
- rf takes new masks and contigs from boss


## Test packaged version

`pip install dist/boss_runs-0.1.0-py3-none-any.whl --force-reinstall --no-deps`


## notes on the readfish wrapper 

- readfish_boss is based on the `targets.py` entry-point of readfish
- instead of launching it as an entry-point, we run _cli_base.main()
- the only mod in that function is to return parser and args instead of sending it
- targets.py is almost the same as the original, bar a few changed lines
- also overwrite the make_decision function -> BossBits.make_decision_boss()
- that function originally lives in _config.py

files/funcs to check after updates

targets.py
_cli_base.main()
_config.py.make_decision()

    

## Example simulation and analysis

https://github.com/W-L/boss_simulation_example


## System-level libs

mappy requires python3-dev
`sudo apt-get install python3-dev`


Only relevant if running on a machine without MinKNOW installed:

To use the python library `minknow_api` certain certificates from ONT need to be present on the System.
These certificates get placed by installing `sudo apt install minknow-core-standalone`, which requires adding the ONT software repos.
(https://nanoporetech.com/document/experiment-companion-minknow#installing-minknow-on-linu)
The certificates are then located at `/var/lib/minknow/data/rpc-certs/minknow/ca.crt`



