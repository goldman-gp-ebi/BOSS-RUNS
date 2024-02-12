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

The superclass Boss initialises and then runs `process_batch()` for each batch of new data. (Or `process_batch_sim()` for simulations)

The `process_batch()` functions take a method of one of the subclasses and execute that as their main body




