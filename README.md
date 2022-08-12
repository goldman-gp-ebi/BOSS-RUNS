# BOSS-RUNS

BOSS-RUNS allows for dynamic, adaptive sampling during nanopore sequencing. It periodically ingests all newly observed sequencing reads to generate updated decision strategies in order to maximise the information gain during the sequencing experiment. This method is implemented to run alongside readfish, which communicates the rejection signals to the sequencing machine and reloads BOSS-RUNS' decision strategies whenever they are updated. 

The method is described in this [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.07.938670v3)


## Requirements

- Linux
- MinKNOW >=5.0 with guppy >=6.0
- working setup of [readfish for BOSS-RUNS](https://github.com/LooseLab/readfish/tree/BossRuns/V0.0.1). 



## Installation

In parallel to BOSS-RUNS, a modified version of readfish needs to be set up, which is available at: [readfish for BOSS-RUNS](https://github.com/LooseLab/readfish/tree/BossRuns/V0.0.1).

Please follow the instructions of the readfish repository for installation instructions.

Then clone this repository:

`git clone https://github.com/goldman-gp-ebi/BOSS-RUNS.git`

BOSS-RUNS can be run in the same virtual environment used for readfish. You simply need to install a few additional dependencies:

`pip install numba natsort scipy toml`



Alternatively, BOSS-RUNS can be run in a separate conda environment. The required dependencies can be installed using the provided `yaml` file

`conda env create -f environment.yml`

To activate the environment:

`conda activate bossruns`


## Usage

BOSS-RUNS can be used in different ways depending on the aim of the sequencing experiment and the initial targets.

Reference fasta sequences are required for all targets. They will be indexed upon launch, or an index can optionally be provided. For creation, you can use:

`scripts/mappy_index_fasta.py references.fa`


By default, all whole genome(s) included in the input fasta file are considered of interest. 

Alternatively, a VCF that matches the input fasta file and contains regions/sites of interest (ROIs) can be provided. Only sites included in the VCF will be considered on-target.


### Starting readfish

readfish needs a single modification to utilize BOSS-RUNS' dynamically updated decision masks:
In the TOML configuration for readfish, add a single line pointing to the location where BOSS-RUNS will deposit new strategies, e.g. `mask = bossruns_name/masks`.
This path will follow the pattern `bossruns_{condition_name}/masks`, where `condition_name` is the name of the condition in the readfish TOML intended to use BOSS-RUNS.
A configuration for such a condition would like this for example:

```
[conditions.0]
name = "select_c20"                     <- condition name 
control = false
min_chunks = 0
max_chunks = 12
targets = ["20"]
single_on = "stop_receiving"
multi_on = "stop_receiving"
single_off = "unblock"
multi_off = "unblock"
no_seq = "proceed"
no_map = "proceed"
mask = "bossruns_select_c20/masks"      <- path to dynamic strategies (bossruns_{condition_name}/masks)
```

### Starting BOSS-RUNS

After sequencing has started and readfish is operating, the minimal command to launch BOSS-RUNS is:

`./bossruns.py --ref REF --device DEVICE --run_name CONDITION_NAME`

where `DEVICE` needs to be the name of the 'position' on the sequencer (displayed in MinKNOW overview),
and `CONDITION_NAME` is the same as in the readfish TOML.


BOSS-RUNS will initialise and start to periodically generate new decision strategies from the sequencing reads deposited by the sequencer.
If readfish is configured properly, the strategies will be reloaded automatically.
This triggers a message in readfish's logfile similar to: `Reloaded mask dict for FASTA_HEADER`.

After sequencing, BOSS-RUNS needs to be stopped by a keyboard interrupt (Ctrl+C).


### Arguments

Arguments can either be specified on the command line, by providing a parameter file, or a mixture or both (command-line arguments take priority).

```
usage: bossruns.py [-h] --ref REF [--ref_idx REF_IDX] [--vcf VCF] [--run_name RUN_NAME] [--ploidy PLOIDY] [--conditions] --device DEVICE
                   [--host HOST] [--port PORT] [--wait WAIT]

optional arguments:
  -h, --help           show this help message and exit
  --ref REF            Path to reference
  --ref_idx REF_IDX    Optional minimap index of reference
  --vcf VCF            Path to vcf file for ROIs
  --run_name RUN_NAME  Experiment identifier. If multiple conditions: must match name in channels.toml
  --ploidy PLOIDY      1 == haploid, 2 == diploid
  --conditions         Multiple conditions on a single flowcell, used to assign channels
  --device DEVICE      Name of device/sequencing position in MinKNOW
  --host HOST          hostname of sequencing device
  --port PORT          port of sequencing device
  --wait WAIT          Period between strategy updates (sec.)


```

A parameter file can be supplied using `@` on the command line. E.g.: `python bossruns.py @params.txt` with the structure (one space-delimited argument per line):

```
--argumentX valueX
--argumentY valueY
```

## Walkthrough


TODO

Sample experiment using Playback mode (under construction ..)


## Issues, questions, suggestions ...

Please use the issue tracker in this repository to get in touch!


## Citation

add bibtex entry 


## License

Licensed under GPLv3

