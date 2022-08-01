# BOSS-RUNS

BOSS-RUNS allows for dynamic, adaptive sampling during nanopore sequencing. It periodically ingests all newly observed sequencing reads to generate updated decision strategies in order to maximise the information gain during the sequencing experiment. This method is implemented to run alongside readfish, which communicates the rejection signals to the sequencing machine and reloads BOSS-RUNS' decision strategies whenever they are updated. 

The method is described in this [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.07.938670v3)


## Requirements

- Linux
- MinKNOW >=5.0 with guppy >=6.0
- working setup of [readfish for BOSS-RUNS](https://github.com/LooseLab/readfish/tree/BossRuns/V0.0.1). 



## Installation

Clone this repository:

`git clone https://github.com/goldman-gp-ebi/BOSS-RUNS.git`

Install a conda environment incl. required dependencies using the provided `yaml` file

`conda env create -f environment.yml`

Activate the environment:

`conda activate bossruns`


In parallel, a modified version of readfish needs to be set up, available at: [readfish for BOSS-RUNS](https://github.com/LooseLab/readfish/tree/BossRuns/V0.0.1).

Please follow the instructions of the readfish repository for installation.


## Usage

BOSS-RUNS can be used in different ways depending on the aim of the sequencing experiment and the initial targets.

Reference fasta sequences are required for all targets. They will be indexed upon launch, or an index can optionally be provided. For creation, you can use:

`scripts/mappy_index_fasta.py references.fa`


By default, all whole genome(s) included in the input fasta file are considered of interest. 

Alternatively, a VCF that matches the input fasta file and contains regions/sites of interest (ROIs) can be provided. Only sites included in the VCF will be considered on-target.


After sequencing has started and readfish is operating, the minimal command to launch BOSS-RUNS is:

`./bossruns.py --ref REF --device DEVICE`

where `DEVICE` needs to be the name of the 'position' on the sequencer (displayed in MinKNOW overview).

BOSS-RUNS will initialise and start to periodically generate new decision strategies from the sequencing reads deposited by the sequencer. If readfish is configured properly, the strategies will be reloaded automatically (triggering a message in readfish's logfile). BOSS-RUNS needs to be stopped by a keyboard interrupt (Ctrl+C) after the sequencing is completed.


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

Sample experiment using Playback mode (under construction ..)


## Issues, questions, suggestions ...

Please use the issue tracker in this repository to get in touch!


## Citation

add bibtex entry 


## License

Licensed under GPLv3

