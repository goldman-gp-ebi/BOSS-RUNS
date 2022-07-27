# BOSS-RUNS

BOSS-RUNS allows for dynamic, adaptive sampling during nanopore sequencing. It periodically ingests all newly observed sequencing reads to generate updated decision strategies in order to maximise the information gain during the sequencing experiment. This method is implemented to run alongside Readfish, which communicates the rejection signals to the sequencing machine and reloads BOSS-RUNS' decision strategies whenever they are updated. 

The method is described in this [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.07.938670v3)

## Installation

For now, the easiest way is to clone this repository

`git clone https://github.com/goldman-gp-ebi/BOSS-RUNS.git`

and install a conda environment with the required dependencies using the provided .yaml file

`conda env create -f environment.yml`

the environment is then activated with

`conda activate bossruns`


In parallel, a modified version of Readfish needs to be set up. This version is available at: [Readfish for BOSS-RUNS](https://github.com/LooseLab/readfish/tree/BossRuns/V0.0.1).

Please follow the instructions of the Readfish repository to install this branch (specific detailed instructions to follow soon).


## Usage

Using BOSS-RUNS requires reference sequences (e.g. `references.fa`) to map reads against. These need to be indexed with minimap2 first, for which we provide a script:

`scripts/mappy_index_fasta.py references.fa`

Secondly, a mask file that indicates the initial sites or regions of interest (ROIs) must be provided. This can be generated for different scenarios.

- If one or multiple whole genomes (e.g. collected in `references.fa`) are of interest, the mask file can be generated using:

	`scripts/ROI_wgs.py --fa references.fa --out out_name`

	this will create 2 files: `out_name.mask.npz` and `out_name.ref.fa`, which are required for BOSS-RUNS


- If the targets are ROIs (regions or individual sites), the mask file can instead be generated from a VCF:

	`scripts/ROI_vcf.py --vcf VCF --out out_name --species input_species`

	in this case only `out_name.mask.npz` is required for BOSS-RUNS


A working setup of adaptive sampling using the BOSS-RUNS branch of Readfish is needed to use BOSS-RUNS. Once Readfish is operating, BOSS-RUNS can be launched by:

`python bossruns.py --roi_mask ROI_MASK --ref_idx REF_IDX --run_name RUN_NAME --device DEVICE`

where `RUN_NAME` can be any identifier and `DEVICE` needs to be the name of the 'position' on the sequencer (displayed in MinKNOW overview).

BOSS-RUNS will initialise and start to periodically generate new decision strategies from the sequencing reads deposited by the sequencer. If Readfish is configured properly, the strategies will be reloaded automatically (triggering a message in Readfish's logfile). After the experiment has finished BOSS-RUNS needs to be stopped by a keyboard interrupt (Ctrl+C).


### Arguments

Command-line arguments can be specified as usual, by providing a parameter file, or a mixture or both (command-line arguments take priority).

```
usage: bossruns.py [-h] --roi_mask ROI_MASK --ref_idx REF_IDX --run_name RUN_NAME [--ref REF] --device DEVICE [--host HOST] [--port PORT] [--conditions] [--whole_genome] [--wait WAIT] [--ploidy PLOIDY]

optional arguments:
  -h, --help           show this help message and exit
  --roi_mask ROI_MASK  Path to ROI mask
  --ref_idx REF_IDX    Index of reference to map against
  --run_name RUN_NAME  Any identifier. If multiple conditions: must match name in channels.toml
  --ref REF            Path to reference fasta, used for priors. If not provided, priors will be uniform
  --device DEVICE      Name of device/sequencing position in MinKNOW
  --host HOST          hostname of sequencing device
  --port PORT          port of sequencing device
  --conditions         Multiple conditions on a single flowcell, used to assign channels
  --whole_genome       switch for several algorithms, if whole genome(s) of interest
  --wait WAIT          Period between strategy updates (sec.)
  --ploidy PLOIDY      1 == haploid, 2 == diploid

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


## License

Licensed under GPLv3

