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

readfish can then be started with the subcommand `boss-runs` (instead of `targets`). E.g. like this:
```
readfish boss-runs --device DEVICE 
                   --experiment-name NAME
                   --toml TOML
                   --log-file LOGFILE
```


### Starting BOSS-RUNS

After sequencing has started and readfish is operating, the minimal command to launch BOSS-RUNS is:

```
./bossruns.py --ref REF
               --device DEVICE
               --run_name CONDITION_NAME
```

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



This walkthrough is adapted from the [readfish repository](https://github.com/LooseLab/readfish)

### Setting up a playback run

To test your configuration of readfish and BOSS-RUNS we recommend first running a playback experiment.
Download an open access [bulk FAST5 file](http://s3.amazonaws.com/nanopore-human-wgs/bulkfile/PLSP57501_20170308_FNFAF14035_MN16458_sequencing_run_NOTT_Hum_wh1rs2_60428.fast5). (Attention 21Gb!)

To configure a run for playback, you need to edit a sequencing TOML file located at `/opt/ont/minknow/conf/package/sequencing`.
Edit `sequencing_MIN106_DNA.toml` and under the entry `[custom_settings]` add:

`simulation = "/full/path/to/your_bulk.FAST5"`

and set the parameter `break_reads_after_seconds = 1.0` to `break_reads_after_seconds = 0.4`

In the MinKNOW GUI select Reload Scripts (vertical dot menu at `Select positions` screen when starting a sequencing run).

Insert a configuration test flowcell into the sequencing device and start a sequencing run (selecting the corresponding flow cell type to the edited script, i.e. `FLO-MIN106`).

The run should start and immediately begin a mux scan. Let it run for a few minutes.

### Starting readfish 


This requires access to a guppy basecall server and a TOML file. Here's an example TOML file:

```
[caller_settings]
config_name = "dna_r9.4.1_450bps_hac"
host = "127.0.0.1"
port = 5555

[conditions]
reference = "/path/to/reference.mmi"

[conditions.0]
name = "select_c20"
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
mask = "bossruns_select_c20/masks"
```

This configures `readfish` to target all reads from chromosome 20 and to continuously read the dynamically updated decision strategies from `BOSS-RUNS`.

You simply need to modify the reference field to a minimap2 index of the human genome.

Modify the targets fields for each condition to reflect the naming convention used in your index. 
This is the sequence name only, up to but not including any whitespace. e.g. `>20 human chromosome 20` would become `20`.

`readfish` can then be launched with 

```
readfish boss-runs --device DEVICE
                   --experiment-name "EXPERIMENT_NAME"
                   --toml example.toml
                   --log-file readfish.log
```


### Starting BOSS-RUNS


After `readfish` is running, you can launch `BOSS-RUNS`. In this walkthrough we are using a reference that contains only chromosome 20. 
This is purely for testing purposes and does not make for an interesting biological question. 

First download and unzip the reference for chromosome 20:

```
mkdir -p ref && wget -O ref/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz && gunzip ref/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
```

Then launch BOSS-RUNS with the following parameters:

```
--run_name select_c20
--ref /ref/Homo_sapiens.GRCh38.dna.chromosome.20.fa
--device DEVICE
--testing
```

`--run_name` needs to match the condition name in the readfish TOML file.


Let the playback sequencing run for a few minutes. 


### Check that it's working


There are 2 things to verify that the setup works:

- (i) `readfish` is rejecting reads from all chromosomes, except for #20. For this, we look at the observed read lengths:

`readfish summary example.toml /path/to/sequencing/output/fastq_pass/`


Check that the mean read length for the enriched chromosome is larger than for the remaining chromosomes.   


TODO replace readfish summary output

```
    contig  number      sum   min     max    std   mean  median     N50
      chr1    1326  4187614   142  224402  14007   3158     795   48026
     chr10     804  2843010   275  248168  15930   3536     842   47764
     chr11     672  2510741   184  310591  18572   3736     841   73473
     chr12     871  2317742   292  116848   9929   2661     825   37159
     chr13     391  1090012   227  189103  12690   2788     781   41292
     chr14     469  2323329   275  251029  20107   4954     830   68887
     chr15     753  2189326   180  154830  12371   2907     812   40686
     chr16     522  1673329   218  166941  12741   3206     862   39258
     chr17     484  1609208   191  169651  15777   3325     816   73019
     chr18     483  1525953   230  252901  14414   3159     813   40090
     chr19     664  1898289   249  171742  13181   2859     820   46271
      chr2    1474  4279420   234  222310  13090   2903     820   43618
     chr20     489  1622910   229  171322  13223   3319     887   33669
     chr21      32  1221224  1053  223477  56923  38163   13238  112200
     chr22      47   724863   244  184049  28113  15423    6781   33464
      chr3    1142  3554814   243  247771  15173   3113     760   62683
      chr4    1224  4402210   210  221084  15769   3597     820   66686
      chr5    1371  4495150   205  330821  16699   3279     801   65394
      chr6     978  2725891   246  146169  10995   2787     791   37791
      chr7    1039  3027136   166  263043  14705   2914     798   56567
      chr8     848  2581406   238  229150  15618   3044     772   44498
      chr9     893  3028224   259  247975  16011   3391     802   54953
      chrM     144   216047   215   20731   2562   1500     864    1391
      chrX     868  3124552   238  192451  15594   3600     832   49047
      chrY       8    47071   510   31654  10743   5884    1382   31654    
```


- (ii) `readfish` is using dynamically updated decision strategies

for this, we can simply grep the log-file of `readfish` for all reloading events of updated strategies.

```
grep "Reloaded" readfish.log
```

This should produce something like this, with updates after 90 seconds each (by default):


```
TODO add grep output from readfish log
```



### Deactivating playback behaviour

After testing, remove the `simulation = ` line from the `sequencing_MIN106_DNA.toml` file and reload scripts in MinKNOW GUI (as above).



## Issues, questions, suggestions ...

Please use the issue tracker in this repository to get in touch!


## Citation

add bibtex entry 


## License

Licensed under GPLv3

