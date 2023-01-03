# <img src=doc/BR_logo.png width="500" />


*B*enefit-*O*ptimising *S*hort-term *S*trategies for *R*ead*U*ntil *N*anopore *S*equencing


BOSS-RUNS allows for dynamic, adaptive sampling during nanopore sequencing. It periodically ingests all newly observed sequencing reads to generate updated decision strategies in order to maximise the information gain during the sequencing experiment. This method is implemented to run alongside readfish, which communicates the rejection signals to the sequencing machine and reloads BOSS-RUNS' decision strategies whenever they are updated. 

The method is described in this [open access paper](https://doi.org/10.1038/s41587-022-01580-z) and this [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.02.07.938670v3)



<p align="center">
	<img src=doc/BR_flowchart.png width="700" />
</p>




## Requirements

- Linux
- MinKNOW >=5.0 with guppy >=6.0
- working setup of [readfish for BOSS-RUNS](https://github.com/LooseLab/readfish/tree/BossRuns/V0.0.2). 



## Installation

In parallel to BOSS-RUNS, a modified version of readfish needs to be set up,
which is available at: [readfish for BOSS-RUNS](https://github.com/LooseLab/readfish/tree/BossRuns/V0.0.2).

To install both readfish and BOSS-RUNS:

- create a virtual environment
    
`python3 -m venv bossruns`

`. bossruns/bin/activate`

- install dependencies and readfish

`pip install --upgrade pip`

`pip install git+https://github.com/nanoporetech/read_until_api@v3.0.0`

`pip install git+https://github.com/LooseLab/readfish@BossRuns/V0.0.2`

- `ont_guppy_client_lib` needs a version specific to your guppy version. On OSes using `apt` you can find your guppy version using:

`apt list --installed ont-guppy* | tail -n 1 | cut -f2 -d' ' | cut -f1 -d'-' >guppy_version`

`cat guppy_version`

`pip install ont_pyguppy_client_lib==$(cat guppy_version)`

- install a few more dependencies for BR and clone the repository

`pip3 install natsort scipy numba toml`

`git clone https://github.com/goldman-gp-ebi/BOSS-RUNS.git`




Alternatively, if you already have the correct version of readfish set up, BOSS-RUNS can be run in a separate conda environment.
The required dependencies can then be installed using the provided `yaml` file

`conda env create -f environment.yml`

To activate the environment:

`conda activate bossruns`


## Usage

BOSS-RUNS can be used in different ways depending on the aim of the sequencing experiment and the initial targets.

Reference fasta sequences are required for all targets. They will be indexed upon launch, or an index can optionally be provided. For creation, you can use:

`BOSS_RUNS/scripts/mappy_index_fasta.py references.fa`


- By default, all whole genome(s) included in the input fasta file are considered of interest. 

- You can also choose to reject all reads from specific sequences in your fasta file. For this, provide fasta headers of the reference file, e.g.: `--reject_refs 1,2,3,X,Y,MT`

- Alternatively, a VCF that matches the input fasta file and contains regions/sites of interest (ROIs) can be provided. Only sites included in the VCF will be considered on-target.


### Starting readfish

readfish needs a single modification to utilize BOSS-RUNS' dynamically updated decision masks:
In the TOML configuration for readfish, add a single line pointing to the location where BOSS-RUNS will deposit new strategies, e.g. `mask = bossruns_name/masks`.
This path will follow the pattern `bossruns_{condition_name}/masks`, where `condition_name` is the name of the condition in the readfish TOML intended to use BOSS-RUNS.
A configuration for such a condition might look like this:

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
readfish boss-runs --device DEVICE \
                   --experiment-name NAME \
                   --toml TOML \
                   --log-file LOGFILE
```


### Starting BOSS-RUNS

After sequencing has started and readfish is operating, the minimal command to launch BOSS-RUNS is:

```
./bossruns.py --ref REF --device DEVICE --run_name CONDITION_NAME
```

where `DEVICE` needs to be the name of the 'position' on the sequencer (displayed in MinKNOW overview),
and `CONDITION_NAME` is the same as in the readfish TOML.


BOSS-RUNS will initialise and start to periodically generate new decision strategies from the sequencing reads deposited by the sequencer.
If readfish is configured properly, the strategies will be reloaded automatically.
This triggers a message in readfish's logfile similar to: `Reloaded mask dict for FASTA_HEADERS`.

After sequencing, BOSS-RUNS needs to be stopped by a keyboard interrupt (Ctrl+C).


### Arguments

Arguments can either be specified on the command line, by providing a parameter file, or a mixture or both (command-line arguments take priority).

`python BOSS-RUNS/bossruns.py --help`

```
usage: bossruns.py [-h] --ref REF [--ref_idx REF_IDX] [--run_name RUN_NAME] [--vcf VCF] [--reject_refs REJECT_REFS] [--ploidy PLOIDY] [--conditions] --device DEVICE [--host HOST] [--port PORT] [--wait WAIT]
                   [--ckp CKP]

optional arguments:
  -h, --help            show this help message and exit
  --ref REF             Path to reference
  --ref_idx REF_IDX     Optional minimap index of reference
  --run_name RUN_NAME   Experiment identifier. Must match name of [conditions.X] in readfish toml file
  --vcf VCF             Path to vcf file for ROIs
  --reject_refs REJECT_REFS
                        reject all reads of some entries in reference file, i.e. chromosomes or species. Multiple headers can be comma separated
  --ploidy PLOIDY       1 == haploid, 2 == diploid
  --conditions          Multiple conditions on a single flowcell, used to assign channels
  --device DEVICE       Name of device/sequencing position in MinKNOW
  --host HOST           hostname of sequencing device
  --port PORT           port of sequencing device
  --wait WAIT           Period between strategy updates (sec.)

```

A parameter file can be supplied using `@` on the command line. E.g.: `python bossruns.py @params.txt` with the structure (one space-delimited argument per line):

```
--argumentX valueX
--argumentY valueY
```

## Walkthrough for testing



This walkthrough is adapted from the [readfish repository](https://github.com/LooseLab/readfish)

### Setting up a playback run

To test your configuration of readfish and BOSS-RUNS we recommend first running a playback experiment.
Download an open access [bulk FAST5 file](http://s3.amazonaws.com/nanopore-human-wgs/bulkfile/PLSP57501_20170308_FNFAF14035_MN16458_sequencing_run_NOTT_Hum_wh1rs2_60428.fast5). (Attention 21Gb!)

To configure a run for playback, you need to edit a sequencing TOML file located at `/opt/ont/minknow/conf/package/sequencing`.

- Edit `sequencing_MIN106_DNA.toml` and under the entry `[custom_settings]` add:

`simulation = "/full/path/to/your_bulk.FAST5"`

- and set the parameter `break_reads_after_seconds = 1.0` to `break_reads_after_seconds = 0.4`

- In the MinKNOW GUI select Reload Scripts (vertical dot menu at `Select positions` screen when starting a sequencing run).

- Insert a configuration test flowcell into the sequencing device and start a sequencing run (selecting the corresponding flow cell type to the edited script, i.e. `FLO-MIN106`).

- The run should start and immediately begin a mux scan. Let it run for a few minutes.

### Starting readfish 


This requires access to a guppy basecall server and a TOML file. Here's an example TOML file for this walkthrough:

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

This configures `readfish` to target all reads from chromosome 20 and to continuously read the dynamically updated decision strategies from `BOSS-RUNS` (`mask = "bossruns_select_c20/masks"`).

You simply need to modify the reference field to a minimap2 index of the human genome.

Modify the `targets` field to reflect the naming convention used in your index. 
This is the sequence name only, up to but not including any whitespace. e.g. the fasta header `>20 human chromosome 20` would become `20`.

`readfish` can then be launched with 

```
readfish boss-runs --device DEVICE \
                   --experiment-name "EXPERIMENT_NAME" \
                   --toml example.toml \
                   --log-file readfish.log
```


### Starting BOSS-RUNS


After `readfish` is running, you can launch `BOSS-RUNS` using the same reference file and indicate which chromosomes are not of interest a priori.

```
./bossruns.py --run_name select_c20 \
              --ref /data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
              --ref_idx /data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.mmi \
              --reject_refs 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,X,Y,MT \
              --device DEVICE \
              --testing
```

`run_name` needs to match the condition name in the readfish TOML file.


Let the playback sequencing run for a few minutes. 


### Check that it's working


There are 2 things to verify that the setup works:

1) `readfish` is rejecting reads from all chromosomes, except for #20. For this, we look at the observed read lengths:

`readfish summary example.toml /path/to/sequencing/output/fastq_pass/`


Check that the mean read length for the enriched chromosome is larger than for the remaining chromosomes.   
(In this example the read lengths of depleted chromosomes are still rather long due to slow base calling) 



```
contig  number      sum  min     max    std   mean  median    N50
     1    1151  1577570  224  197871   7002   1371     581   7096
    10    1012  1147739  211  105585   5128   1134     578   1149
    11     836  1207274  225  212740   8385   1444     590   7727
    12     741   940284  192   81667   5419   1269     569   2143
    13     477   412434  182   54642   2694    865     556    898
    14     771   988124  176  114908   5616   1282     573   7301
    15     540   882612  223   99327   7459   1634     576   8258
    16     425   486124  209  107335   5579   1144     533   1341
    17     732  1099546  195  219260   9307   1502     614   8725
    18     180   310512  228   40407   4565   1725     638   8213
    19     592   824374  222  121745   6374   1393     634   5811
     2    1477  1572749  190   89358   3780   1065     561   1154
    20      50  1163599  227  145149  37887  23272    2214  86773   <---
    21     392   424448  199  118785   6053   1083     548   1213
    22     178   286534  188   49058   4788   1610     655   9884
     3    1198  1638753  201  172141   7473   1368     578   6548
     4    1370  1807013  174  160366   6708   1319     582   6917
     5    1408  2144345  162  212394   9884   1523     544   8397
     6     656  1013424  231  118194   5819   1545     599   7268
     7    1026   932717  185   66384   2972    909     563    914
     8     906  1133194  210  162732   5930   1251     564   2210
     9    1046  1533653  200  248867   9656   1466     552   8059
    MT      19   143721  591   16467   6140   7564    6809  13257
     X    1515  1490398  199  132776   4652    984     526   1014
     Y       9     6531  427    1895    475    726     517    628

```


2) `readfish` is using dynamically updated decision strategies

for this, we can simply grep the log-file of `readfish` for all reloading events of updated strategies.

```
grep "Reloaded" readfish.log
```

This should produce an output similar to this, with updates every ~90 seconds (by default):


```
2022-08-19 22:58:12,327 ru.ru_gen_boss_runs Reloaded mask dict for dict_keys(['2', '14', '21', '18', '8', '20', '16', '15', '6', '17', '19', '1', '22', '13', '5', '7', '9', '4', '12', 'Y', '10', '11', '3', 'X'])
2022-08-19 22:59:32,695 ru.ru_gen_boss_runs Reloaded mask dict for dict_keys(['2', '14', '21', '18', '8', '20', '16', '15', '6', '17', '19', '1', '22', '13', '5', '7', '9', '4', '12', 'Y', '10', '11', '3', 'X'])
2022-08-19 23:01:02,290 ru.ru_gen_boss_runs Reloaded mask dict for dict_keys(['2', '14', '21', '18', '8', '20', '16', '15', '6', '17', '19', '1', '22', '13', '5', '7', '9', '4', '12', 'Y', '10', '11', '3', 'X'])
2022-08-19 23:02:31,781 ru.ru_gen_boss_runs Reloaded mask dict for dict_keys(['2', '14', '21', '18', '8', '20', '16', '15', '6', '17', '19', '1', '22', '13', '5', '7', '9', '4', '12', 'Y', '10', '11', '3', 'X'])

```



### Deactivating playback behaviour

After testing, remove the `simulation =` line from the `sequencing_MIN106_DNA.toml` file and reload scripts in MinKNOW GUI (as above).



## Issues, questions, suggestions ...

Please use the issue tracker in this repository to get in touch!


## Citation


```
@article{weilgunyDynamicAdaptiveSampling2023,
  title = {Dynamic, adaptive sampling during nanopore sequencing using {{Bayesian}} experimental design},
  author = {Weilguny, Lukas and De Maio, Nicola and Munro, Rory and Manser, Charlotte and Birney, Ewan and Loose, Matthew and Goldman, Nick},
  year = {2023},
  journal = {Nature Biotechnology},
  publisher = {{Nature Publishing Group}},
  doi = {10.1038/s41587-022-01580-z}
}

```


## License

Licensed under GPLv3

