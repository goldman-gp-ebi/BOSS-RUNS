# <img src=doc/BR_logo.png width="350" />


*B*enefit-*O*ptimising *S*hort-term *S*trategies for *R*ead*U*ntil *N*anopore *S*equencing


BOSS* strategies allow for dynamic, adaptive sampling during nanopore sequencing. New data is periodically ingested to 
generate updated decision strategies in order to maximise the information gain during the sequencing experiment.
[Readfish](https://github.com/LooseLab/readfish) is launched to run alongside and communicates the rejection signals 
to the sequencing machine. 


The method is described in this [article](https://doi.org/10.1038/s41587-022-01580-z)




<p align="center">
	<img src=doc/BR_flowchart.png width="500" />
</p>



## Changelog

- 2025/02/24 0.3.1 improved pytests, small bug fixes, updated documentation
- 2024/09/23 0.3.0 compatibility with minknow 6, improved performance of paring coverage counts, improved pytests and test coverage
- 2024/05/21 0.2.0 Updated internal readfish functionality to ensure compatibility with newest dorado versions



## Requirements

- Linux
- MinKNOW >=6.0 with dorado ~=7.4




## Installation


Recommended way of installing our software is in a conda/mamba environment. 
These commands will create an environment and install BOSS* alongside the necessary dependencies. 


```shell
mamba create -n boss python=3.10 pip bioconda::gfatools bioconda::minimap2 bioconda::miniasm && mamba activate boss   
pip install boss_runs
```



## Usage

BOSS* can be used in different ways depending on the aim of the sequencing experiment.
When providing a reference (and an optional index) BOSS-RUNS is executed.

The experiment is configured using 2 toml files, one for this software and one for readfish.
The toml file and defaults for BOSS* are described here:

```toml
[general]
name = "boss"                   # experiment name
wait = 60                       # waiting time between periodic updates
ref = ""                        # reference fasta file. Not specifying a file switches to BOSS-AEONS
mmi = ""                        # index of reference (will be built if ref is given but not mmi)

[live]
device = "X1"                   # position on sequencer
host = "localhost"              # host of sequencing device
port = 9502                     # port of sequencing device
data_wait = 100                 # wait for X Mb of data before first update
prom = false                    # switch for using a PromethION flowcell (experimental)

[optional]
reject_refs = ""                # comma-separated list of headers in reference from which to always reject
ploidy = 1                      # 1 or 2
```


By default, all whole genome(s) included in the input fasta file are considered of interest at the beginning of the experiment. 
It is possible to reject all reads from specific sequences in a fasta file. 
For this, provide fasta headers of the reference file, e.g.: `--reject_refs 1,2,3,X,Y,MT`



### Configuring readfish

A separate toml file needs to be given to BOSS that contains the configuration of readfish (according to [their instructions](https://github.com/LooseLab/readfish/blob/main/docs/toml.md))

Here's an example with two regions on a flowcell, where one is running this method and the second half acts as a control:

```toml
[caller_settings.dorado]  
config = 'dna_r10.4.1_e8.2_400bps_5khz_fast'
address = 'ipc:///tmp/.guppy/5555'
debug_log = 'live_reads.fq'

[mapper_settings.mappy_rs]
fn_idx_in = "../data/test.fasta"
debug_log = 'live_alignments.paf'
n_threads = 4

[[regions]]
name = "runs"
min_chunks = 0
max_chunks = 2
targets = []
single_on = "stop_receiving"
multi_on = "stop_receiving"
single_off = "unblock"
multi_off = "unblock"
no_seq = "unblock"
no_map = "unblock"
above_max_chunks = "stop_receiving"
below_min_chunks = "proceed"

[[regions]]
name = "control"
control = true
min_chunks = 0
max_chunks = 2
targets = []
single_on = "stop_receiving"
multi_on = "stop_receiving"
single_off = "stop_receiving"
multi_off = "stop_receiving"
no_seq = "stop_receiving"
no_map = "stop_receiving"
above_max_chunks = "stop_receiving"
below_min_chunks = "stop_receiving"
```

readfish is launched from within BOSS* for ease of use



### Starting BOSS*

After sequencing has started launch BOSS* with:

```shell
boss --toml path/to/toml --toml_readfish path/to/readfish/toml
```


BOSS* will initialise and start to periodically generate new decision strategies from the sequencing reads deposited by the sequencer.
If readfish is configured properly, the strategies will be reloaded automatically.
This triggers a message in readfish's logfile similar to: `Reloaded strategies for X sequences`.

When enough data is collected, BOSS* can be stopped by a keyboard interrupt (Ctrl+C).



## Testing


It is highly recommended to follow the [walkthrough provided in the readfish repository](https://github.com/LooseLab/readfish/tree/main?tab=readme-ov-file#testing)
to set up a playback experiment and test the functionality and interplay of the software and sequencing machine.

As soon as playback is running, BOSS* can be executed using toml files located in `tests/config`:


```shell
boss --toml tests/config/boss_ch20.toml --toml_readfish tests/config/boss_ch20_readfish.toml
``` 
              
This configures targeting of chromosome 20 with continuously updated decision strategies.
Let the playback sequencing run for a few minutes, then verify that the setup works:

1) `readfish` is rejecting reads from all chromosomes, except for #20. For this, look at the observed read lengths:

`readfish summary tests/config/boss_ch20_readfish.toml /path/to/sequencing/output/fastq_pass/`


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

for this, grep the log-file of `readfish` for all reloading events of updated strategies.

```shell 
grep "Reloaded" readfish.log
#2024-03-06 17:06:39,179 readfish.targets Reloaded strategies for 24 sequences
#2024-03-06 17:07:07,994 readfish.targets Reloaded strategies for 24 sequences
#2024-03-06 17:07:37,818 readfish.targets Reloaded strategies for 24 sequences
#2024-03-06 17:08:38,902 readfish.targets Reloaded strategies for 24 sequences
#2024-03-06 17:09:07,365 readfish.targets Reloaded strategies for 24 sequences
# ...
```


## Issues, questions, suggestions ...

Please use the issue tracker in this repository to get in touch!
Notes for development and code organisation can be found in `doc/developer_noted.md` 


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

