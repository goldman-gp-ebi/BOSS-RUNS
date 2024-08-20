from shutil import which
from sys import executable
from pathlib import Path
import sys
sys.path.insert(0,"/")
sys.path.insert(0,"../")

from boss.sampler import FastqStream_mmap, PafStream


"""
Snakemake file to prepare data for BOSS simulations. Produces:
- offsets of full-length reads (This is the only file needed for AEONS simulations, the rest is only needed for RUNS)
- truncated reads (mu-sized)
- mapping of full-length reads
- mapping of truncated reads 
- offsets of truncated mappings
- offsets of full-length mappings

- TODO add source tagging  

input reads must end with .fq for now -> TODO improve this
call with:
snakemake -s ../BOSS-RUNS/scripts/prepare_simulation_data.smk --config ref=ref.fa reads=reads.fq -n 
"""


fastq = config['reads']
ref = config['ref']

mu = 400
thread_no = 8
memory = 4000

try:
    mm2 = which("minimap2", path='/'.join(executable.split('/')[0:-1])).strip()
except:
    raise ValueError("minimap2 not found in PATH")

r = str(Path(fastq).stem)



rule all:
    input:
        # truncate fq
        f"{r}_trunc.fq",
        # offsets fq
        f"{r}.fq.offsets.npy",
        # full mappings
        f"{r}.paf",
        # offsets full mapping
        f"{r}.paf.offsets",
        # trunc mappings
        f"{r}_trunc.paf",
        # offsets trunc mappings
        f"{r}_trunc.paf.offsets"


rule truncate_fq:
    input: "{reads}.fq"
    output: "{reads}_trunc.fq"
    params:
        mu = mu
    shell:
        "cut -c -{params.mu} {input} > {output} "


rule scan_offsets_fq:
    input: "{reads}.fq"
    output: "{reads}.fq.offsets.npy"
    run:
        FastqStream_mmap(source=str(input))


rule map_paf:
    input: "{reads}.fq"
    output: "{reads}.paf"
    resources: mem_mb=memory
    threads: thread_no
    params:
        mm2 = mm2,
        ref = ref,
        threads = thread_no
    shell:
        "{params.mm2} -x map-ont -t {params.threads} --secondary=no -c {params.ref} {input} >{output}"


rule scan_offsets_paf:
    input:
        full = "{reads}.paf",
        truncated = "{reads}_trunc.paf"
    output:
        full_offsets = "{reads}.paf.offsets",
        trunc_offsets= "{reads}_trunc.paf.offsets"
    resources: mem_mb = 16000
    run:
        PafStream(paf_full=input.full, paf_trunc=input.truncated)






