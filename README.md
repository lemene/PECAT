# Introduction

PECAT is a phased error correction and assembly tool for long reads. It includes a haplotype-aware correction method and an efficient diploid assembly method. 

# Dependency

+ minimap2 
+ racon
+ perl (v5.22.1+)


# Installation

### Build from source codes

```shell
$ git clone https://github.com/lemene/PECAT.git
$ cd PECAT
$ make
```

After installation, all the executable files can be found in `PECAT/build/bin`. 


# Quick Start

### Step 1: Create a config file

Create a config file template using the following command:

```shell
$ PECAT/build/bin/pecat.pl config cfgfile
```

The template looks like

``` shell
project=
reads=
genome_size=1
threads=4
memory=0
cleanup=0
compress=0
grid=auto:0
prep_min_length=3000
prep_output_coverage=80
corr_iterate_number=1
corr_block_size=4000000000
corr_correct_options=
corr_filter_options=--filter0=l=2000:al=2000:alr=0.5
corr_rd2rd_options=-x ava-pb
corr_output_coverage=80
align_block_size=4000000000
align_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500
align_filter_options=--filter0=l=2000:aal=4000:oh=3000:ohr=0.3 --task=extend --filter1=oh=100:ohr=0.01
asm1_assemble_options=
phase_method=
phase_rd2ctg_options=-x map-pb -c -p 0.5 -r 1000
phase_use_reads=1
phase_phase_options=
asm2_assemble_options=
polish_map_options=-x map-pb
polish_filter_options=--filter0 oh=1000:ohr=0.1
polish_cns_options=
```

Filling and modifying the relative information, we have

``` shell
project=ecoli
reads=ecoli-reads.fasta
genome_size=4600000
  ......
```

### Step 2: Correct raw reads
Correct the raw noisy reads using the following command:
``` Shell
$ PECAT/build/bin/pecat.pl correct cfgfile
```
The pipeline extracts `prep_output_coverage` longest raw reads for correction and extracts `corr_output_coverage` longest corrected reads for assembly, which are in the file `./ecoli/1-correct/corrected_reads.fasta`

### Step 3: Assemble contigs

After correcting the raw reads, we assemble the contigs using the following command. 

```Shell
$ PECAT/build/bin/pecat.pl assemble cfgfile
```
The assembled contigs are in the file `./ecoli/3-assemble/primary.fasta`.

### Step 4: unzip contigs

After assembling the contigs, we run the bridging-step using the following command.

```Shell
$ PECAT/build/bin/pecat.pl unzip cfgfile
```
It generates primary/alternate contigs in the files `./ecoli/5-assemble/primary.fasta` and `./ecoli/5-assemble/alternate.fasta`. In this step, the pipeline invokes `minimap2` and `racon` to polishes the contigs. The polished contigs are in the files `./ecoli/6-polishe/primary.fasta` and `./ecoli/6-polishe/alternate.fasta`

All commands checks and runs the preceding steps first. So `PECAT/build/bin/pecat.pl unzip cfgfile` will run step 1-4 in turn.

# Running with multiple computation nodes

The pipeline script is written with [plgd](https://github.com/lemene/plgd). It supports PBS, SGE, LSF and Slurm systems. The follow parameter in the config file need to be set:
```shell
grid= auto:4
```

In the above example, `auto` means the pipeline automatically detects the type of cluster system. `pbs`, `sge`, `lsf` and `slurm` represent the corresponding systems, respectively. the pipeline. `4` computation nodes will be used and each computation node will run with `threads` CPU threads.



