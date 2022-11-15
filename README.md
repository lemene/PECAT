# Introduction

PECAT is a phased error correction and assembly tool for long reads. It includes a haplotype-aware correction method and an efficient diploid assembly method. 

# Dependency

+ minimap2 
+ racon
+ perl (v5.22.1+)


# Installation

### Build from source codes

```shell
$ git clone --recursive https://github.com/lemene/PECAT.git
$ cd PECAT
$ make
```
or
```shell
$ git clone  https://github.com/lemene/PECAT.git
$ cd PECAT
$ git submodule init
$ git submodule update
$ make
```

After installation, all the executable files can be found in `PECAT/build/bin`. 

The executale files can be packaged with `bash scripts\package.sh`. `pecat_xxx_yyy_.tar.gz` is created in `PECAT/build`.

# Quick Start
*Before running PECAT please do not forget to add the paths of minimap2, racon and perl to the system PATH.*

### Step 1: Create a config file

Create a config file template using the following command:

```shell
$ PECAT/build/bin/pecat.pl config cfgfile
```

The template looks like

``` shell
project=
reads=
genome_size=
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

#### Config for PacBio CLR diploid datasets
```Shell
project=
reads= 
genome_size=
threads=48
cleanup=0
compress=0
grid=local:1
prep_min_length=3000
prep_output_coverage=80
corr_iterate_number=1
corr_block_size=4000000000
corr_filter_options=--filter0=:al=5000:alr=0.5:aalr=0.5:oh=1000:ohr=0.1
corr_correct_options=--score=weight --aligner diff --min_coverage 4
corr_rd2rd_options=-x ava-pb
corr_output_coverage=80
align_block_size=4000000000
align_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500
align_filter_options=--filter0=l=5000:al=5000:aal=6000:aalr=0.5:oh=1000:ohr=0.1 --task=extend --filter1=oh=100:ohr=0.01
asm1_assemble_options=
phase_method=0
phase_rd2ctg_options=-x map-pb -c -p 0.5 -r 1000
phase_use_reads=1
phase_phase_options=
phase_filter_options=
asm2_assemble_options=
polish_use_reads=1
polish_map_options = -x asm20
polish_filter_options= --filter0 oh=1000:ohr=0.1
polish_cns_options =
```

See [config.md](doc/config.md).

#### Config for Nanopore diploid datasets
```shell
project=
reads= 
genome_size= 
threads=48
cleanup=0
compress=0
grid=local:1
prep_min_length=3000
prep_output_coverage=80
corr_iterate_number=2
corr_block_size=4000000000
corr_filter_options=--filter0=l=5000:al=2500:alr=0.5:aal=5000:oh=3000:ohr=0.3
corr0_correct_options=--score=weight:lc=16 --aligner edlib  --min_coverage 2
corr_correct_options=--score=weight:lc=10 --aligner edlib   --min_coverage 4
corr_rd2rd_options=-x ava-ont -f 0.001
corr1_rd2rd_options=-x ava-ont -w30 -k19 -f 0.001
corr_output_coverage=80
align_block_size=4000000000
align_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500
align_filter_options=--filter0=l=5000:aal=6000:aalr=0.5:oh=3000:ohr=0.3 --task=extend --filter1=oh=300:ohr=0.03
asm1_assemble_options=
phase_method=0
phase_rd2ctg_options=-x map-ont -c -p 0.5 -r 1000
phase_use_reads=1
phase_phase_options= --coverage lc=30 --phase_options icr=0.2:icc=6:sc=8
phase_filter_options= --threshold=1000
asm2_assemble_options=
polish_use_reads=0
polish_filter_options= --filter0 oh=1000:ohr=0.1
polish_map_options = -x map-ont -k19 -w10
polish_cns_options =
```
**Note:** For large genomes such as cattle, we strongly suggest adding the parameter `-f 0.005` or `-f 0.002` to `corr_rd2rd_options`, `corr0_rd2rd_options`, `corr1_rd2rd_options` and `align_rd2rd_options`. The parameter is passed to `minimap2`, which means to filter out top 0.005 or 0.002 fraction of repetitive minimizers. It outputs less candidate overlaps, which reduces disk usage and speeds up error correction step and assembling step.

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

### Step 4: Unzip contigs

After assembling the contigs, we run the bridging-step using the following command.

```Shell
$ PECAT/build/bin/pecat.pl unzip cfgfile
```
It generates primary/alternate contigs in the files `./ecoli/5-assemble/primary.fasta` and `./ecoli/5-assemble/alternate.fasta`. In this step, the pipeline invokes `minimap2` and `racon` to polishes the contigs. The polished contigs are in the files `./ecoli/6-polishe/primary.fasta` and `./ecoli/6-polish/alternate.fasta`

All commands checks and runs the preceding steps first. So `PECAT/build/bin/pecat.pl unzip cfgfile` will run step 1-4 in turn.

# Running with multiple computation nodes

The pipeline script is written with [plgd](https://github.com/lemene/plgd). It supports PBS, SGE, LSF and Slurm systems. The follow parameter in the config file need to be set:
```shell
grid= auto:4
```
In the above example, `auto` means the pipeline automatically detects the type of cluster system. `pbs`, `sge`, `lsf` and `slurm` represent the corresponding systems, respectively. the pipeline. `4` computation nodes will be used and each computation node will run with `threads` CPU threads.


# Contact
+ Nie Fan, niefan@csu.edu.cn

