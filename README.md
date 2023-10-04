# Introduction

PECAT is a phased error correction and assembly tool for long reads. It includes a haplotype-aware correction method and an efficient diploid assembly method. 

# Dependency

+ [python3](https://www.python.org) (3.6+)
+ [minimap2](https://github.com/lh3/minimap2) (2.17+)
+ [racon](https://github.com/lbcb-sci/racon) (v1.4.21+)
+ [perl](https://github.com/Perl) (v5.22.1+)
+ [samtools](https://github.com/samtools/samtools)  (1.7+)
+ [clair3](https://github.com/HKU-BAL/Clair3)  (v0.1-r12+) (optional)
+ [medaka](https://github.com/nanoporetech/medaka)  (1.7.2+) (optional)

# Installation

## Installing PECAT from source codes

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

After building, all the executable files can be found in `PECAT/build/bin`. We can run `PECAT/build/bin/pecal.pl` or add the path to the system PATH and run `pecal.pl`.

### zlib not found
```Shell
wget -c http://www.zlib.net/zlib-1.2.13.tar.gz
tar -xzf zlib-1.2.13.tar.gz
cd zlib-1.2.13
./configure && make
cd ..
export C_INCLUDE_PATH=`pwd`/zlib-1.2.13:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=`pwd`/zlib-1.2.13:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=`pwd`/zlib-1.2.13:$LIBRARY_PATH
make
```


## Installing PECAT using conda
Use [Bioconda](https://bioconda.github.io/).

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Install PECAT.
```
conda create -n pecat-env
conda activate pecat-env
conda install pecat
```
Then we can run `pecal.pl`.

## Installing third-party tools using conda
PECAT depends on other tools, and their paths need to be added to the system PATH. We recommend using conda to install the third-party tools.
```
conda create -n pecat-env
conda activate pecat-env
conda install minimap2=2.24 racon=1.5 perl=5.32 samtools=1.17 python=3.11 # clair3 medaka
```

### Installing and configuring clair3 and medaka
When we installed clair3 and medaka using conda, we encountered a conflict between clair3(v0.1-r12) and medaka (1.7.2). Only one of them can be installed. If you also fail to install the tools, we recommend using singularity or docker to invoke them. 


#### Using singularity
Download the images
```Shell
singularity pull docker://hkubal/clair3:v0.1-r12
singularity pull docker://ontresearch/medaka:v1.7.2
```

Add the following parameters to the config file. See [cfg_cattle_ont](demo/configs/cfg_cattle_ont)
```
phase_clair3_command = singularity exec -B `pwd -P`:`pwd -P` -B /tmp:/tmp clair3_v0.1-r12.sif /opt/bin/run_clair3.sh
polish_medaka_command = singularity exec -B `pwd -P`:`pwd -P`  -B /tmp:/tmp medaka_v1.7.2.sif medaka
```
+ ``` `pwd -P`:`pwd -P` ```: It maps current working directory from the host to the container, so clair3 and medaka can access the files generated by PECAT.
+ ``` /tmp:/tmp ```: prevents `/tmp` in the container from becoming full.
+ `clair3_v0.1-r12.sif` and `medaka_v1.7.2.sif` should be replaced with the paths of corresponding images. Or the images are placed to the current path.

#### Using docker
Add the following parameters to the config file. 
```
phase_clair3_command =  docker run -i -v `pwd -P`:`pwd -P` -v /tmp:/tmp hkubal/clair3:latest /opt/bin/run_clair3.sh
polish_medaka_command = docker run -i -v `pwd -P`:`pwd -P` -v /tmp:/tmp ontresearch/medaka:v1.7.2 medaka
```
+ ``` `pwd -P`:`pwd -P` ```: see [Using singularity](#using-singularity).
+ ``` /tmp:/tmp ```： see [Using singularity](#using-singularity)

#### Without medaka and clair3
PECAT can run and achieve genome assembly without clair3 and medaka. Set the following parameters to the config file.
```
phase_method = 0
polish_medaka = 0
``` 
See [cfg_cattle_clr](demo/configs/cfg_cattle_clr)

## Docker pre-built image
There is a pre-built [docker image](https://hub.docker.com/r/lemene/pecat). Use the following commands to run pecat.

```
docker run -i -v $CWD:/mnt -v /var/run/docker.sock:/var/run/docker.sock lemene/pecat:v0.0.3 pecat.pl unzip cfg
```
+ ```-v $CWD:/mnt```: Map current working directory (`$CWD`) of the host to current working directory (`/mnt`) of the container, so PECAT can access the config file `cfg`.
+ ```-v /var/run/docker.sock:/var/run/docker.sock```: [Docker in Docker](https://devopscube.com/run-docker-in-docker/). By adding this parameter, pecat in the container can run the docker images (clair3 and medaka) of the host.
+ The directory of datasets should also be mapped carefully to ensure that PECAT in the container can access them.

Add the following parameters to the config file, so that pecat can call `clair3` and `medaka` in the container.
```
phase_clair3_command =  docker run -i -v $CWD:/mnt -v /tmp:/tmp hkubal/clair3:latest /opt/bin/run_clair3.sh
polish_medaka_command = docker run -i -v $CWD:/mnt -v /tmp:/tmp ontresearch/medaka:v1.7.2 medaka
```
+ `$CWD`: should be set to an absolute path of current working directory in the host.


### Using singularity
Download PECAT image.
```
singularity pull docker://lemene/pecat:v0.0.3
```
Run PECAT using the following command.
```
singularity exec -B `pwd -P`:`pwd -P` pecat_v0.0.3.sif pecat.pl unzip cfg
```
+ ``` `pwd -P`:`pwd -P` ```: It maps current working directory from the host to the container, so PECAT can access the config file `cfg`.
+ The directory of datasets should also be mapped carefully to ensure that PECAT in the container can access them.

We did not successfully run the singularity image in the container. It reports 
```ERROR  : Failed to create user namespace: user namespace disabled```.
So in this mode PECAT cannot run `clair3` and `medaka`. See [Without medaka and clair3](#without-medaka-and-clair3)


## Testing
We can run the demo to test whether PECAT has been succesfully installed. See [demo/README.md](demo/README.md).
```Shell
cd demo
pecat.pl unzip cfgfile
```

# Quick Start
Create a config file using the following command,

```shell
$ pecat.pl config cfgfile
```

Fill in the necessary parameters.

``` shell
project=S1
reads=./demo/reads.fasta.gz
genome_size=1500000
  ......
```

Run PECAT to assemble the reads.
```
$ pecat.pl unzip cfgfile
```

+ The corrected reads are in the file `S1/1-correct/corrected_reads.fasta`. 
+ The primary/alternate-format contigs are in the files `S1/6-polish/racon/{primary.fasta,alternate.fasta}`. 
+ The dual-format contigs are in the files `S1/6-polish/racon/{haplotype_1.fasta,haplotype_2.fasta}`. 
+ If the paramter `polish_medaka=1` is set, PECAT uses Medaka to further polish the above results, and the contigs are placed in `S1/6-polish/medaka`.

In the `demo` directory, there is a small example (`demo/{cfgfile,reads.fasta.gz}`) and several config files (`demo/configs`). When assembling a dataset, you can choose a config file of a similar species as a template and modify its parameters. See [config.md](doc/config.md). 

## Notes
***Note:*** We strongly recommend setting the parameter `cleanup=1`. PECAT deletes temporary files, otherwise it take up a lot of disk space.

***Note:*** For large genomes such as cattle and human, we strongly suggest adding the parameter `-f 0.005` or `-f 0.002` to `corr_rd2rd_options` and `align_rd2rd_options`. See [cfg_cattle_clr](demo/configs/cfg_cattle_clr), [cfg_cattle_ont](demo/configs/cfg_cattle_ont) and [cfg_hg002_ont](demo/configs/cfg_hg002_ont). The parameter is passed to `minimap2`, which means to filter out top 0.005 or 0.002 fraction of repetitive minimizers. It outputs less candidate overlaps, which reduces disk usage and speeds up error correction step and assembling step. 


# Resource usage
|Dataset|Size|Cov.|Config|CPU time|Peak memory usage|Peak disk space usage|
|:---------------|-----:|-----:|:---------------------|-------:|--------:|----:|
|[Yeast-CLR](demo/README.md)|12Mb   |200|[cfg_yeast_clr](demo/configs/cfg_yeast_clr)  |11h  |18G |4G |
|[Arab-CLR](demo/README.md)|130Mb   |135|[cfg_arab_clr](demo/configs/cfg_arab_clr)    |167h |71G |80G |
|[Dro-CLR](demo/README.md) |140Mb   |146|[cfg_dro_clr](demo/configs/cfg_dro_clr)      |142h |41G |49G |
|[Cattle-CLR](demo/README.md)|2.7Gb |135|[cfg_cattle_clr](demo/configs/cfg_cattle_clr)|4437h|219G|1099G|
|[Arab-ONT](demo/README.md) |130Mb  |106|[cfg_arab_ont](demo/configs/cfg_arab_ont)    |359h |179G|142G|
|[Cattle-ONT](demo/README.md) |2.7Gb|200|[cfg_cattle_ont](demo/configs/cfg_cattle_ont)|8869h|381G|1574G|
|[HG002-ONT](demo/README.md) |3Gb   |59|[cfg_hg002_ont](demo/configs/cfg_hg002_ont)   |7456h|348G|1211G|

The assemblies are available at https://doi.org/10.5281/zenodo.8380113

# More details
PECAT follows the correct-then-assemble strategy, including an error correction module and a two-round string-graph-based assembly module. Here, we describe some important steps and parameters. See [config.md](doc/config.md)

## Correcting raw reads
PECAT first extracts `prep_output_coverage` longest raw reads for correction. It uses minimap2 with `corr_rd2rd_options` to find the candidate overlaps between the extracted reads. PECAT corrects the raw reads with `corr_correct_options`. It implements `corr_iterate_number` rounds of error correction. After correcting, it extracts `corr_output_coverage` longest corrected reads for assembly, which are in the file `$PRJECT/1-correct/corrected_reads.fasta`. We can use the following scripts to correct raw reads.

```Shell
$ pecat.pl correct cfgile
```

## The first round of assembly
In the first round of assembly, PECAT first uses minimap2 with `align_rd2rd_options` to detect the overlaps between corrected reads. Minimap2 uses the seed-based method to find the overlaps, so the overlaps may have long overhangs. To reduce overhangs of overlaps, PECAT (`align_filter_options`) performs local alignment to extend overlaps to the ends of the reads and filter out the overlaps still with long overhangs. Then, PECAT (`asm1_assmeble_options`)  assembles the overlaps to haplotype-collapsed contigs. The contigs file is `$PROJECT/3-assemble/primary.fasta`. We can use the following scripts to run this step. 

```Shell
$ pecat.pl assemble cfgfile
```

## The second round of assembly
In the second round of assembly, PECAT first use minimap2 to map the reads (`phase_use_reads=0` for raw reads `phase_use_reads=1` for corrected reads) to `$PROJECT/3-assemble/primary.fasta` with `phase_rd2ctg_options`. PECAT calls the heterozygous SNP sites based on the base frequency of the alignments and identifies the inconsistent overlaps with `phase_phase_options`. PECAT removes the inconsistent overlaps with `phase_filter_options`. 

For Nanopore reads, we recommend using clair3 to call heterozygous SNPs from the raw reads. This is a similar process. You can use similar parameters above, but the parameters start with `phase_clair3_`.

After filtering out inconsistent overlaps, PECAT use `asm2_assemble_options` to assemble the filtered overlaps. The contigs files are placed in `$PROJECT/5-assemble`.


## Polishing
After generating the contigs, PECAT use minimap2 with `polish_map_options` to map reads  (`polish_use_reads=0` for raw reads `polish_use_reads=1` for corrected reads) to the `$PROJECT/5-assemble/{primary.fasta,alternate.fasta}` or `$PROJECT/5-assemble/haplotype_1.fasta,haplotype_2.fasta}` and uses racon with `polish_cns_options` to polish the contigs. The polished contigs are placed in `$PROJECT/6-polish/racon`.

If `polish_medaka=1` is set, PECAT use medaka to further improve the quality of the assembly. The parameters are similar and start with `polish_medaka_`. The contigs are placed in `$PROJECT/6-polish/medaka`.


## Running with multiple computation nodes
The pipeline script is written with [plgd](https://github.com/lemene/plgd). It supports PBS, SGE, LSF and Slurm systems. The follow parameter in the config file need to be set:
```shell
grid= auto:4
```

In the above example, `auto` means the pipeline automatically detects the type of cluster system. `pbs`, `sge`, `lsf` and `slurm` represent the corresponding systems, respectively. `4` computation nodes are used and each computation node run with `threads` CPU threads.



# Contact
+ Nie Fan, niefan@csu.edu.cn
+ QQ 316859622
