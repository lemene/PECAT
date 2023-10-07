### 
+ `project`: project name. The pipeline makes a folder named `project` to contain all output files.
+ `reads`: the read file. It should be a `.fasta` or `.fastq` file. It can also be a `.txt` file which contains the *full paths* of all read files. 
+ `genome_size`: estimation of genome size.
+ `threads`: CPU threads for each job.
+ `memory`: memory for each job.
+ `cleanup`: `1` delete temporary files, `0` do not delete temporary files
+ `grid`: cluster system. It supports PBS, SGE, LSF and Slurm systems.
    + `auto`: automatically detecting the type of cluster system.
    + `pbs`: PBS systems.
    + `sge`: SGE systems.
    + `lsf`: LSF systems.
    + `slurm`:Slurm systems.
    + `local`: do not use cluster system.
+ `grid_options`: It is used to add additional options for submitting jobs.

### Preparing
The pipeline extracts `genome_size`*`prep_output_coverage` longest reads from `reads` for the next steps.
+ `prep_min_length`: minimum length of reads
+ `prep_output_coverage`: coverage of extracted reads from raw reads.

### Correcting
Long nosiy reads are corrected in this step.
+ `corr_iterate_number`: number of rounds for correction.
+ `corr_block_size`: the reads are split into multiple blocks for parallel processing.
+ `corr_rd2rd_options`: parameters for `minimap2` to find the candiate overlaps between raw reads. 
+ `corr_filter_options`: parameters for filtering low-quality overlaps.
    + `filter0`: filter parameters.
        + `l`: minimun read length
        + `al`: minimum overlap length
        + `alr`: minimum ratio of overlap length to read length. 
        + `oh`: maximum overhang length.
        + `ohr`: maximum ratio of overhang length to read length. 
        + `aal`: the overlap is not filtered if the length exceeds this value.
        + `aalr`: the overlap is not filtered if ratio of overlap length to read length exceeds this value.
+ `corr_correct_options`: parameters for correcting raw reads. 
    + `aligner`: local alignment algorithm. It can be set to `diff` and `edlib`.
    + `min_coverage`: minimum coverage of base for correcting.
    + `score`: scoring method for selecting supporting reads. `weight` for diploid datasets. `count` for haploid datasets. `lc`: minimum coverage for detecting haplotype branch in the POA graph.
    + `candidate`: parameters for selecting candidate overlaps for local alignment. `n`: maxinum number for implementing local alignment. `f`:  maximum number of consecutive local alignment failures.
    + `min_identity`: minimum identity of alignment between the tempalte read and the supporting read.
    + `min_local_identity`: minimum local identity where the local window length is 1000.
    + `filter0`: parameters for filtering overlaps before implementing local alignment. See `corr_filter_options`.
    + `filter1`: parameters for filtering overlaps after implementing local alignment. See `corr_filter_options`.

+ `corr_output_coverage`: coverage for extracting reads from corrected reads. The pipeline extracts `genome_size`*`corr_output_coverage` longest reads for the next steps.

### Algining
The pipeline finds the the overlaps between corrected reads.
+ `align_block_size`: the reads are split into multiple blocks for parallel processing.
+ `align_rd2rd_options`: parameters for `minimap2` to find the candiate overlaps between corrected reads.
+ `align_filter_options`: parameters for filtering low-quality overlaps and extend the overlaps to the ends of the reads. 
    + `filter0`: filter parameter beforing local alignments. (see `corr_filter_options`)
    + `task`: `extend` for performing local alignment on overhangs.
    + `filter1`: filter parameters aftering local alignment. (see `corr_filter_options`)

### The first round of assembly
The pipeline generated haplotype-collapsed contigs.
+ `asm1_assemble_options`: parameters for the first assembly. 
    + `max_trivial_length`: pipeline uses the value to determine which edges in the complex region of the graph are trivial.
    + `min_identity`: minimum identity of the edges in the assembly graph.
    + `min_coverage`: minimum coverage of the reads. 
    + `min_contig_length`: mininum contig length for output.
    + `contig_format`: `prialt` output primary/alternate format contigs, `dual` output dual format contigs, `prialt,dual` output both format contigs. 
    + `reducer0`: parameters for graph reduction
         `best`: best overlap graph algorithm. `cmp`: parameters for comparing two edges.
         `phase`: using SNP alleles to phase adjacent bubble structures. `sc`: minimum support for a SNP allele.
    + `reducer1`: parameters for graph reduction
        + `spur`: removing spur paths. `length` maximum spur length, `nodesize`: maximum node size in graph of spur

### Phasing
The pipeline identifies the inconsistent overlaps. The reads are mapped to the haplotype-collapsed contigs. PECAT use built-in method or clair3 to identity the SNPs. According to SNP allels in reads, PECAT identifies and remove inconsistent overlaps.
+ `phase_method`: using built-in method or clair3. `0` built-in method, `1` clair3, `2` both methods.
+ `phase_use_reads`: using raw reads(`0`) or corrected reads(`1`) for phasing. 
+ `phase_rd2ctg_options`: parameters for `minimap2` to find the alignments bewteen the haplotype-collapsed contigs and the reads.
+ `phase_phase_options`: parameters for detecting inconsistent overlaps.
    + `coverage`: parameters for coverage for detecting SNP sites.
        + `lc`: minimum coverage.
    + `phase_options` parameters for dececting inconsistent overlaps.
        + `icr`: minimum ratio of different SNP alleles to all SNP alleles for identifing the pair of inconsistent reads.
        + `icc`: minimum different SNP alleles for identifing the pair of inconsistent reads.
        + `sc`: minimum support for a SNP allele.
+ `phase_filter_options`: parameters for filtering out inconsistent overlaps.
    + `threshold`: maxinum offset between overlaps detected by `minimap2` and inconsistent overlaps detected in the phasing step. `1000` for corrected reads and `-1` for raw reads.
    + `rate`: maximum rate of offset to read length.

+ `phase_clair3_command`: if using `singularity`, it can be set `singularity exec --containall -B ``pwd -P``:``pwd -P`` clair3_v0.1-r12.sif /opt/bin/run_clair3.sh`
+ `phase_clair3_options`: parameters passed to clair3.
+ `phase_clair3_rd2ctg_options`: similar to `phase_rd2ctg_options`.
+ `phase_clair3_phase_options`: similar to `phase_phase_options`.
+ `phase_clair3_use_reads`: similar to `phase_use_reads`.
+ `phase_clair3_filter_options`: similar to `phase_filter_options`.


### The second round of assembly
After removing the inconsistent overlaps, the pipeline reassembles the corrected reads again.
+ `asm2_assemble_options`: parameters for the second round of assembly. See `asm1_assemble_options`.

### Polishing
By default, we only use `racon` to polish contigs. If `polish_medaka=1` is set, we further polish the output of racon with `medaka`.

Parameters for using `racon`.
+ `polish_use_reads`:  using raw reads(`0`) or corrected reads(`1`) for polishing. 
+ `polish_map_options`: parameters for `minimap2` to find the alignments bewteen the second assembly and the reads.
+ `polish_filter_options`: parameters for filtering low-quality alignments (see `corr_filter_options`). 
+ `polish_cns_options`: parameters for `racon` to polish contigs.

Parameters for using `medaka`.
+ `polish_medaka`: whether to use medaka to polish contigs.
+ `polish_medaka_command`: medaka command. If using `singularity`, it can be set ```singularity exec --containall -B `pwd -P`:`pwd -P` medaka_1.7.2--aa54076.sif medaka``` 
+ `polish_medaka_map_options`: parameters for `minimap2` to find the alignments bewteen the outputs of `racon` and the reads.
+ `polish_medaka_cns_options`: parameters for `medaka`.


