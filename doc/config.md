### 
+ `project`: project name. The pipeline makes a folder named `project` to contain all output files.
+ `reads`: the read file. It should be a `.fasta` or `.fastq` file. It can also be a `.txt` file which contains the *full paths* of all read files. 
+ `genome_size`: estimation of genome size.
+ `threads`: CPU threads for each job.
+ `memory`: memory for each job.
+ `cleanup`: whether to delete temporary files.
+ `compress`: whether to compress output files. The function is not completed.
+ `grid`: cluster system. It supports PBS, SGE, LSF and Slurm systems.
    + `auto`: automatically detecting the type of cluster system
    + `pbs`: PBS systems.
    + `sge`: SGE systems.
    + `lsf`: LSF systems.
    + `slurm`:Slurm systems.
    + `local`: do not use cluster system.

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
    + `aligner`: local alignment algorithm. 
    + `mini_coverage`: minimum coverage of base for correcting.
    + `score`: scoring method for selecting supporting reads. `weight` for diploid datasets. `count` for haploid datasets.
+ `corr_output_coverage`: coverage for extracting reads from corrected reads. The pipeline extracts `genome_size`*`corr_output_coverage` longest reads for the next steps.

### Algining
The pipeline finds the the overlaps between corrected reads.
+ `align_block_size`: the reads are split into multiple blocks for parallel processing.
+ `align_rd2rd_options`: parameters for `minimap2` to find the candiate overlaps between corrected reads.
+ `align_filter_options`: parameters for filtering low-quality overlaps and extend the overlaps to the ends of the reads. 
    + `filter0`: filter parameter beforing local alignments. (see `corr_filter_options`)
    + `task`: `extend` for performing local alignment on overhangs.
    + `filter1`: filter parameters aftering local alignment. (see `corr_filter_options`)

### First assembly
The pipeline generated haplotype-collapsed contigs.
+ `asm1_assemble_options`: parameters for the first assembly. 
    + `max_trivial_length`: pipeline uses the value to determine which edges in the complex region of the graph are trivial.
    + `min_identity`: minimum identity of the edges in the assembly graph.
    + `min_coverage`: minimum coverage of the reads. 
    + `min_contig_length`: mininum contig length for output.

### Phasing
The pipeline identifies the inconsistent overlaps. 
+ `phase_method`: the function is not completed.
+ `phase_use_reads`: using raw reads(`0`) or corrected reads(`1`) for phasing. 
+ `phase_rd2ctg_options`: parameters for `minimap2` to find the alignments bewteen the first assembly and the reads.
+ `phase_phase_options`: parameters for detecting inconsistent overlaps.
    + `coverage`: parameters for coverage for detecting SNP sites.
        + `lc`: minimum coverage.
    + `phase_options` parameters for dececting inconsistent overlaps.
        + `icr`: minimum ratio of different SNP alleles to all SNP alleles for identifing the pair of inconsistent reads.
        + `icc`: minimum different SNP alleles for identifing the pair of inconsistent reads.
        + `sc`: minimum support for a SNP allele.
+ `phase_filter_options`: parameters for filtering out inconsistent overlaps.
    + `threshold`: maxinum offset between overlaps detected by `minimap2` and inconsistent overlaps detected in the phasing step. `1000` for corrected reads and `-1` for raw reads.


### Second assembly
After removing the inconsistent overlaps, the pipeline generates primary/alterate-style contigs.
+ `asm2_assemble_options`: parameters for the second assembly. See `asm1_assemble_options`.

### Polishing
The pipeline polishes primary/alterate-style contigs.
+ `polish_use_reads`:  using raw reads(`0`) or corrected reads(`1`) for polishing. 
+ `polish_map_options`: parameters for `minimap2` to find the alignments bewteen the second assembly and the reads.
+ `polish_filter_options`: parameters for filtering low-quality alignments (see `corr_filter_options`). 
+ `polish_cns_options`: parameters for `racon` to polish contigs.



