
use strict;

use Plgd::Pipeline;

package FsaPipeline;
our @ISA = qw(Plgd::Pipeline);  


sub new { 
    my ($cls, $default) = @_; 
    my $self = $cls->SUPER::new($default); 

    bless $self, $cls; 
    return $self; 
} 


# copy reads to project and filter out some reads
sub getjob_prepare($$$$$) {
    my ($self, $name, $workDir, $ifile, $ofile) = @_;

    mkdir $workDir;
    
    my $isGz = ($ofile =~ /\.gz$/);
    my $ofileTemp = $ofile;
    $ofileTemp =~ s/\.gz//;

    my $binPath = $self->get_env("BinPath");
    my $threads = $self->get_config("THREADS");

    my $minLength = $self->get_config("PREP_MIN_LENGTH") + 0;
    my $baseSize = ($self->get_config("PREP_OUTPUT_COVERAGE") + 0) * ($self->get_config("GENOME_SIZE") + 0);
    my $id2name = ::dirname($ofile) . "/id2name.gz";

    my $job = $self->newjob(
        name => "${name}_job",
        ifiles => [$ifile],
        ofiles => [$ofile],
        gfiles => [$ofile],
        mfiles => [],
        threads => 1,
        cmds => ["$binPath/fsa_rd_tools longest  $ifile $ofileTemp --base_size $baseSize --min_length $minLength --id2name $id2name"],
        msg => "preparing reads",
    );

    # if ($isGz) {
    #     push @{$job->cmds}, "$binPath/pigz -f -p $threads $ofileTemp";
    # }

    return $job

}


sub run_prepare($) {
    my ($self) = @_;

    my $name = "prp";
    my $wrkdir = $self->get_work_folder("0-prepare");

    my $ifile = $self->get_config("reads");

    my $isGz = $self->get_config("compress");
    my $ofile = $isGz ? "$wrkdir/prepared_reads.fasta.gz" : "$wrkdir/prepared_reads.fasta";

    $self->run_jobs($self->getjob_prepare($name, $wrkdir, $ifile, $ofile));
}


sub getjob_make_softlink($$$$) {
    my ($self, $name, $ifile, $ofile) = @_;
    
    my $job = $self->newjob (
        name => "${name}_skip",
        ifiles => [$ifile],
        ofiles => [$ofile],
        mfiles => [],
        funcs => [sub ($$) {
            `ln -s -f $ifile $ofile`;
        }],        
        msg => "make softlink, $name",
    );
    return $job;
}

## 
sub getjob_extract_longest_reads($$$$$) {
    my ($self, $name, $ifname, $ofname, $bsize) = @_;
    
    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");

    my $job = $self->newjob(
        name => "${name}_extract",
        ifiles => [$ifname],
        ofiles => [$ofname],
        gfiles => [$ofname],
        mfiles => [],
        threads => 1,
        cmds => ["$bin_path/fsa_rd_tools longest $ifname $ofname --base_size=$bsize"],
        msg => "extracting longest reads, $name"
    );
    return $job;
}


sub getjob_map_read_to_read_0($$$$$$) {
    my ($self, $name, $reads, $options, $rd_2_rd, $wrkdir) = @_;

    my $isgz = ($rd_2_rd =~ /\.gz$/);

    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");

    my $opts_map = $options->[0];
    my $opts_flt = $options->[1];

    my $cmdSub = $isgz ? " | $bin_path/pigz -c -p $threads > $rd_2_rd" : " > $rd_2_rd";
    my $job = $self->newjob(
        name => "${name}_rd_2_rd_0",
        ifiles => [$reads], 
        ofiles => [$rd_2_rd], 
        mfiles => [],
        cmds => ["minimap2 -t $threads $opts_map $reads $reads | $bin_path/fsa_ol_cut - - --itype paf --otype paf --thread_size 4 $opts_flt $cmdSub"],
    
        msg => "mapping reads",
    );
    return $job;
}


sub jobRead2ReadParallelly($$$$$$$) {
    my ($self, $name, $workDir, $reads, $rd2rd, $options, $blockSize) = @_;


    mkdir $workDir;

    my $isGz = ($rd2rd =~ /\.gz$/);

    my $binPath = $self->get_env("BinPath");
    my $threads = $self->get_config("THREADS");

    my $blockInfo = "$workDir/al_block_info";
    my $align_options = $options->[0];
    my $refine_options = $options->[1];

    my $readsIndex = "$workDir/reads.index";
    my $jobIndex = $self->newjob(
        name => "${name}_rd2rd_index",
        ifiles => [$reads],
        ofiles => [$readsIndex],
        gfiles => [$readsIndex],
        mfiles => [],
        cmds => ["minimap2 -t $threads $align_options -d $readsIndex $reads"],
        msg => "indexing reads, $name"
    );

    my $jobSplit = $self->newjob(
        name => "${name}_rd2rd_split",
        ifiles => [$reads],
        ofiles => [$blockInfo],  
        gfiles => [$blockInfo, "$workDir/cc*.fasta"],
        mfiles => [],
        threads => 1,
        cmds => ["$binPath/fsa_rd_tools split $reads $workDir/cc{}.fasta --block_size=$blockSize",
                 "ls $workDir/cc*.fasta > $blockInfo"],
        msg => "spliting reads, $name",
    );

    my $jobAlign = $self->newjob(
        prefunc => sub ($) {
            my ($job) = @_;
            @{$job->{ifiles}} = @{$jobSplit->{ofiles}};
            my $size = `wc -l $blockInfo`;

            for (my $i=0; $i < $size; $i=$i+1) {
                my $f = "$workDir/cc${i}.fasta";
                
                #my $ofile = $isGz ? "$f.paf.gz" : "$f.paf";
                my $ofile = "$f.paf";
                my $jobSub = $self->newjob(
                    name => "${name}_rd2rd_align_$i",
                    ifiles => [$blockInfo],
                    ofiles => [$ofile],
                    gfiles => [$ofile],
                    mfiles => [$f],
                    cmds => ["minimap2 $align_options -t $threads $readsIndex $f | $binPath/fsa_ol_refine - - --itype paf --otype paf  --read_fname $reads --thread_size $threads $refine_options > $ofile"],
                    msg => "aligning reads $i, $name"
                );

                push @{$job->{pjobs}}, $jobSub;
                push @{$job->{ofiles}}, $ofile;
            }
        my $n = scalar @{$job->{pjobs}};

        },
        name => "${name}_rd2rd_align",
        ifiles => [$blockInfo], #prefunc
        ofiles => [], # prefunc
        mfiles => [],
        pjobs => [], # prefunc
        msg => "aligning reads, $name",
    );


    my $is_list = ($rd2rd =~ /\.txt$/);
    my $jobCat = $self->newjob(
        prefunc => sub ($) {
            my ($job) = @_;


            @{$job->{ifiles}} = @{$jobAlign->{ofiles}};
            if ($is_list) {
                my $s = join("\n", @{$jobAlign->{ofiles}});
                $job->{funcs} = [sub ($$) { Plgd::Utils::echo_file($rd2rd, $s); }];
            } else {
                $job->{cmds} = ["cat @{$jobAlign->{ofiles}} > $rd2rd && rm @{$jobAlign->{ofiles}}"];
            }

        },
        name => "${name}_rd2rd_cat",
        ifiles => [], #prefunc
        ofiles => [$rd2rd], 
        gfiles => [$rd2rd],
        mfiles => [],
        #cmds => [],     # prefunc
        funcs => [],
        threads => 1,
        msg => "cat overlaps, $name",
    );

    return $self->newjob(
        name => "${name}_rd2rd_job",
        ifiles => [$reads],
        ofiles => [$rd2rd],
        mfiles => ["$workDir/cc*.fasta", $readsIndex, $blockInfo],
        jobs => [$jobSplit, $jobIndex, $jobAlign, $jobCat],
        msg => "aligning reads to reads, $name",

    );
}

sub job_map_reads_to_contigs($$$$$$) {
    my ($self, $name, $wrkdir, $reads, $contigs, $options) = @_;

    my $threads = $self->get_config("threads");

    my $prictg = $contigs->[0];
    my $read2ctg = "$wrkdir/rd2ctg.paf";

    my $job = $self->newjob(
        name => "${name}_rd2ctg",
        ifiles => [$reads, $prictg],
        ofiles => [$read2ctg],
        gfiles => [$read2ctg],
        mfiles => [],
        cmds => ["minimap2  -t $threads $options $prictg $reads > $read2ctg"],
        msg => "mapping reads to contigs, ${name}",
    );

    return $job;
}


sub getjob_map_read_to_reference($$$$$$$$) {
    my ($self, $name, $reads, $ref, $rd_2_ref, $options, $wrkdir) = @_;

    my $threads = $self->get_config("threads");
    my $job = $self->newjob(
        name => "${name}_rd_2_ref",
        ifiles => [$reads, $ref],
        ofiles => [$rd_2_ref],
        gfiles => [$rd_2_ref],
        mfiles => [],
        cmds => ["minimap2  -t $threads $options $ref $reads > $rd_2_ref"],
        msg => "mapping reads to reference, ${name}",
    );

    return $job;
}

sub getjob_map_read_to_ref_sam($$$$$$$$) {
    my ($self, $name, $reads, $ref, $rd_2_ref, $options, $wrkdir) = @_;

    my $threads = $self->get_config("threads");
    my $job = $self->newjob(
        name => "${name}_rd_2_ref",
        ifiles => [$reads, $ref],
        ofiles => [$rd_2_ref],
        gfiles => [$rd_2_ref],
        mfiles => [],
        cmds => ["minimap2 -t $threads $options $ref $reads -a | samtools sort -@ $threads > $rd_2_ref",
                 "samtools index -@ $threads $rd_2_ref"],
        msg => "mapping reads to reference(sam), ${name}",
    );
}

sub job_map_hic_reads_to_contigs($$$$$$) {
    my ($self, $name, $wrkdir, $reads, $contigs, $options) = @_;

    my $threads = $self->get_config("threads");

    my $prictg = $contigs->[0];
    my $hic1_reads = $reads->[0];
    my $hic2_reads = $reads->[1];

    my $hic1_2_ctg = "$wrkdir/hic1_2_ctg.paf";
    my $hic2_2_ctg = "$wrkdir/hic2_2_ctg.paf";

    my $job1 = $self->newjob(
        name => "${name}_hic1_2_ctg",
        ifiles => [$hic1_reads, $prictg],
        ofiles => [$hic1_2_ctg],
        gfiles => [$hic1_2_ctg],
        mfiles => [],
        cmds => ["minimap2  -t $threads -x sr -c $prictg $hic1_reads > $hic1_2_ctg"],
        msg => "mapping hic1 reads to contigs, ${name}_hic1",
    );

    my $job2 = $self->newjob(
        name => "${name}_hic2_2_ctg",
        ifiles => [$hic2_reads, $prictg],
        ofiles => [$hic2_2_ctg],
        gfiles => [$hic2_2_ctg],
        mfiles => [],
        cmds => ["minimap2  -t $threads -x sr -c $prictg $hic2_reads > $hic2_2_ctg"],
        msg => "mapping hic2 reads to contigs, ${name}_hic2",
    );

    #return [$job1, $job2];
    return $self->newjob(
        name => "${name}_hic_2_ctg",
        ifiles => [$hic1_reads, $hic2_reads, $prictg], # [$reads, $prictg, $altctg],
        ofiles => [$hic1_2_ctg, $hic2_2_ctg], #TODO 有两份独立依赖关系
        mfiles => [],
        jobs => [$job1, $job2],
        msg => "mapping hic reads to contigs, ${name}_hic2",
    );

}

sub job_identify_snps_in_hic($$$$) {
    
    my ($self, $name, $wrkdir, $reads) = @_;

    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");

    my $hic1_reads = $reads->[0];
    my $hic2_reads = $reads->[1];

    my $hic1_2_ctg = "$wrkdir/hic1_2_ctg.paf";
    my $hic2_2_ctg = "$wrkdir/hic2_2_ctg.paf";
    my $variants = "$wrkdir/fsa/variants";
    my $snp_in_hic = "$wrkdir/hic_infos";

    my $job = $self->newjob(
        name => "${name}_snp_in_hic",
        ifiles => [$hic1_reads, $hic2_reads, $hic1_2_ctg, $hic2_2_ctg, $variants ],
        ofiles => [$snp_in_hic],
        gfiles => [$snp_in_hic],
        mfiles => [],
        cmds => ["$bin_path/fsa_misc_tools snp_in_hic $hic1_reads $hic1_2_ctg $hic2_reads $hic2_2_ctg $variants $snp_in_hic"],
        msg => "identify snps in hic reads, ${name}_snp_in_hic",
    );


    return $job;
}

sub job_phase($$$$$$$) {
    my ($self, $name, $wrkdir, $reads, $contigs, $ol_r2c, $options) = @_;
 
    my $inconsistent = "$wrkdir/inconsistent";

    my $binPath = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");

    my $job = $self->newjob(
        name => "${name}_phase",
        ifiles => [$reads, $contigs, $ol_r2c],
        ofiles => [$inconsistent],
        gfiles => [$inconsistent],
        mfiles => [],
        cmds => ["$binPath/fsa_rd_haplotype $options --thread_size=$threads --ctg_fname=$contigs --rd_fname=$reads --ol_fname=$ol_r2c --output_directory=$wrkdir"],
        msg => "phasing reads with contigs, $name",
    );

    return $job;
}


sub runAssemble($$$$$$$) {
    my ($self, $name, $workDir, $reads, $overlaps, $options) = @_;

    mkdir $workDir;

    my $wrkdir = $workDir;

    my $primary = "$workDir/primary.fasta";
    my $alternate = "$workDir/alternate.fasta";
    my $haplotype_1 = "$workDir/haplotype_1.fasta";
    my $haplotype_2 = "$workDir/haplotype_2.fasta";

    my $hic_info = $self->get_config("hic_reads") ne "" ? "--hic_info $wrkdir/../4-phase/hic_infos" : "";

    my $phased = "";
    my $variants = "";
    if ($name eq "asm2") {
        
        my $phase_method = $self->get_config("phase_method") + 0;
        if ($phase_method == 0) {
            $phased = "--phased $workDir/../4-phase/fsa/inconsistent";
            $variants = "--variants $workDir/../4-phase/fsa/readinfos";
        } elsif ($phase_method == 1 || $phase_method == 2) {
            $phased = "--phased $workDir/../4-phase/clair3/inconsistent";
            $variants = "--variants $workDir/../4-phase/clair3/readinfos";
        }
    }
 
    my $dual = index($options, "dual") != -1;
    my $prialt = index($options, "prialt") != -1 || (!$dual);
    my @contigs = ();
    if ($dual) {
        push @contigs, $haplotype_1;
    }
    if ($prialt) {
        push @contigs, $primary;
    }
    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("THREADS");

    my $job = $self->newjob(
        name => "${name}_assemble",
        ifiles => [$overlaps],
        ofiles => [@contigs],
        gfiles => [$primary, $alternate, $haplotype_1, $haplotype_2],
        mfiles => ["$wrkdir/debug_*", "$wrkdir/filter*", "$wrkdir/readinfos", "$wrkdir/*.gz"],
        cmds => ["$bin_path/fsa_ol_assemble $overlaps --thread_size=$threads --output_directory=$workDir --read_file=$reads  $phased $variants $hic_info $options"],
        msg => "assembling overlaps, $name",
    );

    $self->run_jobs($job);
}

sub job_polish($$$$$$$) {
    my ($self, $name, $wrkdir, $contigs, $reads, $options, $sub) = @_;

    mkdir $wrkdir;

    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");

    my $map_options = $options->[0];
    my $polish_options = $options->[1];

    my $rd2ctg = "$wrkdir/rd2ctg$sub.paf";
    my $rd2ctg_flt = "$wrkdir/rd2ctg.flt$sub.paf";
    my $polished = "$wrkdir/polished$sub.fasta";

    my $job = $self->newjob(
        name => "${name}_polish$sub",
        ifiles => [$contigs, $reads],
        ofiles => [$polished],
        gfiles => [$polished],
        mfiles => [],
        cmds => ["minimap2  -t $threads $map_options $contigs $reads > $rd2ctg",
                 "$bin_path/fxtools.py fx_purge_paf $rd2ctg $wrkdir/contig_tiles$sub $wrkdir/id2name.txt.gz $wrkdir/../4-phase/phased > $rd2ctg_flt",
                 "racon -t $threads $reads $rd2ctg_flt $contigs > $polished"],
        msg => "polishing contigs$sub, $name",
    );

    return $job;

}

1;
