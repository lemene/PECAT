#!/usr/bin/env perl

use FindBin;
use lib $FindBin::RealBin;

use Cwd;
use File::Basename;
use Carp;
use POSIX;

use Plgd::Utils;
use Plgd::Logger;
use Plgd::Config;
use Plgd::Pipeline;
use Plgd::Job;
use FsaUtils;

use Env qw(PATH);

use strict;

package FsaPipeline;

our @ISA = qw(Plgd::Pipeline);  

sub new { 
    my ($cls, $default) = @_; 
    my $self = $cls->SUPER::new($default); 

    bless $self, $cls; 
    return $self; 
} 

sub initialize($$) {
    my ($self, $fname) = @_;

    $self->SUPER::initialize($fname);
    
}

sub get_fsa_wrkdir($$) {
    my ($self, $name) = @_;

    if ($name eq "prp") {
        return $self->get_work_folder("0-prepare");
    } elsif ($name eq "crr") {
        return $self->get_work_folder("1-prepare");
    } elsif ($name eq "al") {
        return $self->get_work_folder("2-align");
    } elsif ($name eq "asm1") {
        return $self->get_work_folder("3-assemble");
    } elsif ($name eq "phs") {
        return $self->get_work_folder("4-phase");
    } elsif ($name eq "asm2") {
        return $self->get_work_folder("5-assemble");
    } elsif ($name eq "pol") {
        return $self->get_work_folder("6-polish");
    } elsif ($name eq "gcpp") {
        return $self->get_work_folder("7-pbgcpp");
    } else {
        Plgd::Logger::error("There is no step:  $name");
    }
}

# copy reads to project and filter out some reads
sub job_prepare($$$$$) {
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

    if ($isGz) {
        push @{$job->cmds}, "$binPath/pigz -f -p $threads $ofileTemp";
    }

    return $job

}

sub run_prepare($) {
    my ($self) = @_;

    my $name = "prp";
    my $wrkdir = $self->get_work_folder("0-prepare");

    my $ifile = $self->get_config("reads");

    my $isGz = $self->get_config("compress");
    my $ofile = $isGz ? "$wrkdir/prepared_reads.fasta.gz" : "$wrkdir/prepared_reads.fasta";

    $self->run_jobs($self->job_prepare($name, $wrkdir, $ifile, $ofile));
}


sub run_align($) {
    my ($self) = @_;
 
    my $name = "al";
    my $wrkdir = $self->get_work_folder("2-align");
    
    my $workdir_crr = $self->get_work_folder("1-correct");
    my $isGz = $self->get_config("compress");

    my $corrReads = $isGz ? "$workdir_crr/corrected_reads.fasta.gz" : "$workdir_crr/corrected_reads.fasta";
    my $overlaps = "$wrkdir/overlaps.txt";

    mkdir $wrkdir;

    $self->run_jobs($self->jobRead2ReadParallelly($name, $wrkdir, $corrReads, $overlaps,
                        [$self->get_config("ALIGN_RD2RD_OPTIONS"), $self->get_config("ALIGN_FILTER_OPTIONS")],  
                        $self->get_config("ALIGN_BLOCK_SIZE")));
}

sub run_assemble1($) {
    my ($self,) = @_;
    
    my $name = "asm1";
    my $wrkdir = $self->get_work_folder("3-assemble");
    
    my $wrkdir_al = $self->get_work_folder("2-align");
    my $workdir_crr = $self->get_work_folder("1-correct");
    my $isGz = $self->get_config("compress");

    my $reads = $isGz ? "$workdir_crr/corrected_reads.fasta.gz" : "$workdir_crr/corrected_reads.fasta";
    my $overlaps = "$wrkdir_al/overlaps.txt";
        
    $self->runAssemble($name, $wrkdir, $reads, $overlaps, [$self->get_config("ASM1_FILTER_OPTIONS"), $self->get_config("ASM1_ASSEMBLE_OPTIONS")]);
}


sub jobCorrect($$$$$) {
    my ($self, $name, $rawreads, $corrected, $workDir) = @_;
    
    
    my $isGz = ($corrected =~ /\.gz$/);
    
    my $rd2rd = "$workDir/rd2rd.txt";

    my $binPath = $self->get_env("BinPath");
    my $threads = $self->get_config("THREADS");
    my $readName = "$workDir/readname";

    my $step = substr($name, 3, length($name)-3);   #  ”crr0 crr1"

    my $blockSize = $self->get_config("CORR_BLOCK_SIZE");
    my $blockInfo = "$workDir/block_info";
    my $base_size = $self->get_config("genome_size") * $self->get_config("prep_output_coverage") * 1.3;

    my $rd2rdOptions = $self->get_config2("CORR${step}_RD2RD_OPTIONS", "CORR_RD2RD_OPTIONS");
    my $filterOptions = $self->get_config2("CORR${step}_FILTER_OPTIONS", "CORR_FILTER_OPTIONS");
    my $correctOptions = $self->get_config2("CORR${step}_CORRECT_OPTIONS","CORR_CORRECT_OPTIONS");
 
    my $jobCand = $self->jobRead2ReadParallelly($name, $workDir, $rawreads, $rd2rd, [$rd2rdOptions, $filterOptions], $blockSize);

    my $jobSplit = $self->newjob(
        name => "${name}_split",
        ifiles => [$rawreads,$rd2rd],
        ofiles => [$blockInfo],
        gfiles => [$blockInfo, "$readName.*"],
        mfiles => [],
        # cmds => ["rm -rf $readName.*", 
        #          #"$binPath/fsa_misc_tools split_ols $rd2rd $rd2rd.sub.{}.paf --rdfname0 $readName.core.{} --rdfname1 $readName.all.{} --block_size $blockSize",
        #          "$binPath/fsa_rd_tools split_name $rawreads $readName.core.{}  --block_size $blockSize --base_size $base_size",
        #          "ls $readName.core.* > $blockInfo",
        #          "SUBSIZE=`wc -l $blockInfo | awk '{print \$1}'` && $binPath/fsa_misc_tools split_ols2 $rd2rd $rd2rd.sub.{}.paf $rawreads --rdfname0 $readName.core.{} --rdfname1 $readName.all.{} --sub_size \$SUBSIZE --thread_size 10",
        #          #"$binPath/fsa_ol_purge $rd2rd $readName.{} --read_file $rawreads --thread_size $threads",
        #          ],
        cmds => ["rm -rf $readName.*",
                 "$binPath/fsa_rd_tools split_name $rawreads $readName.core.{}  --block_size $blockSize --base_size $base_size --overlaps  $rd2rd  --sub_overlaps $rd2rd.sub.{}.paf --thread_size $threads",
                 "ls $readName.core.* > $blockInfo"],
        msg => "spliting read names, $name",
    );

    my $jobCorr = $self->newjob(
        prefunc => sub($) {
            my ($job) = @_;
            my $size = `wc -l $blockInfo`;
            for (my $i=0; $i < $size; $i=$i+1) {

                my $corrSub = $isGz ? "$corrected.$i.gz" : "$corrected.$i";
                
                my $jobSub = $self->newjob(
                    name => "${name}_correct_$i",
                    ifiles => [$rawreads, $rd2rd, "$readName.core.$i"],
                    ofiles => [$corrSub],
                    gfiles => [$corrSub],
                    mfiles => [],
                    cmds => ["$binPath/fsa_rd_correct $rd2rd.sub.$i.paf $rawreads $corrected.$i --output_directory=$workDir --thread_size=$threads " . 
                                " --read_name_fname=$readName.core.$i --infos_fname $corrected.$i.infos $correctOptions"],
                    msg => "correcting reads $i, $name"
                );
                if ($isGz) {
                    push @{$jobSub->cmds}, "$binPath/pigz -p $threads $corrected.$i";
                }
                push @{$job->{ofiles}}, $corrSub;
                push @{$job->{pjobs}}, $jobSub;
            }

        },
        name => "${name}_correct_all",
        ifiles => [$blockInfo, $rawreads, $rd2rd],
        ofiles => [],                   # prefunc
        mfiles => [],
        pjobs => [],                    # prefunc
        msg => "correcting rawreads, $name",
    );

    my $jobCat = $self->newjob(
        prefunc => sub($) {
            my ($job) = @_;
            my $size = `wc -l $blockInfo`;

            my @correctedSub = ();
            for (my $i=0; $i < $size; $i=$i+1) {
                $correctedSub[$i] = $isGz ? "$corrected.$i.gz" : "$corrected.$i";
            }

            push @{$job->{ifiles}}, @correctedSub;

            push @{$job->{cmds}}, "cat @correctedSub > $corrected && rm @correctedSub";


        },
        name => "${name}_cat",
        ifiles => [],      # prefunc
        ofiles => [$corrected], 
        gfiles => [$corrected], 
        mfiles => [],
        cmds => [],                     # prefunc
        threads => 1,
        msg => "cat corrected reads, $name",

    );
    
    return $self->newjob(
        name => "${name}_correct",
        ifiles => [$rawreads],
        ofiles => [$corrected], # prefunc
        mfiles => ["$rd2rd.sub.*.paf", "$corrected.*.gz", "$corrected.*"],
        jobs => [$jobCand, $jobSplit, $jobCorr, $jobCat],
        msg => "correcting rawreads, $name");
    
}

sub run_correct($) {
    my ($self) = @_;

    my $name = "crr";
    my $wrkdir = $self->get_work_folder("1-correct");

    my $isGz = $self->get_config("compress");
    my $wrkdir_prp = $self->get_work_folder("0-prepare");

    my $reads = $isGz ? "$wrkdir_prp/prepared_reads.fasta.gz" : "$wrkdir_prp/prepared_reads.fasta";
    my $corrReads = $isGz ? "$wrkdir/corrected_reads.fasta.gz" : "$wrkdir/corrected_reads.fasta";

    mkdir $wrkdir;
    
    my @jobs = ();

    my $baseSize = ($self->get_config("CORR_OUTPUT_COVERAGE") + 0) * ($self->get_config("GENOME_SIZE") + 0);
    my $iterNum = $self->get_config("CORR_ITERATE_NUMBER") + 0;

    my $corrInput = $reads;
    my $corrOutput = $reads;

    for (my $i=0; $i<$iterNum; $i=$i+1) {
        $corrOutput = $isGz ? "$wrkdir/$i/corrected_reads.fasta_$i.gz" : "$wrkdir/$i/corrected_reads_$i.fasta";
        push @jobs, $self->jobCorrect("crr" . $i, $corrInput, $corrOutput, "$wrkdir/$i");
        $corrInput = $corrOutput;
    }
    
    if ($baseSize > 0) {
        push @jobs, $self->jobExtract($name, $corrOutput, $corrReads, $baseSize);  
    } else {
        push @jobs, $self->jobSkip("crr", $corrOutput, $corrReads);
    }

    $self->run_jobs($self->newjob(
        name => "${name}_job",
        ifiles => [$reads],
        ofiles => [$corrReads],
        mfiles => [],
        jobs => [@jobs],
        msg => "correcting reads, $name",
    ));

}


sub job_filter_inconsistent($$$$$) {
    my ($self, $name, $wrkdir, $overlaps, $filtered) = @_;
 
    
    my $binPath = $self->get_env("BinPath"); 
    my $threads = $self->get_config("threads");
    my $consistent = "$wrkdir/consistent";
    my $inconsistent = "$wrkdir/inconsistent";
    my $threshold = $self->get_config("phase_use_reads") + 0 == 1 ? "" : "--threshold=-1";
    my $options = $self->get_config("phase_filter_options");

    my $job = $self->newjob(
        name => "${name}_filter",
        ifiles => [$overlaps, $consistent, $inconsistent],
        ofiles => [$filtered],
        gfiles => [$filtered],
        mfiles => [],
        #cmds => ["$binPath/fsa_ol_tools filter $overlaps $filtered --thread_size=$threads --consistent=$consistent --inconsistent=$inconsistent $threshold $options"],
        cmds => ["$binPath/fsa_ol_tools filter $overlaps $filtered --thread_size=$threads --inconsistent=$inconsistent $threshold $options"],
        msg => "filtering inconsistent overlaps, ${name}",
    );

    return $job;
}

sub job_phase_method0($$$) {
    my ($self, $name, $wrkdir) = @_;


    mkdir $wrkdir;

    my $isGz = $self->get_config("COMPRESS");
    my $prjDir = $self->get_project_folder();

    my $corrReads = $isGz ? "$prjDir/1-correct/corrected_reads.fasta.gz" : "$prjDir/1-correct/corrected_reads.fasta";
    my $prepReads = $isGz ? "$prjDir/0-prepare/prepared_reads.fasta.gz" : "$prjDir/0-prepare/prepared_reads.fasta";

    my $reads = $self->get_config("phase_use_reads") + 0 == 1 ? $corrReads : $prepReads;
    my $prictg = "$prjDir/3-assemble/primary.fasta";
    my $altctg = "$prjDir/3-assemble/alternate.fasta";
    my $rd2ctg = "$wrkdir/rd2ctg.paf";

    my $rd2rd = "$prjDir/2-align/overlaps.txt";
    my $rd2rd_flt = "$wrkdir/filtered_overlaps.paf";

    my $job_map = $self->job_map_reads_to_contigs("phs", $wrkdir, $reads, [$prictg, $altctg], $self->get_config("PHASE_RD2CTG_OPTIONS"));
    my $job_phase = $self->job_phase("phs", $wrkdir, $reads, $prictg, "$prjDir/4-phase/rd2ctg.paf");
    my $job_filter = $self->job_filter_inconsistent("phs", $wrkdir, $rd2rd, $rd2rd_flt);
    
    return $self->newjob(
        name => "phs_step",
        ifiles => [$reads, $prictg, $altctg],
        ofiles => [$rd2rd_flt],
        mfiles => [$rd2ctg],
        jobs => [$job_map, $job_phase, $job_filter],
        msg => "phasing reads",
    );

}


sub run_phase_with_contig($) {
    my ($self) = @_;
    
    my $name = "phs";
    my $wrkdir = $self->get_work_folder("4-phase");

    mkdir $wrkdir;
    
    my $method = $self->get_config("phase_method");

    my $job = undef;
    if ($method == 1) {
        $job = $self->job_phase_method1($name, $wrkdir)
    } else {
        $job = $self->job_phase_method0($name, $wrkdir);
    }
    
    $self->run_jobs($job);
}


sub job_phase_method1($$$) {
    my ($self, $name, $wrkdir) = @_;


    mkdir $wrkdir;

    my $isGz = $self->get_config("COMPRESS");
    my $prj_dir = $self->get_project_folder();

    my $corrReads = $isGz ? "$prj_dir/1-correct/corrected_reads.fasta.gz" : "$prj_dir/1-correct/corrected_reads.fasta";
    my $prepReads = $isGz ? "$prj_dir/0-prepare/prepared_reads.fasta.gz" : "$prj_dir/0-prepare/prepared_reads.fasta";
    #my $reads = $prepReads; # TODO corrReads
    my $reads = $corrReads;
    my $contigs = "$prj_dir/3-assemble/primary.fasta";
    my $ol_r2r = "$prj_dir/2-align/overlaps.txt";
    my $ol_r2r_s = "$wrkdir/ol_c2c_s.paf";

    my $job_call = $self->job_calling_with_deepvariant($name, $wrkdir, $reads, $contigs,  "-x ccs -g 140m");
    my $job_phase = $self->job_phasing_with_whatshap($name, $wrkdir, "$wrkdir/rd2ctg.bam", $contigs, "$wrkdir/rd2ctg.filtered.vcf");
    my $job_filter = $self->job_filter_inconsistent("phs", $wrkdir, $ol_r2r, "$wrkdir/phased", $ol_r2r_s);

    return $self->newjob(
        name => "phs_step",
        ifiles => [$reads, $contigs],
        ofiles => [$ol_r2r_s],
        mfiles => [],
        jobs => [$job_call, $job_phase, $job_filter],
        msg => "phasing reads",
    );

}

sub job_calling_with_deepvariant($$$$$) {
    my ($self, $name, $wrkdir, $reads, $contigs, $vcf) = @_;

    mkdir $wrkdir;
    my $prj_name = $self->get_config("project"); 
    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");
    my $reads_fastq = "$wrkdir/reads.fastq";
    my $rd2ctg = "$wrkdir/rd2ctg";
    my $map_options = "-k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k --secondary=no";
    my $s = "\"\@RG\\tSM:$prj_name\\tID:$prj_name\"";

    my $job_al = $self->newjob(
        name => "${name}_al",
        ifiles => [$reads, $contigs],
        ofiles => ["$rd2ctg.bam"],
        gfiles => ["$rd2ctg.bam"],
        mfiles => [],
        cmds => ["$bin_path/fxtools.py fx_fa2fq $reads $reads_fastq",
                 "minimap2 $map_options -a --eqx -R \"\@RG\\tSM:$prj_name\\tID:$prj_name\"  -t $threads $contigs $reads_fastq | samtools sort -\@$threads --output-fmt BAM -o $rd2ctg.bam",
                 "samtools index -\@$threads $rd2ctg.bam"],
        msg => "mapping reads to contigs",
    );

    my $prj_dir = $self->get_project_folder();
    my $rd2ctg_rel = "4-phase/rd2ctg";
    my $contigs_rel = "3-assemble/primary.fasta";
    my $job_call = $self->newjob(
        name => "${name}_call",
        ifiles => ["$rd2ctg.bam", $contigs],
        ofiles => ["$rd2ctg.filtered.vcf"],
        gfiles => ["$rd2ctg.filtered.vcf"],
        mfiles => [],
        cmds => ["samtools faidx $contigs",
                 "docker run -v \"$prj_dir\":\"/input\" -v \"$prj_dir\":\"/output\" google/deepvariant:0.8.0  /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=/input/$contigs_rel  --reads=/input/$rd2ctg_rel.bam --output_vcf=/output/$rd2ctg_rel.vcf.gz  --num_shards=$threads",
                 #"/public/software/singularity/bin/singularity run --cleanenv -B \"$prj_dir\":\"/input\" -B \"$prj_dir\":\"/output\" ~/test/bin/deepvariant.img  /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=/input/$contigs_rel  --reads=/input/$rd2ctg_rel.bam --output_vcf=/output/$rd2ctg_rel.vcf.gz  --num_shards=$threads",
                 "bgzip -cd $rd2ctg.vcf.gz > $rd2ctg.vcf",
                 "grep -E '^#|0/0|1/1|0/1|1/0|0/2|2/0' $rd2ctg.vcf > $rd2ctg.filtered.vcf"],
        msg => "calling variant",
    );
    
    return $self->newjob(
        name => "${name}_job",
        ifiles => [$reads, $contigs],
        ofiles => ["$rd2ctg.filtered.vcf"],
        mfiles => [],
        jobs => [$job_al, $job_call],
        msg => "phasing reads",
    );
}

sub job_phasing_with_whatshap($$$$$) {
    my ($self, $name, $wrkdir, $rd2ctg, $contigs, $vcf) = @_;

    mkdir $wrkdir;
    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");

    my $phased = "$wrkdir/phased";

    #my $haplotag = "$wrkdir/haplotag.bam";
    
    # my $job_phase = $self->newjob(
    #     name => "${name}_phase",
    #     ifiles => ["$rd2ctg", $contigs],
    #     ofiles => ["$phased.vcf.gz"],
    #     gfiles => ["$phased.vcf.gz"],
    #     mfiles => [],
    #     cmds => ["whatshap phase --reference $contigs $vcf $rd2ctg -o $phased.vcf",
    #              "bgzip -c $phased.vcf > $phased.vcf.gz",
    #              "tabix -p vcf $phased.vcf.gz"],
    #     msg => "phase",
    # );
    
    # my $job_hap = $self->newjob(
    #     name => "${name}_hap",
    #     ifiles => ["$rd2ctg", $contigs],
    #     ofiles => [$phased],
    #     gfiles => [$phased],
    #     mfiles => [],
    #     cmds => ["whatshap haplotag --reference $contigs $phased.vcf.gz $rd2ctg -o $haplotag",
    #              "samtools view $haplotag | $bin_path/fxtools.py fx_get_phased_reads - $phased"],

    #     msg => "phase",
    # );
    
    my $job_phase = $self->newjob(
        name => "${name}_phase",
        ifiles => ["$rd2ctg", $contigs],
        ofiles => [$phased],
        gfiles => [$phased],
        mfiles => [],
        cmds => ["$bin_path/phase_using_whatshap.py $contigs $vcf $rd2ctg $phased --threads $threads --wrkdir $wrkdir/wrkdir"],
        msg => "phase",
    );

    return $self->newjob(
        name => "${name}_job",
        ifiles => [$rd2ctg, $contigs, $vcf],
        ofiles => [$phased],
        mfiles => [],
        jobs => [$job_phase],
        msg => "phasing reads",
    );
}

sub run_assemble2($$$$$) {
    my ($self) = @_;

    my $name = "asm2";
    my $wrkdir = $self->get_work_folder("5-assemble");

    my $prjDir = $self->get_project_folder();
    my $isGz = $self->get_config("COMPRESS");

    my $reads = $isGz ? "$prjDir/1-correct/corrected_reads.fasta.gz" : "$prjDir/1-correct/corrected_reads.fasta";
    my $overlaps = "$prjDir/4-phase/filtered_overlaps.paf";


    my $rd2ctg_option = " --rd2ctg_fname $wrkdir/../4-phase/rd2ctg.paf";

    $self->runAssemble($name, $wrkdir, $reads, $overlaps, [$self->get_config("ASM2_FILTER_OPTIONS") . $rd2ctg_option, $self->get_config("ASM2_ASSEMBLE_OPTIONS")]);
    
}

sub run_polish($) {
    my ($self) = @_;
    my $name = "pol";
    my $prjdir = $self->get_project_folder();
    my $wrkdir = $self->get_work_folder("6-polish");
    my $wrkdir_asm = $self->get_work_folder("5-assemble");
    
    my $ctg_pri = "$wrkdir_asm/primary.fasta";
    my $ctg_alt = "$wrkdir_asm/alternate.fasta";
    my $ctg_all = "$wrkdir/ctgall.fasta";

    my $tile_pri = "$wrkdir_asm/primary_tiles";
    my $tile_alt = "$wrkdir_asm/alternate_tiles";
    my $tile_all = "$wrkdir/all_tiles";

    my $readinfo = "$wrkdir/../4-phase/readinfos";

    my $rd2ctg = "$wrkdir/rd2ctg.paf";
    my $rd2ctg_flt = "$wrkdir/rd2ctg_flt.paf";
    my $rd2ctg_pri = "$wrkdir/rd2ctg_pri.paf";
    my $rd2ctg_alt = "$wrkdir/rd2ctg_alt.paf";


    my $pol_pri = "$wrkdir/primary.fasta";
    my $pol_alt = "$wrkdir/alternate.fasta";
    
    my $read_crr = "$prjdir/1-correct/corrected_reads.fasta";
    my $read_prp = "$prjdir/0-prepare/prepared_reads.fasta";
    #my $reads = $prepReads; # TODO corrReads
    my $reads = $self->get_config("polish_use_reads") + 0 == 1 ? $read_crr : $read_prp;

    my $map_options = $self->get_config("polish_map_options");
    my $filter_options = $self->get_config("polish_filter_options");
    my $cns_options = $self->get_config("polish_cns_options");

    print($map_options);
    if ($map_options=~/\-[a-zAz]*a/) {
        $rd2ctg = "$wrkdir/rd2ctg.sam";
        $rd2ctg_flt = "$wrkdir/rd2ctg_flt.sam";
        $rd2ctg_pri = "$wrkdir/rd2ctg_pri.sam";
        $rd2ctg_alt = "$wrkdir/rd2ctg_alt.sam";
    }

    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");
    
    my $filter_cmd = "";
    if ($filter_options ne "") {
        $filter_cmd = " | $bin_path/fsa_ol_refine - - --itype paf --otype paf --thread_size 4 $filter_options ";
    }

    mkdir $wrkdir;

    my $job_map = $self->newjob(
        name => "${name}_map",
        ifiles => [$ctg_pri, $ctg_alt, $reads],
        ofiles => [$rd2ctg],
        gfiles => [$ctg_all, $rd2ctg],
        mfiles => [$ctg_all],
        cmds => ["cat $ctg_pri $ctg_alt > $ctg_all",
                 "minimap2  -t $threads $map_options $ctg_all $reads $filter_cmd > $rd2ctg"],
        msg => "mapping reads to contigs",
    );

    my $job_flt = $self->newjob(
        name => "${name}_filter",
        ifiles => [$rd2ctg, $tile_pri, $tile_alt],
        ofiles => [$rd2ctg_pri, $rd2ctg_alt],
        gfiles => [$rd2ctg_flt, $rd2ctg_pri, $rd2ctg_alt],
        mfiles => [$rd2ctg_flt],
        cmds => ["$bin_path/fxtools.py fx_purge_overlaps $rd2ctg $readinfo --tile $tile_pri --tile $tile_alt > $rd2ctg_flt",
                 "$bin_path/fxtools.py fx_split_mappings $rd2ctg_flt $ctg_pri,$ctg_alt $rd2ctg_pri,$rd2ctg_alt"],
        msg => "filter overlaps in which the reads are different haplotype",
    );
    
    my $job_pol = $self->newjob(
        name => "${name}_polish",
        ifiles => [$ctg_pri, $ctg_alt, $reads, $rd2ctg_pri, $rd2ctg_alt],
        ofiles => [$pol_pri, $pol_alt],
        mfiles => [],
        pjobs => [$self->newjob(
                    name => "${name}_polish_prj",
                    ifiles => [$ctg_pri, $reads, $rd2ctg_pri],
                    ofiles => [$pol_pri],
                    gfiles => [$pol_pri],
                    mfiles => [],
                    cmds => ["if [ -s $ctg_pri ]; then racon $cns_options -t $threads $reads $rd2ctg_pri $ctg_pri > $pol_pri; else touch $pol_pri; fi"],
                    msg => "polishing primary contigs",
                  ), 
                  $self->newjob(
                    name => "${name}_polish_alt",
                    ifiles => [$ctg_alt, $reads, $rd2ctg_alt],
                    ofiles => [$pol_alt],
                    gfiles => [$pol_alt],
                    mfiles => [],
                    cmds => ["if [ -s $ctg_alt ]; then racon $cns_options -t $threads $reads $rd2ctg_alt $ctg_alt > $pol_alt; else touch $pol_alt; fi"],
                    msg => "polishing alternate contigs",

                  )
        ],
        msg => "polishing contigs",
    );

    $self->run_jobs($self->newjob(
        name => "${name}_job",
        ifiles => [$ctg_pri, $ctg_alt, $reads],
        ofiles => [$pol_pri, $pol_alt],
        mfiles => [$rd2ctg_flt, $rd2ctg],
        jobs => [$job_map, $job_flt, $job_pol],
        msg => "polishing contigs",
    ));

}

sub run_pbgcpp($) {
    my ($self) = @_;
    my $name = "gcpp";
    my $prjdir = $self->get_project_folder();
    my $wrkdir = $self->get_work_folder("7-pbgcpp");
    my $wrkdir_asm = $self->get_work_folder("5-assemble");
    my $wrkdir_pol = $self->get_work_folder("6-polish");

    my $ctg_pri = "$wrkdir_pol/polished.fasta";
    my $ctg_alt = "$wrkdir_pol/polished1.fasta";
    my $ctg_all = "$wrkdir/ctgall.fasta";

    my $tile_all = "$wrkdir_pol/all_tiles";
    my $id2name = "$wrkdir/../0-prepare/id2name.gz";

    my $readinfo = "$wrkdir/../4-phase/readinfos";

    my $rd2ctg = "$wrkdir/rd2ctg.bam";
    my $rd2ctg_flt = "$wrkdir/rd2ctg_flt.bam";

    my $reads = $self->get_config("gcpp_bam_fofn");

    my $map_options = $self->get_config("gcpp_map_options");
    my $cns_options = $self->get_config("gcpp_cns_options");

    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");
    
    my $pol_pri = "$wrkdir/polished.fasta";
    my $pol_alt = "$wrkdir/polished1.fasta";

    mkdir $wrkdir;

    my $job_map = $self->newjob(
        name => "${name}_map",
        ifiles => [$ctg_pri, $ctg_alt, $reads],
        ofiles => [$rd2ctg],
        gfiles => [$ctg_all, $rd2ctg],
        mfiles => [$ctg_all],
        cmds => ["cat $ctg_pri $ctg_alt > $ctg_all",
                 "pbmm2 align -j $threads --sort $map_options $ctg_all $reads $rd2ctg"],
        msg => "mapping reads to contigs",
    );

    my $job_flt = undef;
    if ($self->get_config("gcpp_filter") + 0 == 1) {
        $job_flt = $self->newjob(
            name => "${name}_filter",
            prefunc => sub($) {
                if (-l $rd2ctg_flt) {
                    unlink $rd2ctg_flt;
                } 
                if (-l "$rd2ctg_flt.bai") {
                    unlink "$rd2ctg_flt.bai";
                }
            },
            ifiles => [$rd2ctg, $tile_all],
            ofiles => [$rd2ctg_flt],
            gfiles => [$rd2ctg_flt],
            mfiles => [$tile_all],
            cmds => ["$bin_path/fxtools.py fx_purge_rd2ctg $rd2ctg $rd2ctg_flt $tile_all $readinfo $id2name --threads $threads",
                     "samtools index -@ $threads $rd2ctg_flt"],
            msg => "filter overlaps in which the reads are different haplotype",
        );
    } else {
        $job_flt = $self->newjob(
            name => "${name}_filter",
            ifiles => [$rd2ctg],
            ofiles => [$rd2ctg_flt],
            gfiles => [$rd2ctg_flt, "$rd2ctg_flt.bai"],
            mfiles => [$tile_all],
            funcs => [ sub ($$) {
                `ln -s -f $rd2ctg $rd2ctg_flt`;
                `ln -s -f $rd2ctg.bai $rd2ctg_flt.bai`;
            }],
            msg => "filter overlaps in which the reads are different haplotype",
        );
    }
    
    my $job_pol = $self->newjob(
        name => "${name}_polish",
        ifiles => [$ctg_pri, $ctg_alt, $reads, $rd2ctg_flt],
        ofiles => [$pol_pri, $pol_alt],
        mfiles => [],
        pjobs => [$self->newjob(
                    name => "${name}_polish_prj",
                    ifiles => [$ctg_pri, $reads, $rd2ctg_flt],
                    ofiles => [$pol_pri],
                    gfiles => [$pol_pri. "$ctg_pri.fai"],
                    mfiles => [],
                    cmds => ["gcpp $cns_options -j $threads -r $ctg_pri $rd2ctg_flt -o $pol_pri"],
                    msg => "polishing primary contigs",
                  ), 
                  $self->newjob(
                    name => "${name}_polish_alt",
                    ifiles => [$ctg_alt, $reads, $rd2ctg_flt],
                    ofiles => [$pol_alt],
                    gfiles => [$pol_alt. "$ctg_alt.fai"],
                    mfiles => [],
                    cmds => ["gcpp $cns_options -j $threads -r $ctg_alt $rd2ctg_flt -o $pol_alt"],
                    msg => "polishing alternate contigs",
                  )
        ],
        msg => "polishing contigs",
    );

    $self->run_jobs($self->newjob(
        name => "${name}_job",
        ifiles => [$ctg_pri, $ctg_alt, $reads],
        ofiles => [$pol_pri, $pol_alt],
        mfiles => [$rd2ctg_flt, $rd2ctg],
        jobs => [$job_map, $job_flt, $job_pol],
        msg => "polishing contigs",
    ));

}

sub stat_read_n50($$$) {
    my ($self, $seq, $msg) = @_;

    my $bin_path = $self->get_env("BinPath");

    Plgd::Logger::info("N50 of $msg: $seq");
    my $cmd = "$bin_path/fsa_rd_tools n50  $seq";
    print $cmd;
    system($cmd);
}


package main;

my @defaultConfig = (
    ["project", "", 1, "project name"],
    ["reads", "", 1, "reads path"],
    ["genome_size", "", 1, "genome size"],
    ["threads", "4", 0],
    ["memory", "0", 0],
    ["cleanup", "0", 0],
    ["compress", "0", 0],
    ["grid", "auto:0", 0],
    ["prep_min_length", "3000", 0],
    ["prep_output_coverage", "80", 0],
    ["corr_iterate_number", "1", 0],
    ["corr_block_size", "4000000000", 0],
    ["corr_correct_options", "", 0],
    ["corr_filter_options", "--filter0=l=2000:al=2000:alr=0.5", 0],
    ["corr_rd2rd_options", "-x ava-pb", 0],
    ["corr_output_coverage", "80", 0],
    ["align_block_size", "4000000000", 0],
    ["align_rd2rd_options", "-X -g3000 -w30 -k19 -m100 -r500", 0],
    ["align_filter_options", "--filter0=l=2000:aal=4000:oh=3000:ohr=0.3 --task=extend --filter1=oh=100:ohr=0.01", 0],
    ["asm1_assemble_options", "", 0],
    ["phase_method", "", 0],
    ["phase_rd2ctg_options", "-x map-pb -c -p 0.5 -r 1000", 0],
    ["phase_use_reads", "1", 0],
    ["phase_phase_options", "", 0],
    ["asm2_assemble_options", "", 0],
    ["polish_map_options", "-x map-pb", 0], #  --secondary=no
    ["polish_filter_options", "--filter0 oh=1000:ohr=0.1", 0],
    ["polish_cns_options", "", 0],
#    ["gcpp_bam_fofn", "", 0], 
#    ["gcpp_map_options", "", 0],
#    ["gcpp_cns_options", "--algorithm=arrow -x 5 -X 120 -q 0", 0],
#    ["gcpp_filter", "1", 0],
    
);


my $pipeline = FsaPipeline->new(\@defaultConfig);


sub cmd_correct($) {
    my ($fname) = @_;

    $pipeline->initialize($fname);

    $pipeline->run_prepare();
    $pipeline->run_correct();
}

sub cmd_assemble($) {
    my ($fname) = @_;

    cmd_correct($fname);

    $pipeline->initialize($fname);
    $pipeline->run_align();
    $pipeline->run_assemble1();
    
}

sub cmd_unzip($) {
    my ($fname) = @_;

    cmd_assemble($fname);
    $pipeline->initialize($fname);

    $pipeline->run_phase_with_contig();
    $pipeline->run_assemble2();
    $pipeline->run_polish();
}

sub cmd_gcpp($) {
    my ($fname) = @_;

    cmd_unzip($fname);
    
    $pipeline->initialize($fname);
    $pipeline->run_pbgcpp();
}

sub cmd_config($) {
    my ($fname) = @_;

    open(F, "> $fname") or die; 
    foreach my $item (@defaultConfig) {
        print F "$item->[0]=$item->[1]\n";
    }

    close(F);

}


sub usage() {
    print "Usage: fsa.pl correct|assemble|bridge|config cfg_fname\n".
          "    correct:     correct rawreads\n" .
          "    assemble:    generate contigs\n" .
          "    unzip:      \n" .
          "    gcpp:        polish assembly with gcpp" .
          "    config:      generate default config file\n" 
}

sub main() {
    if (scalar @ARGV >= 2) {
        my $cmd = @ARGV[0];
        my $cfgfname = @ARGV[1];

        if ($cmd eq "correct") {
            cmd_correct($cfgfname);
        } elsif ($cmd eq "assemble") {
            cmd_assemble($cfgfname);
        } elsif ($cmd eq "unzip") {
            cmd_unzip($cfgfname);
        } elsif ($cmd eq "gcpp") {
            cmd_gcpp($cfgfname);
        } elsif ($cmd eq "test") {
            cmd_test($cfgfname);
        } elsif ($cmd eq "config") {
            cmd_config($cfgfname);
        } else {
            usage();
        }
    } else {
        usage();
    }
}


$SIG{TERM}=$SIG{INT}=\& catchException;
sub catchException { 
    Plgd::Logger::info("Catch an Exception, and do cleanup");
    $pipeline->stop_running();
    exit -1; 
} 

#eval {
    main();
#};

if ($@) {
    catchException();
}

END {
    $pipeline->stop_running();
}
