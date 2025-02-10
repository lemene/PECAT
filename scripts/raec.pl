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

package RaecPipeline;
our @ISA = qw(FsaPipeline);  

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


sub job_get_unmaped_reads($$$$$$) {
    my ($self, $name, $reads, $rd_2_ref, $unmapped, $wrkdir) = @_;

    my $threads = $self->get_config("threads");
    my $bin_path = $self->get_env("BinPath");

    my $mapped_names = "$wrkdir/mapped.txt";
    return $self->newjob(
        name => "${name}_get_unmmaped",
        ifiles => [$rd_2_ref, $reads],
        ofiles => [$unmapped],
        gfiles => [$unmapped],
        mfiles => [$mapped_names],
        cmds => ["$bin_path/fsa_ol_tools location --thread_size $threads $rd_2_ref  $mapped_names",
                 "$bin_path/fsa_rd_tools sub $reads $unmapped --names_fname $mapped_names --complement"],
        msg => "get unmapped reads, ${name}",
    );
}

sub getjob_align_unmapped_reads($$$$$$) {
    my ($self, $name, $unmapped, $overlaps, $options, $wrkdir) = @_;

    my $threads = $self->get_config("threads");
    my $bin_path = $self->get_env("BinPath");

    # unmapped reads vs unmapped reads
    
    my $map_options = $options->[0];
    my $flt_options = $options->[1];

    return $self->newjob(
        name => "${name}_align_unmmaped",
        ifiles => [$unmapped],
        ofiles => [$overlaps],
        gfiles => [$overlaps],
        mfiles => [],
        cmds => ["minimap2 $map_options -t $threads $unmapped $unmapped | $bin_path/fsa_ol_refine - - --itype paf --otype paf  --read_fname $unmapped --thread_size $threads $flt_options > $overlaps"],
        msg => "aligning unmapped reads, ${name}",
    );
}


sub getjob_correct_all_reads() {
    my ($self, $name, $rreads, $rd_2_ref, $unmapped, $creads, $options, $wrkdir) = @_;

    #my $isGz = ($corrected =~ /\.gz$/);
    

    my $readname = "$wrkdir/readname";

    my $block_size = $self->get_config("corr_block_size");
    my $threads = $self->get_config("threads");
    my $bin_path = $self->get_env("BinPath");
    my $block_info = "$wrkdir/block_info";

    my $job_split = $self->newjob(
        name => "${name}_split",
        ifiles => [$rd_2_ref],
        ofiles => [$block_info],
        gfiles => [$block_info, "$readname.*"],
        mfiles => [],
        cmds => ["$bin_path/fsa_sam_tools group $rd_2_ref $readname.core.{}  --block_size $block_size  --thread_size $threads",
                 "ls $readname.core.* > $block_info"],
        msg => "spliting read names, $name",
    );

    
    my $job_corr = $self->newjob(
        prefunc => sub($) {
            my ($job) = @_;
            my $size = `wc -l $block_info`;
            for (my $i=0; $i < $size; $i=$i+1) {

                my $corr_sub = "$creads.$i";
                
                my $job_sub = $self->newjob(
                    name => "${name}_correct_$i",
                    ifiles => [$rreads, $block_info, $rd_2_ref, $unmapped],
                    ofiles => [$corr_sub],
                    gfiles => [$corr_sub],
                    mfiles => ["$readname.core.$i"],
                    cmds => ["$bin_path/fsa_rd_correct $unmapped $rreads $creads.$i --output_directory=$wrkdir --thread_size=$threads " . 
                                " --read_name_fname=$readname.core.$i --infos_fname $creads.$i.infos $options --rd_2_ref $rd_2_ref"],
                    msg => "correcting reads $i, $name"
                );
                push @{$job->{ofiles}}, $corr_sub;
                push @{$job->{pjobs}}, $job_sub;
            }

        },
        name => "${name}_correct_all",
        ifiles => [$rreads, $block_info, $rd_2_ref, $unmapped],
        ofiles => [],                   # prefunc
        mfiles => [],
        pjobs => [],                    # prefunc
        msg => "correcting rawreads, $name",
    );



    my $job_cat = $self->newjob(
        prefunc => sub($) {
            my ($job) = @_;
            my $size = `wc -l $block_info`;

            my @curr_sub = ();
            for (my $i=0; $i < $size; $i=$i+1) {
                $curr_sub[$i] = "$creads.$i";
            }

            push @{$job->{ifiles}}, @curr_sub;
            push @{$job->{cmds}}, "cat @curr_sub > $creads && rm @curr_sub";

        },
        name => "${name}_cat",
        ifiles => [],      # prefunc
        ofiles => [$creads], 
        gfiles => [$creads], 
        mfiles => [],
        cmds => [],                     # prefunc
        threads => 1,
        msg => "cat corrected reads, $name",

    );
    
    return $self->newjob(
        name => "${name}_correct",
        ifiles => [$rreads],
        ofiles => [$creads], # prefunc
        mfiles => ["$readname.core.*"],
        jobs => [$job_split, $job_corr, $job_cat],
        msg => "correcting rawreads, $name");

    # return $self->newjob(
    #     name => "${name}_correct",
    #     ifiles => [$rreads, $unmapped, $rd_2_ref],
    #     ofiles => [$creads],
    #     gfiles => [$creads],
    #     mfiles => [],
    #     cmds => ["$bin_path/fsa_rd_correct $unmapped $rreads $creads --output_directory=$wrkdir --thread_size=$threads " . 
    #                 "--infos_fname $creads.infos $options --rd_2_ref $rd_2_ref"],
    #     msg => "correcting reads0, $name"
    # );
}

sub getjob_correct_with_reference($$$$$$) {
    my ($self, $name, $rreads, $ref, $creads, $options, $wrkdir) = @_;
    
    my $threads = $self->get_config("threads");
    my $bin_path = $self->get_env("BinPath");

    my $opts_rd_2_ref = $options->[0];
    my $opts_rd_2_rd = $options->[1];
    my $opts_rd_2_rd_flt = $options->[2];
    my $opts_correct = $options->[3];

    # map reads to reference
    my $rd_2_ref = "$wrkdir/rd_2_ref.bam";
    my $job_rd_2_ref = $self->getjob_map_read_to_ref_sam($name, $rreads, $ref, 
        $rd_2_ref, $opts_rd_2_ref, $wrkdir);

    # 
    my $unmapped = "$wrkdir/unmapped.fasta";
    my $job_get_unmapped = $self->job_get_unmaped_reads($name, $rreads, $rd_2_ref, $unmapped, $wrkdir);
    
    my $unmapped_overlaps = "$wrkdir/unmapped.paf";
    my $job_align_ummapped = $self->getjob_align_unmapped_reads($name, $unmapped, $unmapped_overlaps, 
        [$opts_rd_2_rd, $opts_rd_2_rd_flt], $wrkdir);

    # correct
    my $creads = "$wrkdir/corrected_reads.fasta";
    my $job_correct_reads = $self->getjob_correct_all_reads($name, $rreads, $rd_2_ref, $unmapped_overlaps, 
        $creads, $opts_correct, $wrkdir);

    return $self->newjob(
        name => "${name}_job",
        ifiles => [$rreads, $ref],
        ofiles => [], 
        mfiles => [],
        jobs => [$job_rd_2_ref, $job_get_unmapped, $job_align_ummapped, $job_correct_reads],
        msg => "correcting reads assisting with reference, $name"
    );
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
        
    $self->runAssemble($name, $wrkdir, $reads, $overlaps, $self->get_config("ASM1_ASSEMBLE_OPTIONS"));
}


sub run_correct($) {
    my ($self) = @_;

    my $name = "crr";
    my $wrkdir = $self->get_work_folder("1-correct");
    mkdir $wrkdir;

    my $is_gz = $self->get_config("compress");
    my $wrkdir_prp = $self->get_work_folder("0-prepare");

    my $rreads = $is_gz ? "$wrkdir_prp/prepared_reads.fasta.gz" : "$wrkdir_prp/prepared_reads.fasta";
    my $creads = $is_gz ? "$wrkdir/corrected_reads.fasta.gz" : "$wrkdir/corrected_reads.fasta";
    my $ref = Cwd::abs_path($self->get_config("reference"));

    my $opts_rd_2_ref = $self->get_config("corr_rd2ref_options");
    my $opts_rd_2_rd = $self->get_config("corr_rd2rd_options");
    my $opts_rd_2_rd_flt = $self->get_config("corr_filter_options");
    my $opts_correct = $self->get_config("corr_correct_options");

    $self->run_jobs($self->getjob_correct_with_reference($name, $rreads, $ref, $creads, 
        [$opts_rd_2_ref, $opts_rd_2_rd, $opts_rd_2_rd_flt, $opts_correct], $wrkdir));
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
    ["reference", "", 0, "reference path"],
    ["threads", "4", 0],
#    ["memory", "0", 0],
    ["cleanup", "1", 0],
#    ["compress", "0", 0],
    ["grid", "local", 0],

    ["corr_iterate_number", "1", 0],
    ["corr_block_size", "4000000000", 0],
    ["corr_correct_options", "", 0],
    ["corr_filter_options", "--filter0=l=5000:al=2500:alr=0.5:aal=5000:oh=3000:ohr=0.3", 0],
    ["corr_rd2rd_options", "-x ava-ont", 0],
);


my $pipeline = RaecPipeline->new(\@defaultConfig);


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


sub cmd_config($) {
    my ($fname) = @_;

    open(F, "> $fname") or die; 
    foreach my $item (@defaultConfig) {
        print F "$item->[0]=$item->[1]\n";
    }

    close(F);
}


sub usage() {
    print "Usage: raec.pl correct|config cfg_fname\n".
          "    correct:     correct rawreads\n" .
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
    #$pipeline->stop_running();
    exit -1; 
} 

#eval {
    main();
#};

if ($@) {
    catchException();
}

END {
    #$pipeline->stop_running();
}
