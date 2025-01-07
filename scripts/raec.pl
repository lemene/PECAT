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

    my $threads = $self->get_config("threads");
    my $bin_path = $self->get_env("BinPath");

    return $self->newjob(
        name => "${name}_correct",
        ifiles => [$rreads, $unmapped, $rd_2_ref],
        ofiles => [$creads],
        gfiles => [$creads],
        mfiles => [],
        cmds => ["$bin_path/fsa_rd_correct $unmapped $rreads $creads --output_directory=$wrkdir --thread_size=$threads " . 
                    "--infos_fname $creads.infos $options --rd_2_ref $rd_2_ref"],
        msg => "correcting reads0, $name"
    );
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
    my $creads = "$wrkdir/corrected.fasta";
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
