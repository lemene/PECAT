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
        return $self->get_work_folder("1-correct");
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
        
    $self->runAssemble($name, $wrkdir, $reads, $overlaps, $self->get_config("ASM1_ASSEMBLE_OPTIONS"));
}


sub jobCorrect($$$$$) {
    my ($self, $name, $rawreads, $corrected, $workDir) = @_;
    
    my $wrkdir = $workDir;
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
        ifiles => [$rawreads, $rd2rd],
        ofiles => [$blockInfo],
        gfiles => [$blockInfo, "$readName.*"],
        mfiles => ["$wrkdir/cc*.fasta.paf"],
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
                    #ifiles => [$rawreads, "$rd2rd.sub.$i.paf", "$readName.core.$i"],
                    ifiles => [$rawreads, $blockInfo],
                    ofiles => [$corrSub],
                    gfiles => [$corrSub],
                    mfiles => ["$rd2rd.sub.$i.paf", "$readName.core.$i"],
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
        mfiles => ["$rd2rd.sub.*.paf", "$rd2rd", "$readName.core.*"],
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


sub newjob_filter_inconsistent($$$$$$) {
    my ($self, $name, $wrkdir, $overlaps, $filtered, $options) = @_;
 
    
    my $binPath = $self->get_env("BinPath"); 
    my $threads = $self->get_config("threads");
    my $consistent = "$wrkdir/consistent";
    my $inconsistent = "$wrkdir/inconsistent";
    my $crr_reads = "$wrkdir/../1-correct/corrected_reads.fasta";

    my $job = $self->newjob(
        name => "${name}_filter",
        ifiles => [$overlaps, $inconsistent],
        ofiles => [$filtered],
        gfiles => [$filtered],
        mfiles => [],
        cmds => ["$binPath/fsa_ol_tools filter $overlaps $filtered --thread_size=$threads --inconsistent=$inconsistent $options"],
        msg => "filtering inconsistent overlaps, ${name}",
    );

    return $job;
}


sub job_phase_method0($$$$$) {
    my ($self, $name, $wrkdir, $rd2rd , $rd2rd_flt) = @_;

    mkdir $wrkdir;

    my $isGz = $self->get_config("COMPRESS");
    my $prjDir = $self->get_project_folder();

    my $corrReads = $isGz ? "$prjDir/1-correct/corrected_reads.fasta.gz" : "$prjDir/1-correct/corrected_reads.fasta";
    my $prepReads = $isGz ? "$prjDir/0-prepare/prepared_reads.fasta.gz" : "$prjDir/0-prepare/prepared_reads.fasta";

    my $reads = $self->get_config("phase_use_reads") + 0 == 1 ? $corrReads : $prepReads;
    my $prictg = "$prjDir/3-assemble/primary.fasta";
    my $altctg = "$prjDir/3-assemble/alternate.fasta";
    my $rd2ctg = "$wrkdir/rd2ctg.paf";

    
    #my $rd2rd = "$prjDir/2-align/overlaps.txt";
    #my $rd2rd_flt = "$wrkdir/filtered_overlaps.paf";

    my $job_map = $self->job_map_reads_to_contigs("phs", $wrkdir, $reads, [$prictg, $altctg], $self->get_config("PHASE_RD2CTG_OPTIONS"));
    my $job_phase = $self->job_phase("phs", $wrkdir, $reads, $prictg, "$wrkdir/rd2ctg.paf", "");
    my $job_filter = $self->newjob_filter_inconsistent("phs", $wrkdir, $rd2rd, $rd2rd_flt);

    return $self->newjob(
        name => "phs_step",
        ifiles => [$reads, $prictg, $altctg],
        ofiles => [$rd2rd_flt],
        mfiles => [$rd2ctg, "$wrkdir/consistent", "$wrkdir/variants"],
        jobs => [$job_map, $job_phase, $job_filter],
        msg => "phasing reads",
    );
}

## about Hi-C
# {
#     my @hic_reads = split(";", $self->get_config("hic_reads"));
#     my $hic1_2_ctg = "$wrkdir/hic1_2_ctg.paf";
#     my $hic2_2_ctg = "$wrkdir/hic2_2_ctg.paf";
#     my $snp_in_hic = "$wrkdir/hic_infos";

#     my $job_map_hic = $self->job_map_hic_reads_to_contigs("phs", $wrkdir, \@hic_reads, [$prictg, $altctg], "");
#     my $job_snp_in_hic = $self->job_identify_snps_in_hic("phs", $wrkdir, \@hic_reads);

#     my $job_map_all = $self->newjob(
#         name => "phs_map_step",
#         ifiles => [], # [$reads, @hic_reads, $prictg, $altctg],
#         ofiles => [$rd2ctg, $hic1_2_ctg, $hic2_2_ctg], #TODO 有两份独立依赖关系
#         mfiles => [],
#         pjobs => [$job_map, @$job_map_hic],
#         msg => "phasing reads and hic reads",
#     );

#     return $self->newjob(
#         name => "phs_step1",
#         ifiles => [], # [$reads, $prictg, $altctg],
#         ofiles => [$rd2rd_flt, $snp_in_hic], #TODO 有两份独立依赖关系
#         mfiles => [$rd2ctg, $hic1_2_ctg, $hic2_2_ctg],
#         jobs => [$job_map_all, $job_phase, $job_filter, $job_snp_in_hic],
#         msg => "phasing reads",
#     );
# }

sub run_phase_with_contig($) {
    my ($self) = @_;
    
    my $name = "phs";
    my $wrkdir = $self->get_work_folder("4-phase");
    my $prjdir = $self->get_project_folder();

    mkdir $wrkdir;
    
    my $method = $self->get_config("phase_method") + 0;

    ### input
    my $crr_reads = "$prjdir/1-correct/corrected_reads.fasta";
    my $raw_reads = "$prjdir/0-prepare/prepared_reads.fasta";
    my $prictg = "$prjdir/3-assemble/primary.fasta";
    my $rd2rd = "$prjdir/2-align/overlaps.txt";

    ### output
    my $rd2rd_flt = "$wrkdir/filtered_overlaps.paf";

    my $job = undef;
    if ($method == 0) {

        ### output
        my $inconsistent = "$wrkdir/fsa/inconsistent";
        my $readinfos = "$wrkdir/fsa/readinfos";

        my $job_inconsistent = $self->newjob_identify_inconsistent_overlaps_fsa($name, $wrkdir . "/fsa");

        my $filter_options = $self->get_config("phase_filter_options");
        if ($self->get_config("phase_use_reads") + 0 == 0) {
            $filter_options = $filter_options . " --range $crr_reads";
        }
        my $job_filter = $self->newjob_filter_inconsistent($name."_fsa", $wrkdir. "/fsa", $rd2rd, $rd2rd_flt, $filter_options);

        $job =  $self->newjob(
            name => "${name}_method_0",
            ifiles => [$crr_reads, $raw_reads, $prictg],
            ofiles => [$rd2rd_flt, $inconsistent, $readinfos], 
            mfiles => [],
            jobs => [$job_inconsistent, $job_filter],
            msg => "phasing reads",
        );

    } elsif ($method == 1) {
        ### output
        my $inconsistent = "$wrkdir/clair3/inconsistent";
        my $readinfos = "$wrkdir/clair3/readinfos";

        my $job_inconsistent = $self->newjob_identify_inconsistent_overlaps_clair3($name, $wrkdir . "/clair3");

        my $filter_options = $self->get_config("phase_clair3_filter_options");
        if ($self->get_config("phase_clair3_use_reads") + 0 == 0) {
            $filter_options = $filter_options . " --range $crr_reads";
        }
        my $job_filter = $self->newjob_filter_inconsistent($name."_clair3", $wrkdir . "/clair3", $rd2rd, $rd2rd_flt, $filter_options);

        $job =  $self->newjob(
            name => "${name}_method_1",
            ifiles => [$crr_reads, $raw_reads, $prictg],
            ofiles => [$rd2rd_flt, $inconsistent, $readinfos], 
            mfiles => [],
            jobs => [$job_inconsistent, $job_filter],
            msg => "phasing reads"
        );

    } elsif ($method == 2) {
        ### output
        my $inconsistent_fsa = "$wrkdir/fsa/inconsistent";
        my $readinfos_fsa = "$wrkdir/fsa/readinfos";
        my $inconsistent_clair3 = "$wrkdir/clair3/inconsistent";
        my $readinfos_clair3 = "$wrkdir/clair3/readinfos";

        my $job_inconsistent_fsa = $self->newjob_identify_inconsistent_overlaps_fsa($name, $wrkdir . "/fsa");
        my $job_inconsistent_clair3 = $self->newjob_identify_inconsistent_overlaps_clair3($name, $wrkdir . "/clair3");

        my $rd2rd_flt0 = "$wrkdir/clair3/filtered_overlaps.paf";        
        
        my $filter0_options = $self->get_config("phase_filter_options");
        if ($self->get_config("phase_use_reads") + 0 == 0) {
            $filter0_options = $filter0_options . " --range $crr_reads";
        }
        my $filter1_options = $self->get_config("phase_clair3_filter_options");
        if ($self->get_config("phase_clair3_use_reads") + 0 == 0) {
            $filter1_options = $filter1_options . " --range $crr_reads";
        }

        my $job_filter0 = $self->newjob_filter_inconsistent($name."_fsa", $wrkdir. "/fsa", $rd2rd, $rd2rd_flt0, $filter0_options);
        my $job_filter1 = $self->newjob_filter_inconsistent($name."_clair3", , $wrkdir . "/clair3", $rd2rd_flt0, $rd2rd_flt, $filter1_options);

        my $job_inconsistent =  $self->newjob(
            name => "${name}_inconsistent_2",
            ifiles => [$crr_reads, $raw_reads, $prictg],
            ofiles => [$inconsistent_fsa, $readinfos_fsa, $inconsistent_clair3, $readinfos_clair3], 
            mfiles => [],
            pjobs => [$job_inconsistent_fsa, $job_inconsistent_clair3],
            msg => "identifying inconsistent overlaps with fsa and clair3",
        );
        
        $job =  $self->newjob(
            name => "${name}_method_2",
            ifiles => [$crr_reads, $raw_reads, $prictg],
            ofiles => [$rd2rd_flt, $inconsistent_fsa, $readinfos_fsa, $inconsistent_clair3, $readinfos_clair3], 
            mfiles => [],
            jobs => [$job_inconsistent, $job_filter0, $job_filter1],
            msg => "phasing reads",
        );
    } else {
        Plgd::Logger::error("No phase method: $method");
    }
    
    $self->run_jobs($job);
}


sub newjob_identify_inconsistent_overlaps_fsa($$$) {
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
    my $inconsistent = "$wrkdir/inconsistent";
    my $readinfos = "$wrkdir/readinfos";

    my $job_map = $self->job_map_reads_to_contigs($name."_fsa", $wrkdir, $reads, [$prictg, $altctg], $self->get_config("PHASE_RD2CTG_OPTIONS"));
    
    my $job_phase = $self->job_phase($name."_fsa", $wrkdir, $reads, $prictg, "$wrkdir/rd2ctg.paf", $self->get_config("phase_phase_options"));

    return $self->newjob(
        name => "${name}_fsa_inconsitent",
        ifiles => [$reads, $prictg, $altctg],
        ofiles => [$inconsistent, $readinfos],
        mfiles => [$rd2ctg, "$wrkdir/variants", "$wrkdir/consistent", "$wrkdir/load.paf"],
        jobs => [$job_map, $job_phase],
        msg => "identifing inconsistent overlaps",
    );
}


sub newjob_identify_inconsistent_overlaps_clair3($$$) {
    my ($self, $name, $wrkdir) = @_;
    mkdir $wrkdir;

    my $subname = $name . "_clair3";

    my $isGz = $self->get_config("COMPRESS");
    my $prjdir = $self->get_project_folder();

    my $corrReads = $isGz ? "$prjdir/1-correct/corrected_reads.fasta.gz" : "$prjdir/1-correct/corrected_reads.fasta";
    my $prepReads = $isGz ? "$prjdir/0-prepare/prepared_reads.fasta.gz" : "$prjdir/0-prepare/prepared_reads.fasta";

    my $reads = $self->get_config("phase_clair3_use_reads") + 0 == 1 ? $corrReads : $prepReads;
    my $contigs = "$prjdir/3-assemble/primary.fasta";
    my $rd_2_ctg = "$wrkdir/rd_2_ctg.sam";
    my $vcf = "$wrkdir/merge_output.vcf.gz";
    
    my $inconsistent = "$wrkdir/inconsistent";
    my $readinfos = "$wrkdir/readinfos";

    my $job_call = $self->job_calling_with_clair3($subname, $wrkdir, $reads, $contigs,  "");    

    my $options = $self->get_config("phase_clair3_phase_options") . " --vcf_fname $vcf";
    my $job_phase = $self->job_phase($subname, $wrkdir, $reads, $contigs, $rd_2_ctg, $options);

    return $self->newjob(
        name => "${subname}_inconsitent",
        ifiles => [$reads, $contigs],
        ofiles => [$inconsistent, $readinfos],
        mfiles => [$rd_2_ctg, "$wrkdir/variants", "$wrkdir/consistent", "$wrkdir/load.paf"],
        jobs => [$job_call, $job_phase],
        msg => "identifing inconsistent overlaps with clair3",
    );
}


sub job_calling_with_clair3($$$$$) {
    my ($self, $name, $wrkdir, $reads, $contigs, $vcf) = @_;

    mkdir $wrkdir;
    mkdir "$wrkdir/clair3";
    my $prj_name = $self->get_config("project"); 
    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");
    
    my $rd_2_ctg = "$wrkdir/rd_2_ctg";
    my $map_options = $self->get_config("phase_clair3_rd2ctg_options");

    my $prjdir = $self->get_project_folder();
    my $clair3_options = $self->get_config("phase_clair3_options") ;

    my $clair3_cmd = "run_clair3.sh ";
    if ($self->get_config("phase_clair3_command") ne "") {
        $clair3_cmd = $self->get_config("phase_clair3_command");
    #     if (not $clair3_cmd =~ /\".*\"/) {
    #         $clair3_cmd = "\"$clair3_cmd\"";
    #     }
    }

#  --platform="ont" \               ## options: {ont,hifi,ilmn}
#  --model_path=${MODEL_PREFIX} \   ## absolute model path prefix

    my $job_al = $self->newjob(
        name => "${name}_rd_2_ctg",
        ifiles => [$reads, $contigs],
        ofiles => ["$rd_2_ctg.sam", "$rd_2_ctg.bam"],
        gfiles => ["$rd_2_ctg.bam", "$rd_2_ctg.sam", "$rd_2_ctg.bam.tmp.*.bam"],
        mfiles => [],
        cmds => ["minimap2 $map_options -a  -t $threads $contigs $reads > $rd_2_ctg.sam",
                 "samtools sort -\@$threads --output-fmt BAM -o $rd_2_ctg.bam $rd_2_ctg.sam",
                 "samtools index -\@$threads $rd_2_ctg.bam"],
        msg => "mapping reads to contigs",
    );

    my $job_call = $self->newjob(
        name => "${name}_call",
        ifiles => ["$rd_2_ctg.bam", $contigs],
        ofiles => ["$wrkdir/merge_output.vcf.gz"],
        gfiles => ["$wrkdir/clair3"],
        mfiles => ["$wrkdir/clair3"],
        cmds => ["samtools faidx $contigs",
                 "$clair3_cmd --bam_fn=$rd_2_ctg.bam --ref_fn=$contigs --threads=$threads --output=$wrkdir/clair3 $clair3_options", 
                 "cp $wrkdir/clair3/merge_output.vcf.gz $wrkdir/merge_output.vcf.gz"],
        msg => "calling variant with clair3",
    );
    
    return $self->newjob(
        name => "${name}_clair3",
        ifiles => [$reads, $contigs],
        ofiles => ["$wrkdir/merge_output.vcf.gz"],
        mfiles => ["$rd_2_ctg.bam*"],
        jobs => [$job_al, $job_call],
        msg => "calling variant with clair3",
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


    $self->runAssemble($name, $wrkdir, $reads, $overlaps, $self->get_config("ASM2_ASSEMBLE_OPTIONS"));
    
}

sub newjob_racon_prialt($$) {
    my ($self, $wrkdir) = @_;
    my $name = "pol";
    my $prjdir = $self->get_project_folder();
    
    my $ctg_pri = "$prjdir/5-assemble/primary.fasta";
    my $ctg_alt = "$prjdir/5-assemble/alternate.fasta";
    my $ctg_all = "$wrkdir/ctgall.fasta";

    my $tile_pri = "$prjdir/5-assemble/primary_tiles";
    my $tile_alt = "$prjdir/5-assemble/alternate_tiles";

    my $readinfo = "$prjdir/4-phase/readinfos";
    my $phase_method = $self->get_config("phase_method") + 0;
    if ($phase_method == 0) {
        $readinfo = "$prjdir/4-phase/fsa/readinfos";
    } elsif ($phase_method == 1 || $phase_method == 2) {
        $readinfo = "$prjdir/4-phase/clair3/readinfos";
    }


    my $pol_pri = "$wrkdir/primary.fasta";
    my $pol_alt = "$wrkdir/alternate.fasta";
    
    my $read_crr = "$prjdir/1-correct/corrected_reads.fasta";
    my $read_prp = "$prjdir/0-prepare/prepared_reads.fasta";
    my $reads = $self->get_config("polish_use_reads") + 0 == 1 ? $read_crr : $read_prp;

    my $map_options = $self->get_config("polish_map_options");
    my $filter_options = $self->get_config("polish_filter_options");
    my $cns_options = $self->get_config("polish_cns_options");

    my $rd2ctg = "$wrkdir/rd2ctg.paf";
    my $rd2ctg_flt = "$wrkdir/rd2ctg_flt.paf";
    my $rd2ctg_pri = "$wrkdir/rd2ctg_pri.paf";
    my $rd2ctg_alt = "$wrkdir/rd2ctg_alt.paf";
    my $itype = "paf";
    my $otype = "paf";
    if ($map_options=~/([ \t]|^)-[a-zAz]*a/) {
        $rd2ctg = "$wrkdir/rd2ctg.sam";
        $rd2ctg_flt = "$wrkdir/rd2ctg_flt.sam";
        $rd2ctg_pri = "$wrkdir/rd2ctg_pri.sam";
        $rd2ctg_alt = "$wrkdir/rd2ctg_alt.sam";
        $itype = "sam";
        $otype = "sam";
    }
    my $rd_2_pri_names = "$wrkdir/rd_2_pri_names";
    my $rd_2_alt_names = "$wrkdir/rd_2_alt_names";

    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");
    
    my $filter_cmd = "";
    if ($filter_options ne "") {
        $filter_cmd = " | $bin_path/fsa_ol_refine - - --itype $itype --otype $otype --thread_size 4 $filter_options ";
    }

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
        ofiles => [$rd2ctg_pri, $rd2ctg_alt, $rd_2_pri_names, $rd_2_alt_names],
        gfiles => [$rd2ctg_flt, $rd2ctg_pri, $rd2ctg_alt],
        mfiles => [$rd2ctg_flt],
        cmds => ["$bin_path/fxtools.py fx_purge_overlaps $rd2ctg $readinfo --tile $tile_pri --tile $tile_alt > $rd2ctg_flt",
                 "$bin_path/fxtools.py fx_split_mappings $rd2ctg_flt $ctg_pri,$ctg_alt $rd2ctg_pri,$rd2ctg_alt",
                 "awk '{print \$1}' $rd2ctg_pri > $rd_2_pri_names",
                 "awk '{print \$1}' $rd2ctg_alt > $rd_2_alt_names"],
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
                    mfiles => [$rd2ctg_pri],
                    cmds => ["if [ -s $ctg_pri ]; then racon $cns_options -t $threads $reads $rd2ctg_pri $ctg_pri > $pol_pri; else touch $pol_pri; fi"],
                    msg => "polishing primary contigs",
                  ), 
                  $self->newjob(
                    name => "${name}_polish_alt",
                    ifiles => [$ctg_alt, $reads, $rd2ctg_alt],
                    ofiles => [$pol_alt],
                    gfiles => [$pol_alt],
                    mfiles => [$rd2ctg_alt],
                    cmds => ["if [ -s $ctg_alt ]; then racon $cns_options -t $threads $reads $rd2ctg_alt $ctg_alt > $pol_alt; else touch $pol_alt; fi"],
                    msg => "polishing alternate contigs",

                  )
        ],
        msg => "polishing pri/alt format contigs",
    );

    return $self->newjob(
        name => "${name}_job",
        ifiles => [$ctg_pri, $ctg_alt, $reads],
        ofiles => [$pol_pri, $pol_alt, $rd_2_pri_names, $rd_2_alt_names],
        mfiles => [$rd2ctg_flt, $rd2ctg, $rd2ctg_pri, $rd2ctg_alt],
        jobs => [$job_map, $job_flt, $job_pol],
        msg => "polishing job for pri/alt contigs",
    );    
}

sub newjob_racon_dual($$) {
    my ($self, $wrkdir) = @_;
    my $name = "pol";
    my $prjdir = $self->get_project_folder();
    my $wrkdir_asm = $self->get_work_folder("5-assemble");
    
    my $asm2_options = $self->get_config("asm2_assemble_options");

    my $dual = index($asm2_options, "dual") != -1;
    my $prialt = index($asm2_options, "prialt") != -1 || (!$dual);

    my $ctg_pri = "$wrkdir_asm/haplotype_1.fasta";
    my $ctg_alt = "$wrkdir_asm/haplotype_2.fasta";
    my $ctg_all = "$wrkdir/ctgall_hap.fasta";

    my $tile_pri = "$wrkdir_asm/haplotype_1_tiles";
    my $tile_alt = "$wrkdir_asm/haplotype_2_tiles";

    my $readinfo = "$prjdir/4-phase/readinfos";
    my $readinfo = "$prjdir/4-phase/readinfos";
    my $phase_method = $self->get_config("phase_method") + 0;
    if ($phase_method == 0) {
        $readinfo = "$prjdir/4-phase/fsa/readinfos";
    } elsif ($phase_method == 1 || $phase_method == 2) {
        $readinfo = "$prjdir/4-phase/clair3/readinfos";
    }

    my $pol_pri = "$wrkdir/haplotype_1.fasta";
    my $pol_alt = "$wrkdir/haplotype_2.fasta";
    
    my $read_crr = "$prjdir/1-correct/corrected_reads.fasta";
    my $read_prp = "$prjdir/0-prepare/prepared_reads.fasta";
    #my $reads = $prepReads; # TODO corrReads
    my $reads = $self->get_config("polish_use_reads") + 0 == 1 ? $read_crr : $read_prp;

    my $map_options = $self->get_config("polish_map_options");
    my $filter_options = $self->get_config("polish_filter_options");
    my $cns_options = $self->get_config("polish_cns_options");

    my $rd2ctg = "$wrkdir/rd2ctg_hap.paf";
    my $rd2ctg_flt = "$wrkdir/rd2ctg_flt_hap.paf";
    my $rd2ctg_pri = "$wrkdir/rd2ctg_hap1.paf";
    my $rd2ctg_alt = "$wrkdir/rd2ctg_hap2.paf";
    my $itype = "paf";
    my $otype = "paf";

    if ($map_options=~/([ \t]|^)-[a-zAz]*a/) {
        $rd2ctg = "$wrkdir/rd2ctg_hap.sam";
        $rd2ctg_flt = "$wrkdir/rd2ctg_flt_hap.sam";
        $rd2ctg_pri = "$wrkdir/rd2ctg_hap1.sam";
        $rd2ctg_alt = "$wrkdir/rd2ctg_hap2.sam";
        $itype = "sam";
        $otype = "sam";
    }
    my $rd_2_pri_names = "$wrkdir/rd_2_hap1_names";
    my $rd_2_alt_names = "$wrkdir/rd_2_hap2_names";


    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");
    
    my $filter_cmd = "";
    if ($filter_options ne "") {
        $filter_cmd = " | $bin_path/fsa_ol_refine - - --itype $itype --otype $otype --thread_size 4 $filter_options ";
    }

    my $job_map = $self->newjob(
        name => "${name}_map_hap",
        ifiles => [$ctg_pri, $ctg_alt, $reads],
        ofiles => [$rd2ctg],
        gfiles => [$ctg_all, $rd2ctg],
        mfiles => [$ctg_all],
        cmds => ["cat $ctg_pri $ctg_alt > $ctg_all",
                 "minimap2  -t $threads $map_options $ctg_all $reads $filter_cmd > $rd2ctg"],
        msg => "mapping reads to contigs",
    );

    my $job_flt = $self->newjob(
        name => "${name}_filter_hap",
        ifiles => [$rd2ctg, $tile_pri, $tile_alt],
        ofiles => [$rd2ctg_pri, $rd2ctg_alt, $rd_2_pri_names, $rd_2_alt_names],
        gfiles => [$rd2ctg_flt, $rd2ctg_pri, $rd2ctg_alt],
        mfiles => [$rd2ctg_flt],
        cmds => ["$bin_path/fxtools.py fx_purge_overlaps $rd2ctg $readinfo --tile $tile_pri --tile $tile_alt > $rd2ctg_flt",
                 "$bin_path/fxtools.py fx_split_mappings $rd2ctg_flt $ctg_pri,$ctg_alt $rd2ctg_pri,$rd2ctg_alt",
                 "awk '{print \$1}' $rd2ctg_pri > $rd_2_pri_names",
                 "awk '{print \$1}' $rd2ctg_alt > $rd_2_alt_names"],
        msg => "filter overlaps in which the reads are different haplotype",
    );
    
    my $job_pol = $self->newjob(
        name => "${name}_polish_hap",
        ifiles => [$ctg_pri, $ctg_alt, $reads, $rd2ctg_pri, $rd2ctg_alt],
        ofiles => [$pol_pri, $pol_alt],
        mfiles => [],
        pjobs => [$self->newjob(
                    name => "${name}_polish_hap_1",
                    ifiles => [$ctg_pri, $reads, $rd2ctg_pri],
                    ofiles => [$pol_pri],
                    gfiles => [$pol_pri],
                    mfiles => [$rd2ctg_pri],
                    cmds => ["if [ -s $ctg_pri ]; then racon $cns_options -t $threads $reads $rd2ctg_pri $ctg_pri > $pol_pri; else touch $pol_pri; fi"],
                    msg => "polishing dual_1 contigs",
                  ), 
                  $self->newjob(
                    name => "${name}_polish_hap_2",
                    ifiles => [$ctg_alt, $reads, $rd2ctg_alt],
                    ofiles => [$pol_alt],
                    gfiles => [$pol_alt],
                    mfiles => [$rd2ctg_alt],
                    cmds => ["if [ -s $ctg_alt ]; then racon $cns_options -t $threads $reads $rd2ctg_alt $ctg_alt > $pol_alt; else touch $pol_alt; fi"],
                    msg => "polishing dual_2 contigs",

                  )
        ],
        msg => "polishing dual format contigs",
    );

    return $self->newjob(
        name => "${name}_job_hap",
        ifiles => [$ctg_pri, $ctg_alt, $reads],
        ofiles => [$pol_pri, $pol_alt, $rd_2_pri_names, $rd_2_alt_names],
        mfiles => [$rd2ctg_flt, $rd2ctg, $rd2ctg_pri, $rd2ctg_alt],
        jobs => [$job_map, $job_flt, $job_pol],
        msg => "polishing dual format contigs",
    );    
}

sub newjob_racon($) {
    my ($self) = @_;

    my $wrkdir = $self->get_work_folder("6-polish") . "/racon";
    mkdir $wrkdir;

    my $asm2_options = $self->get_config("asm2_assemble_options");

    my $dual = index($asm2_options, "dual") != -1;
    my $prialt = index($asm2_options, "prialt") != -1 || (!$dual);

    my @pjobs = ();
    my @ofiles = ();
    my @ifiles = ();

    if ($prialt) {
        my $job_prialt = $self->newjob_racon_prialt($wrkdir);
        push @pjobs, $job_prialt;
        push @ifiles, @{$job_prialt->{ifiles}};
        push @ofiles, @{$job_prialt->{ofiles}};

    }
    if ($dual) {
        my $job_dual = $self->newjob_racon_dual($wrkdir);
        push @pjobs, $job_dual;
        push @ifiles, @{$job_dual->{ifiles}};
        push @ofiles, @{$job_dual->{ofiles}};
    }
    
    my $job = $self->newjob(
        name => "pol_job_racon",
        ifiles => [@ifiles],
        ofiles => [@ofiles],
        mfiles => [],
        pjobs => [@pjobs],
        msg => "polishing contigs using racon",
    );
    return $job;
}
sub newjob_medaka_oneset($) {
    my ($self, $name, $wrkdir, $ctg, $ctg_tile, $alt_tile, $rd_2_ctg_names, $pol) = @_;

    mkdir $wrkdir;

    my $prjdir = $self->get_project_folder();

    my $readinfo = "$prjdir/4-phase/readinfos";
    my $phase_method = $self->get_config("phase_method") + 0;
    if ($phase_method == 0) {
        $readinfo = "$prjdir/4-phase/fsa/readinfos";
    } elsif ($phase_method == 1 || $phase_method == 2) {
        $readinfo = "$prjdir/4-phase/clair3/readinfos";
    }

    my $reads = "$prjdir/0-prepare/prepared_reads.fasta";

    my $map_options = $self->get_config("polish_medaka_map_options");
    my $filter_options = $self->get_config("polish_medaka_filter_options");
    my $cns_options = $self->get_config("polish_medaka_cns_options");
    if ($cns_options ne "") {
        if (not $cns_options =~ /\".*\"/) {
            $cns_options = " --options \"$cns_options\"";
        } else {
            $cns_options = " --options $cns_options";
        }
    }

    my $bin_path = $self->get_env("BinPath");
    my $threads = $self->get_config("threads");
        
    my $filter_cmd = "";
    if ($filter_options ne "") {
        $filter_cmd = " | $bin_path/fsa_ol_refine - - --itype sam --otype sam --thread_size 4 $filter_options ";
    }

    my $sub_reads = "$wrkdir/sub_reads.fasta";
    my $rd_2_ctg_sam = "$wrkdir/rd_2_ctg.sam";

    my $job_map = $self->newjob(
        name => "${name}_map",
        ifiles => [$ctg, $ctg_tile, $reads, $rd_2_ctg_names],
        ofiles => [$rd_2_ctg_sam],
        gfiles => [$rd_2_ctg_sam, $sub_reads],
        mfiles => [$sub_reads],
        cmds => ["$bin_path/fsa_rd_tools sub $reads $sub_reads --names_fname $rd_2_ctg_names",
                 "minimap2 -a -t $threads $map_options $ctg $sub_reads $filter_cmd > $rd_2_ctg_sam"],
        msg => "mapping reads to contigs ($name)",
    );

    my $rd_2_ctg_flt_sam = "$wrkdir/rd_2_ctg_flt.sam";
    my $rd_2_ctg_flt_bam = "$wrkdir/rd_2_ctg_flt.bam";
    my $rd_2_ctg_flt_sorted_bam = "$wrkdir/rd_2_ctg_flt_sorted.bam";

    my $job_flt = $self->newjob(
        name => "${name}_filter",
        ifiles => [$rd_2_ctg_sam, $ctg_tile, $alt_tile],
        ofiles => [$rd_2_ctg_flt_sorted_bam],
        gfiles => [$rd_2_ctg_flt_sam, $rd_2_ctg_flt_bam, $rd_2_ctg_flt_sorted_bam],
        mfiles => [$rd_2_ctg_flt_sam, $rd_2_ctg_flt_bam],
        cmds => ["$bin_path/fxtools.py fx_purge_overlaps $rd_2_ctg_sam $readinfo --tile $ctg_tile --tile $alt_tile > $rd_2_ctg_flt_sam",
                 "samtools view -bS -@ $threads $rd_2_ctg_flt_sam -o $rd_2_ctg_flt_bam",
                 "samtools sort -@ $threads $rd_2_ctg_flt_bam  -o $rd_2_ctg_flt_sorted_bam",
                 "samtools index -@ $threads $rd_2_ctg_flt_sorted_bam"],
        msg => "filter overlaps in which the reads are different haplotype  ($name)",
    );

    my $medaka_cmd = "medaka";
    if ($self->get_config("polish_medaka_command") ne "") {
        $medaka_cmd = $self->get_config("polish_medaka_command");
        if (not $medaka_cmd =~ /\".*\"/) {
            $medaka_cmd = "\"$medaka_cmd\"";
        }
    }

    mkdir "$wrkdir/hdf";

    my $job_pol = $self->newjob(
        name => "${name}_polish",
        ifiles => [$ctg, $rd_2_ctg_flt_sorted_bam],
        ofiles => [$pol],
        gfiles => [$pol, "$wrkdir/hdf/*", "$ctg.fai"],
        mfiles => ["$wrkdir/hdf", "$pol.*", "$ctg.*"],
        #cmds => ["medaka consensus --threads $threads --bam_workers 24 $cns_options $rd_2_ctg_flt_sorted_bam $wrkdir/ctg_hdf",
        #         "medaka stitch --threads $threads $wrkdir/ctg_hdf $wrkdir/draft/${name}.fasta $wrkdir/${name}.fasta"],
        cmds => ["if [ -s $ctg ]; then $bin_path/parallel_medaka.py $ctg $rd_2_ctg_flt_sorted_bam $pol --threads $threads $cns_options --wrkdir $wrkdir/hdf --medaka $medaka_cmd ; else touch $pol; fi"],
        msg => "running medaka ($name)",
    );

    return $self->newjob(
        name => "${name}_job",
        ifiles => [$ctg, $ctg_tile, $alt_tile, $rd_2_ctg_names],
        ofiles => [$pol],
        mfiles => ["$rd_2_ctg_flt_sorted_bam.*", "$sub_reads", "$rd_2_ctg_sam", "$rd_2_ctg_flt_sorted_bam"],
        jobs => [$job_map, $job_flt, $job_pol],
        msg => "polishing contigs ($name)",
    );
}
sub newjob_medaka($) {
    my ($self) = @_;

    my $prjdir = $self->get_project_folder();
    my $wrkdir = $self->get_work_folder("6-polish") . "/medaka";
    mkdir $wrkdir;

    my $asm2_options = $self->get_config("asm2_assemble_options");
    my $dual = index($asm2_options, "dual") != -1;
    my $prialt = index($asm2_options, "prialt") != -1 || (!$dual);

    my @pjobs = ();
    my @ofiles = ();
    my @ifiles = ();

    if ($prialt) {
        my $ctg_pri = "$prjdir/6-polish/racon/primary.fasta";
        my $ctg_alt = "$prjdir/6-polish/racon/alternate.fasta";
        my $ctg_pri_tile = "$prjdir/5-assemble/primary_tiles";
        my $ctg_alt_tile = "$prjdir/5-assemble/alternate_tiles";
        my $rd_2_pri_name = "$prjdir/6-polish/racon/rd_2_pri_names";
        my $rd_2_alt_name = "$prjdir/6-polish/racon/rd_2_alt_names";
        my $pol_pri = "$prjdir/6-polish/medaka/primary.fasta";
        my $pol_alt = "$prjdir/6-polish/medaka/alternate.fasta";
        my $job_pri = $self->newjob_medaka_oneset("pol_medaka_pri", "$wrkdir/pri", $ctg_pri, $ctg_pri_tile, $ctg_alt_tile, $rd_2_pri_name, $pol_pri);
        my $job_alt = $self->newjob_medaka_oneset("pol_medaka_alt", "$wrkdir/alt", $ctg_alt, $ctg_alt_tile, $ctg_pri_tile, $rd_2_alt_name, $pol_alt);

        push @pjobs, $job_pri, $job_alt;
        push @ifiles, @{$job_pri->{ifiles}}, @{$job_alt->{ifiles}};
        push @ofiles, @{$job_pri->{ofiles}}, @{$job_alt->{ofiles}};

    }
    if ($dual) {
        my $ctg_hap1 = "$prjdir/6-polish/racon/haplotype_1.fasta";
        my $ctg_hap2 = "$prjdir/6-polish/racon/haplotype_2.fasta";
        my $ctg_hap1_tile = "$prjdir/5-assemble/haplotype_1_tiles";
        my $ctg_hap2_tile = "$prjdir/5-assemble/haplotype_2_tiles";
        my $rd_2_hap1_name = "$prjdir/6-polish/racon/rd_2_hap1_names";
        my $rd_2_hap2_name = "$prjdir/6-polish/racon/rd_2_hap2_names";
        my $pol_hap1 = "$prjdir/6-polish/medaka/haplotype_1.fasta";
        my $pol_hap2 = "$prjdir/6-polish/medaka/haplotype_2.fasta";

        my $job_hap1 = $self->newjob_medaka_oneset("pol_medaka_hap1", "$wrkdir/hap1", $ctg_hap1, $ctg_hap1_tile, $ctg_hap2_tile, $rd_2_hap1_name, $pol_hap1);
        my $job_hap2 = $self->newjob_medaka_oneset("pol_medaka_hap2", "$wrkdir/hap2", $ctg_hap2, $ctg_hap2_tile, $ctg_hap1_tile, $rd_2_hap2_name, $pol_hap2);

        push @pjobs, $job_hap1, $job_hap2;
        push @ifiles, @{$job_hap1->{ifiles}}, @{$job_hap2->{ifiles}};
        push @ofiles, @{$job_hap1->{ofiles}}, @{$job_hap2->{ofiles}};
    }

    my $job = $self->newjob(
        name => "pol_job_medaka",
        ifiles => [],           # @ifiles
        ofiles => [@ofiles],    
        mfiles => ["$wrkdir/hap1", "$wrkdir/hap2", "$wrkdir/pri", "$wrkdir/alt"],
        pjobs => [@pjobs],
        msg => "polishing contigs using medaka",
    );
    return $job;
}


sub run_polish($) {
    my ($self) = @_;

    my $wrkdir = $self->get_work_folder("6-polish");
    mkdir $wrkdir;

    my @ofiles = ();
    my @jobs = ();

    my $job_racon = $self->newjob_racon();
    push @ofiles, @{$job_racon->{ofiles}};
    push @jobs, $job_racon;

    my $use_medaka = $self->get_config("polish_medaka") + 0;
    if ($use_medaka == 1) {
        my $job_medaka = $self->newjob_medaka();
        push @ofiles, @{$job_medaka->{ofiles}};
        push @jobs, $job_medaka;
    }



    my $job = $self->newjob(
        name => "pol_job",
        ifiles => [],
        ofiles => [@ofiles],
        mfiles => [],
        jobs => [@jobs],
        msg => "polishing contigs",
    );
    $self->run_jobs($job);
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
#    ["hic_reads", "", 0, "hic reads path"],
    ["genome_size", "", 1, "genome size"],
    ["threads", "4", 0],
#    ["memory", "0", 0],
    ["cleanup", "1", 0],
#    ["compress", "0", 0],
    ["grid", "local", 0],

    ["prep_min_length", "3000", 0],
    ["prep_output_coverage", "80", 0],

    ["corr_iterate_number", "1", 0],
    ["corr_block_size", "4000000000", 0],
    ["corr_correct_options", "", 0],
    ["corr_filter_options", "--filter0=l=5000:al=2500:alr=0.5:aal=5000:oh=3000:ohr=0.3", 0],
    ["corr_rd2rd_options", "-x ava-ont", 0],
    ["corr_output_coverage", "80", 0],

    ["align_block_size", "4000000000", 0],
    ["align_rd2rd_options", "-X -g3000 -w30 -k19 -m100 -r500", 0],
    ["align_filter_options", "--filter0=l=5000:aal=6000:aalr=0.5:oh=3000:ohr=0.3 --task=extend --filter1=oh=300:ohr=0.03", 0],
    
    ["asm1_assemble_options", "", 0],
    
    ["phase_method", "", 0],
    ["phase_rd2ctg_options", "-x map-ont -c -p 0.5 -r 1000", 0],
    ["phase_use_reads", "1", 0],
    ["phase_phase_options", "", 0],
    ["phase_clair3_command", "singularity exec --containall -B `pwd -P`:`pwd -P` clair3_v0.1-r12.sif /opt/bin/run_clair3.sh", 0],
    ["phase_clair3_use_reads", "0", 0],
    ["phase_clair3_options", "--platform=ont --model_path=/opt/models/ont_guppy5/  --include_all_ctgs", 0],
    ["phase_clair3_rd2ctg_options", "-x map-ont -w10 -k19 -c -p 0.5 -r 1000", 0],
    ["phase_clair3_phase_options", "--coverage lc=30 --phase_options icr=0.1:icc=6:sc=10 --filter i=70"],
    ["phase_clair3_filter_options", "--threshold=2500 --rate 0.05", 0],

    ["asm2_assemble_options", "", 0],

    ["polish_map_options", "-x map-ont", 0], #  --secondary=no
    ["polish_filter_options", "--filter0 oh=1000:ohr=0.1", 0],
    ["polish_cns_options", "", 0],
    ["polish_medaka", "", 0], #  
    ["polish_medaka_command", "singularity exec --containall -B `pwd -P`:`pwd -P` medaka_1.7.2--aa54076.sif medaka", 0],
    ["polish_medaka_map_options", "-x map-ont", 0],
    ["polish_medaka_cns_options", "--model r941_prom_sup_g507", 0],
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
          "    unzip:       generate diploid assembly\n" .
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
