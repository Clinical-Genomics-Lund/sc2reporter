#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
use Data::Dumper;

my $PIPELINE = "/data/bnf/sarscov2/pipeline/example_wrapper.pl";
my $OUTPUT_DIR = "/data/bnf/sarscov2/results";
my $SCRIPT_ROOT = dirname($0);
my $FASTQ_SUBFOLDER = "Data/Intensities/BaseCalls/";
my $POST_ANALYSIS_SCRIPT = $SCRIPT_ROOT."/post_analysis_script.pl";

my $run_folder = $ARGV[0];

if( is_sarscov2_run($run_folder) ) {
    LOG("Starting pipeline for $run_folder");
    my $output = $OUTPUT_DIR.'/'.basename($run_folder);
	my $on_hopper =  ( $run_folder =~ /^\/fs2/ ? 1 : 0 );
	my $pipeline_log;
	if  ( $run_folder =~ /^\/fs2/ ) {
	  $pipeline_log = $OUTPUT_DIR."/hopperlog/".basename($run_folder)."-pipeline.log";
	}
	else {
	  $pipeline_log = "$output/sarscov2.log";
	}

    system("$PIPELINE $run_folder/$FASTQ_SUBFOLDER $output 24 $POST_ANALYSIS_SCRIPT|bash > $pipeline_log");
}
else {
    print STDERR "NOT A SARSCOV2 RUN\n";
}

sub is_sarscov2_run {
    my $run_folder = shift;

    return 0 unless -e $run_folder.'/SampleSheet.csv';
    open(my $ss, $run_folder.'/SampleSheet.csv');
    my $part = "pre";
    while(<$ss>) {
	chomp;
	s/\r//;
	if( /^\[(.*?)\]/ ) {
	    $part = $1;
	    next;
	}

	if( $part eq "Data" ) {
	    my @data = split /,/;
	    return 1 if $data[-1] eq "sarscov2";
	}
    }
    return 0
}

sub LOG{
    my $msg = shift;
    open( my $log_fh, ">>". $SCRIPT_ROOT.'/start_pipeline_if_sarscov2.log' );
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
    my $timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
    print $log_fh $timestamp.": $msg\n";
}
