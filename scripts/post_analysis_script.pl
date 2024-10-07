#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
use Data::Dumper;

my $run_folder = $ARGV[0];
my $RESULTS_DIR = "/data/bnf/sarscov2/results";
my $SCRIPT_ROOT = dirname($0);

my $metadata_filemask = "/media/mikro-cmd/SARS-CoV-2/pat_metadata/*.csv";
#my $metadata_filemask = "/data/bnf/sarscov2/scripts/*.csv";

# Check for missing data in run
LOG("Checking if all isolates have been demultiplexed in $run_folder");
# print($SCRIPT_ROOT."/rawdatacheck.sh $run_folder\n");
system("$SCRIPT_ROOT/rawdatacheck.sh $run_folder");

# Add results to reporting database
LOG("Adding $run_folder to database");
print($SCRIPT_ROOT."/load_run_to_db.pl '$run_folder'\n");

# Create results for file KM LIMS import
my $lims_file = "$RESULTS_DIR/csv_to_kmlims/".basename($run_folder).".txt";
LOG("Creating KMLIMS results file for $run_folder in $lims_file");
print($SCRIPT_ROOT."/create_run_summary_tsv.pl '$run_folder' > $lims_file\n");

# Try to find a matching metadata file
my $metadata_file = find_matching_metadata_file($metadata_filemask, $run_folder);
if( $metadata_file ) {

    # Add collection date to reporting database
    LOG("Adding metadata from $metadata_file corresponding to $run_folder");
    print($SCRIPT_ROOT."/add_metadata_to_db.pl '$metadata_file'\n");

    # Create data for GISAID
    my $run_id = basename($run_folder);
    my $gisaid_outdir = $RESULTS_DIR.'/gisaid/';
    LOG("Preparing data for GISAID upload");
    print($SCRIPT_ROOT."/gisaid.pl '$metadata_file' '$run_folder' '$gisaid_outdir' REAL > $gisaid_outdir/gisaid_".$run_id.".csv\n");

    # Upload data to GISAID
    LOG("Uploading data to GISAID");
    print("ml conda/miniconda3; source activate cli3_env; cli3 upload --fasta $RESULTS_DIR/gisaid/all_sequences_".$run_id.".fasta --metadata $RESULTS_DIR/gisaid/gisaid_".$run_id.".csv --failed $RESULTS_DIR/gisaid/gisaid_".$run_id.".failed.csv --log $RESULTS_DIR/gisaid/gisaid_".$run_id.".log\n");

    # Prepare data for FOHM
    LOG("Preparing data for FOHM upload");
    print($SCRIPT_ROOT."/fohm_upload.pl '$run_folder' '$gisaid_outdir/gisaid_".$run_id.".log' '$metadata_file'\n");

    # Running UShER
    LOG("Running UShER on all data");
    print($SCRIPT_ROOT."/dousher.sh\n");


    # TODO: Actually upload data to FOHM
    
} else {
    LOG("Not matching metadata file found for $run_folder!");
}




sub LOG{
    my $msg = shift;
    open( my $log_fh, ">>". $SCRIPT_ROOT.'/post_analysis_script.log' );
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime();
    my $timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
    print $log_fh $timestamp.": $msg\n";
}

sub find_matching_metadata_file {
    my $mask = shift;
    my $run_folder = shift;

    my %samples_in_run_folder = get_sample_ids_from_runfolder($run_folder);

    my @metadata_files = glob $mask;
    my @metadata_pwd = glob("*.qc.csv");
    push(@metadata_files, @metadata_pwd);
	my $maxfn = '';
	my $maxfound = 0;
	my %maxfound_samples;
    foreach my $fn ( @metadata_files ) {
	open(my $fh, $fn);
	my $found = 0;
	my %found_samples;
	# print $fn."\n";
	while(<$fh>) {
	    my @a = split /;/;
	    $found++ if $samples_in_run_folder{$a[1]};
	    $found_samples{$a[1]}=1;
	}
	if( $found >= keys %samples_in_run_folder ) {
	    return $fn;
	}
	elsif($found >= $maxfound) {
	  $maxfound = $found;
	  $maxfn = $fn;
	  %maxfound_samples = %found_samples;
	}
	}
	if( $maxfound > 0 ) {
	    print "in $maxfn, found: $maxfound out of ".scalar(keys %samples_in_run_folder)."\n";
	    foreach( keys %samples_in_run_folder ) {
		print "missing: $_\n" unless $maxfound_samples{$_};
	    }
	    return $maxfn;
    }
    return 0;
}


sub get_sample_ids_from_runfolder {
    my $folder = shift;
    my @qc_files = glob($folder.'/*.qc.csv');
    my %ids;
    foreach my $fn ( @qc_files ) {
	my( $id ) = ( basename($fn) =~ /^(.*?)\.qc.csv/ );
	next if $id eq "NTC" or $id =~ /No_sample/i or $id eq 'No' or $id =~ /NegativeControl/;
	$ids{$id}=1;
    }
#	print Dumper(\%ids);
    return %ids;
}
