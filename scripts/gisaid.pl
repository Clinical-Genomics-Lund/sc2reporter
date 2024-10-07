#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename qw(basename dirname);

my $SCRIPT_ROOT = dirname($0);
my $first_id = get_first_id("$SCRIPT_ROOT/conversion.table");
my %already_uploaded = get_already_uploaded_ids("$SCRIPT_ROOT/conversion.table");
my %config = read_config("$SCRIPT_ROOT/gisaid.config");

print "submitter,fn,covv_virus_name,covv_type,covv_passage,covv_collection_date,covv_location,covv_add_location,covv_host,covv_add_host_info,covv_gender,covv_patient_age,covv_patient_status,covv_specimen,covv_outbreak,covv_last_vaccinated,covv_treatment,covv_seq_technology,covv_assembly_method,covv_coverage,covv_orig_lab,covv_orig_lab_addr,covv_provider_sample_id,covv_subm_lab,covv_subm_lab_addr,covv_subm_sample_id,covv_authors\n";

my $METADATA_FILE = $ARGV[0];
my $IN_DIR = $ARGV[1];
my $OUT_DIR = $ARGV[2];
my $REAL = ($ARGV[3] or "");
my $fasta_fn = $OUT_DIR."/all_sequences_".basename($IN_DIR).".fasta";

open(CONV, ">>"."$SCRIPT_ROOT/conversion.table");
open(CLIN, "iconv -f iso-8859-1 -t UTF-8 '$METADATA_FILE'|");
open(ALL_FASTA, ">".$fasta_fn);
<CLIN>;

my $pseudo_id = $first_id;
while(<CLIN>) {
    chomp;
    s/\r//;
    my( $collection_date, $sample_id, $origin, $city, $gender, $age, $analysis, $lab, $criteria, $fohm ) = split /;/;

    if( $already_uploaded{$sample_id} ) {
	print STDERR "$sample_id already uploaded, skipping!\n";
	next;
    }
    
    if( qc_pass($IN_DIR, $sample_id) ) {
    
	print join(",", ($config{submitter},
			 basename($fasta_fn),
			 long_id($pseudo_id, $collection_date),
			 $config{type},
			 "Original", # FIXME
			 date($collection_date), 
			 location($origin), 
			 "", # Additonal location (e.g. Cruise Ship, Convention, Live animal market)
			 $config{host},
			 "", # Additional host information (e.g. Patient infected while traveling in...)
			 gender($gender),
			 age($age),
			 patient_status(),
			 specimen_source(),
			 outbreak(),
			 last_vaccinated(),
			 treatment(),
			 $config{sequencing_technology}, # FIXME: Take as argument
			 $config{assembly_method},
			 coverage($IN_DIR, $sample_id),
			 $config{originating_lab},
			 $config{originating_lab_address},
			 "", # Sample ID given by the sample provider
			 $config{submitting_lab},
			 $config{submitting_lab_address},
			 "", # Sample ID given by the submitting lab
			 $config{authors})
	    )."\n";
	print ALL_FASTA consensus_sequence($IN_DIR, $sample_id, $pseudo_id, $collection_date);
	print CONV "$pseudo_id\t$sample_id\n" if $REAL eq "REAL";
	$pseudo_id++;
    }

}

sub read_config {
    my $fn = shift;
    open(my $fh, $fn);
    my %config;
    while(<$fh>) {
	chomp;
	my ($key, $value) = split /=/;
	$config{$key} = $value;
    }
    return %config;
}

sub gender {
    my $str = shift;
    return "Male" if uc $str eq "MAN";    
    return "Female" if uc $str eq "KVINNA";
    return "unknown";
}

sub date{
    return (split ' ', $_[0])[0]
}

sub year{
    return (split '-', $_[0])[0]
}

sub age{
    return $_[0];
}

sub patient_status {
    return "unknown";
}

sub specimen_source {
    return "";
}

sub outbreak {
    return "";
}

sub last_vaccinated {
    return "";
}

sub treatment {
    return "";
}

sub coverage {
    my($dir, $id) = @_;
    open(DEPTH, $dir."/".$id.".depth") or return "unknown";
    my($tot_cov, $npos) = (0,0);
    <DEPTH>;
    while(<DEPTH>) {
	my @dp = split /\t/;
	$tot_cov += $dp[2];
	$npos++;
    }
    return sprintf("%dx", $tot_cov/$npos);
    
}

sub consensus_sequence {
    my($dir, $id, $pseudo_id, $date) = @_;
    open(CONS, $dir."/".$id.".consensus.fa") or die "NO SEQUENCE FOUND FOR $id";
#    open(CONS, $dir."/".$id.".consensus.fa") or return ">noseq\n";
    my $header = <CONS>;
    my $seq = "";
    while(<CONS>) {
	chomp;
	$seq .= $_;
    }
    return ">".long_id($pseudo_id, $date)."\n$seq\n";
}

sub long_id {
    my $full_pseudo_id = $config{'pseudo_id_prefix'}.sprintf("%07d", $_[0]);
	
    return "hCoV-19/Sweden/M-$full_pseudo_id/".year($_[1]); # 230627 adding "M-" as a iso code for Skane /J
}

sub location {
    my $region = shift;
    if( $region ) {
	$region =~ s/Å/a/g;
	$region =~ s/Ä/a/g;
	$region =~ s/Ö/o/g;
	$region =  ucfirst(lc $region);
    }
    else {
	$region = "Skane";
    }
    return "Europe / Sweden / $region";
}

sub qc_pass {
    my $dir = shift;
    my $id = shift;
    if( -e $dir."/".$id.".qc.csv" ) {
	my @qc = read_csv($dir."/".$id.".qc.csv");
	return 1 if $qc[0]->{pct_N_bases} <= 10;
	print STDERR "FAILED QC: $id\n";
	return 0;
    }else {
	print STDERR "WARNING: No qc data found for $id\n";
	return 0;
    }

}

sub read_csv {
    my $fn = shift;
    open (my $fh, $fn);
    chomp(my $header = <$fh>);
    $header =~ s/\r//;
    my @header = split /,/, $header;

    my @data;
    while(<$fh>) {
	chomp;
	s/\r//;
	my @a = split /,/;
	my %entry;
	for my $i (0..$#header) {
	    $entry{$header[$i]} = $a[$i];
	}
	push @data, \%entry;
    }
    return @data;

}

sub get_first_id {
    my $fn = shift;
    open(LIST, $fn) or return "1";

    my $max_id = 0;
    while(<LIST>) {
	chomp;
	my( $id, $sample_id) = split /\t/;
	if($id > $max_id) {
	    $max_id = $id;
	}
    }
    return $max_id + 1;
}

sub get_already_uploaded_ids {
    my $fn = shift;
    my %ids;
    open(LIST, $fn) or return %ids;

    while(<LIST>) {
	chomp;
	my( $id, $sample_id) = split /\t/;
	$ids{$sample_id} = 1;
    }
    return %ids;
}
