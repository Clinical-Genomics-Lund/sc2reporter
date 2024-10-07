#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename qw(basename dirname);
use POSIX qw( strftime);

#sftp -i .ssh/id_rsa_gensam se120@gensam-sftp.folkhalsomyndigheten.se

# ARGUMENTS: 1: pipeline output folder, 2: gisaid upload log file
my $in_dir = $ARGV[0];
my $gisaid_log = $ARGV[1];
my $metadata_csv = $ARGV[2]; 

my $SCRIPT_ROOT = dirname($0);

my %config = read_config($SCRIPT_ROOT.'/fohm.config');
my %id_conversion_table = read_conversion_table($SCRIPT_ROOT.'/conversion.table');
my %gisaid_ids = read_gisaid_ids($gisaid_log);
my %metadata = read_pat_metadata($metadata_csv);
my $full_out_dir = $config{out_dir}.'/'.basename($in_dir);
my $prefix_base = $config{region_code}."_".$config{lab_code};
my $prefix = $full_out_dir.'/'.$prefix_base;
system("mkdir -p $full_out_dir");

my @vcfs = glob "$in_dir/*.freebayes.vep.vcf";

# Get pangolin file
my $pangolin_fn = "$in_dir/pangolin_all.csv";
die "No pangolin output!" unless -e $pangolin_fn;

my $date = strftime '%Y-%m-%d', localtime;
   
open(CSV, ">".$config{out_dir}.'/'.$config{region_code}."_".$config{lab_code}."_".$date."_komplettering.csv");
#print CSV "provnummer,urvalskriterium,GISAID_accession\n";
my %excluded_samples;

foreach my $vcf_fn ( @vcfs ) {
    my ($sample_id) = (split /\./, basename($vcf_fn))[0];
    next if $sample_id eq "NTC" or $sample_id =~ "No_Sample" or $sample_id eq "No" or $sample_id =~ /NegativeControl/;
    my $mlu_id = $metadata{$sample_id}->{SID};
    unless ($mlu_id) {
	$excluded_samples{$sample_id} = 1;
	print STDERR "No SID (MLU-ID) found for $sample_id. Won't upload!\n";
	next;
    }
    # Parse QC data
    my $qc_data;
    if( -e "$in_dir/$sample_id.qc.csv") {
	$qc_data = read_qc_data("$in_dir/$sample_id.qc.csv");
    } else {
	die "QC data not found for $sample_id!\n";
    }

    # Check if QC passed
    if( $qc_data->{pct_N_bases} > 40 ) {
	$excluded_samples{$sample_id} = 1;
	next;
    }
	

    print CSV "$mlu_id,".($metadata{$sample_id}->{Resultat} or "Information saknas").",".$gisaid_ids{$id_conversion_table{$sample_id}}."\n";
    print "$mlu_id,".($metadata{$sample_id}->{Resultat} or "Information saknas").",".$gisaid_ids{$id_conversion_table{$sample_id}}."\n";
    
    my $fa_fn = "$in_dir/$sample_id.consensus.fa";
    die "No fasta for $sample_id" unless -e $fa_fn;

    my $fq1 = "$in_dir/${sample_id}_subsample_R1_001.fastq.gz";
    my $fq2 = "$in_dir/${sample_id}_subsample_R2_001.fastq.gz";
    die "Fastq files missing for $sample_id!" if( ! -e "$in_dir/${sample_id}_R1_001.fastq.gz" or ! -e "$in_dir/${sample_id}_R2_001.fastq.gz" );

    copy_and_fix_fasta($fa_fn, "${prefix}_${mlu_id}.consensus.fasta", "${prefix_base}_${mlu_id}");
    system "ln -sf $vcf_fn ${prefix}_${mlu_id}.vcf";
    system "ln -sf $fq1 ${prefix}_${mlu_id}_1.fastq.gz";
    system "ln -sf $fq2 ${prefix}_${mlu_id}_2.fastq.gz";
    
}

my $pangolin_format = detect_pangolin_format($pangolin_fn);
reformat_pangolin($pangolin_fn, "${prefix}_${date}_pangolin_classification$pangolin_format.txt", \%metadata, \%excluded_samples);


#############################################################################################
#############################################################################################
#############################################################################################

sub detect_pangolin_format {
    my $fn = shift;
    open(my $fh, $fn);
    my $header = <$fh>;

    return "_format4" if $header =~ /scorpio_call/;
    return "_format2" if $header =~ /conflict/;
    return "";
}
	

sub copy_and_fix_fasta {
    my ($orig_file, $new_file, $new_id) = @_;
    open(my $orig_fh, $orig_file) or die "cannot read: $orig_file\n";
    open(my $new_fh, '>'.$new_file) or die "cannot create file: $new_file\n";
    while(<$orig_fh>) {
	if( /^>/ ) {
	    print $new_fh ">$new_id\n";
	}
	else {
	    print $new_fh $_;
	}
    }
    close $new_fh;
    close $orig_fh;
}

sub reformat_pangolin {
    my ($in_fn, $out_fn, $metadata, $excluded) = @_;
    open(my $in, $in_fn);
    open(my $out, ">".$out_fn);
    my $header = <$in>;
    print $out $header;
    while(my $line = <$in>) {
	my ($old_id) = ($line =~ /^.*?_(.*?)\./);
	next if $excluded->{$old_id};
	my $new_id;
	if($metadata->{$old_id}->{SID}) {
	    $new_id =  $metadata->{$old_id}->{SID};
	} else {
	    print STDERR "No SID (MLU ID) found for $old_id! Removing from pangolin file\n";
	    next;
	}
    
	$line =~ s/^.*?_(.*?)\..*?,/$new_id,/; # Remove the stuff surrounding the ID.
	print $out $line;
    }
    close $in;
    close $out;
}
    
sub read_qc_data{
    my $fn = shift;
    my @data = read_csv($fn, ',');
    return $data[0];

}

sub read_csv {
    my $fn = shift;
    my $sep = shift;
    open (my $fh, $fn);
    chomp(my $header = <$fh>);
    $header =~ s/\r//;
    my @header = split /$sep/, $header;
    
    my @data;
    while(<$fh>) {
	chomp;
	s/\r//;
	my @a = split /$sep/;
	my %entry;
	for my $i (0..$#header) {
	    $entry{$header[$i]} = $a[$i];
	}
	push @data, \%entry;
    }
    return @data;
    
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

sub read_conversion_table {
    my $fn = shift;
    open(my $fh, $fn);
    my %table;
    while(<$fh>) {
	chomp;
	my ($pseudo, $real) = split /\t/;
	$table{$real} = $pseudo;
    }
    return %table;
}

sub read_gisaid_ids {
    my $fn = shift;
    open(my $fh, $fn);
    my %table;
    while(<$fh>) {
    chomp;
    my @line = split '"';
    next if not @line[3] eq 'epi_isl_id';
    my ($name, $gisaid) = split '; ', @line[7];
	my @a = split /\//, $name;
	my $pseudo_id = $a[2];
	$pseudo_id =~ /(0)*(\d+)$/;
	my $pseudo_num = $2;
	$table{$pseudo_num} = $gisaid;
    }
    return %table;
}

sub read_pat_metadata {
    my $fn = shift;
    my @csv = read_csv("iconv -f iso-8859-1 -t UTF-8 '$fn'|", ";");
    my %csv;
    foreach my $entry (@csv) {
	$csv{$entry->{Labbnummer}} = $entry;
    }
    return %csv;
}
