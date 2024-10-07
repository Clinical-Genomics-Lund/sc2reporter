#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use vcf2;
use File::Basename qw(basename);
use Cwd qw(abs_path);
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;
use DateTime;

# my @VARIANTS_OF_BIOLOGICAL_SIGNIFICANCE = ("S:N501Y","S:N501T", "S:N501S", "S:E484K", "S:K417T", "S:F157L", "S:V367F", "S:Q613H", "S:P681R", "S:Q677H", "S:F888L", "S:H69_V70del", "S:N439K", "S:Y453F", "S:S98F", "S:L452R", "S:D80Y", "S:A626S", "S:V1122L", "S:A222V", "S:S477N");
my @VARIANTS_OF_BIOLOGICAL_SIGNIFICANCE = ("S:L18F","S:H69_V70del","S:D80Y","S:S98F","S:F157L","S:A222V","S:V367F","S:K417N","S:K417T","S:N439K","S:L452R","S:Y453F","S:S477N","S:E484K","S:S494P","S:N501Y","S:N501T","S:N501S","S:Q613H","S:H655Y","S:P681H","S:A626S","S:Q677H","S:P681R","S:F888L","S:V1122L","ORF1ab:S3675_F3677del");

my @POSITIONS_OF_BIOLOGICAL_SIGNIFICANCE = ("S:501", "S:484");


# Connect to database, and create handles for collections
my $client = MongoDB->connect();
my $VARIANT = $client->ns("sarscov2.variant");
$VARIANT = $VARIANT->with_codec( prefer_numeric => 1 );

# Get all variants that are already defined in the database
my %vars_in_db = fetch_variants();

# Find input variant files
my $dir = $ARGV[0];
my @vcfs = glob "$dir/*.freebayes.vep.vcf";

# Read pangolin data if it exists
my %pangolin;
if( -e "$dir/pangolin_all.csv" ) {
    %pangolin = get_pangolin_data("$dir/pangolin_all.csv");
}

print STDERR "Parsing data...\n";
print "id\tanalysis_code\tpango_lineage\tvariants_of_significance\n";

foreach my $vcf_fn ( @vcfs ) {
    my ($sample_id) = (split /\./, basename($vcf_fn))[0];

    # Parse QC data
    my $qc_data;
    if( -e "$dir/$sample_id.qc.csv" and -e "$dir/$sample_id.flagstat" ) {
	$qc_data = read_qc_data("$dir/$sample_id.qc.csv");
	$qc_data->{on_target} = get_ontarget_from_flagstat("$dir/$sample_id.flagstat");
    } else {
	print STDERR "QC data not found for $sample_id!\n";
    }

    # Parse and add variant information
    my @variants;
    my @sample_variants = get_variants_from_vcf_vep($vcf_fn);
    foreach my $var (@sample_variants) {
	my $aa = $var->{var}->{csq}->{MUTATION};
	
	if( grep(/^$aa$/, @VARIANTS_OF_BIOLOGICAL_SIGNIFICANCE) ) {
	    push @variants, $aa;
	}
	else {
	    my ($gene, $change) = split /:/, $aa;

	    while( $change =~ /(\d+)/g ) {
		my $position = $gene.":".$1;
		if( grep(/^$position$/, @POSITIONS_OF_BIOLOGICAL_SIGNIFICANCE )) {
		    push @variants, $aa;
		    last;
		}
	    }
	    
	}
    }
    
    
    # Check if QC passed
    my $qc = ($qc_data->{pct_N_bases} < 10 ? "passed_qc" : "reextraction");

    # Print data
    if ($qc eq "passed_qc") {
	my $var_str = join("|", @variants);
	$var_str = "" unless $var_str;
	print $sample_id."\t"."NGSRESULTS"."\t".$pangolin{$sample_id}->{type}."\t".$var_str."\n";
    }
    print $sample_id."\t"."NGSQC"."\t".$qc."\t\n";
}


#############################################################################################
#############################################################################################
#############################################################################################




sub get_variants_from_vcf_vep {
    my $fn = shift;

    my $vcf = CMD::vcf2->new('file'=>$fn );

    my @vars;
    while ( my $v = $vcf->next_var() ) {
	my $id = $v->{POS}."_".$v->{REF}.">".$v->{ALT};
	my $csq = $v->{INFO}->{CSQ}->[0];
	delete $csq->{'MN908947.3.gff.gz'};
	my $dp = $v->{INFO}->{DP};
	my $alt_freq = $v->{GT}->[0]->{AO} / $dp;
	$csq->{MUTATION} = $csq->{SYMBOL}.":".simplify_hgvs($csq->{HGVSp});
	push( @vars, {'var'=>{'csq' => $csq, '_id'=>$id}, 'sample'=>{'dp'=>$dp, 'alt_freq'=>$alt_freq}} );
    }
    return @vars;
}


sub simplify_hgvs {
    my %short = ('Cys'=> 'C', 'Asp'=> 'D', 'Ser'=> 'S', 'Gln'=> 'Q', 'Lys'=> 'K',
		 'Ile'=> 'I', 'Pro'=> 'P', 'Thr'=> 'T', 'Phe'=> 'F', 'Asn'=> 'N',
		 'Gly'=> 'G', 'His'=> 'H', 'Leu'=> 'L', 'Arg'=> 'R', 'Trp'=> 'W',
		 'Ala'=> 'A', 'Val'=> 'V', 'Glu'=> 'E', 'Tyr'=> 'Y', 'Met'=> 'M', 'Ter'=> '*' );
    
    my $s = shift;
    return "" unless $s;
    my ($tid, $change) = split /:/, $s;
    $change =~ s/^[cp]\.//;
    for my $long (keys %short) {
	$change =~ s/$long/$short{$long}/g;
    }
    return $change;	
}

sub fetch_variants {
    my %vars;
    my $res = $VARIANT->find()->fields({'_id'=>1});
    while(my $var = $res->next) {
	$vars{$var->{_id}} = 1;
    }
    return %vars;
}
 

sub get_pangolin_data {
    my $fn = shift;
    open (my $fh, $fn);
    my $header = <$fh>;
    my %pangolin;
    while(<$fh>) {
	chomp;
	my @a = split /,/;
# 	print @a;
	my($id) = ($a[0] =~ /Consensus_(.*?)\.consensus_threshold/);
	die "could not extract sample ID from $a[0]" unless $id;
	$pangolin{$id}->{type} = $a[1]; 
	$pangolin{$id}->{probability} = $a[2];
	$pangolin{$id}->{pangolearn_version} = $a[3];
    }
    return %pangolin;
}
    
sub read_qc_data{
    my $fn = shift;
    my @data = read_csv($fn);
    return $data[0];

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


sub get_ontarget_from_flagstat {
    my $fn = shift;
    open(my $fh, $fn);
    while(<$fh>) {
	chomp;
	s/\r//;
	
	if (/mapped \((.*?)\% :/) {
	    return $1;
	}
    }	
    return "N/A"
}



