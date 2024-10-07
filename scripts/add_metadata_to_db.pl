#!/usr/bin/perl -w
use strict;
use MongoDB;
use MongoDB::BSON;
use MongoDB::OID;

my $sequencing_facility = "CMD";

my $db = ($ARGV[1] or "sarscov2");
print STDERR "Loading data to $db!\n";

# Connect to database, and create handles for collections
my $client = MongoDB->connect();
my $SAMPLE = $client->ns("$db.sample");

my %metadata = read_pat_metadata($ARGV[0]);

foreach my $sample_id (keys %metadata) {
	
    my $collection_date = $metadata{$sample_id}->{Provtagningsdatum};
    my $selection_raw = ($metadata{$sample_id}->{Urval} or $metadata{$sample_id}->{Resultat} or "Okänd");
    my $selection_criterion = htmlencode($selection_raw);
    my $age = ( $metadata{$sample_id}->{'Ålder'} or "unknown" );

    my $sex = "unknown";
    if( $metadata{$sample_id}->{'Kön'} eq "MAN" ) {
	$sex = "male";
    }
    elsif( $metadata{$sample_id}->{'Kön'} eq "KVINNA" ) {
	$sex = "female";
    }

    my $mlu = ($metadata{$sample_id}->{SID} or "");
#    my $lab = ($metadata{$sample_id}->{Info2} =~ /^DC/ ? "DynamicCode" : "Mikro");
    my $lab = ($metadata{$sample_id}->{Labbnummer} =~ /^DC/ ? "DynamicCode" : "Mikro");
    if ($metadata{$sample_id}->{Info2} =~ /^DC/) {
	  $lab = "DynamicCode" ;
	}

    my $date_fmt = date($collection_date);
    if( $date_fmt =~ /^\d\d\d\d-\d\d-\d\d$/ ) {

	my( $year,$month,$day ) = split /-/, $date_fmt;
	my $date_obj = DateTime->new(
	    year => $year, month => $month, day => $day,
	    hour => 12, minute => 0, second => 0, time_zone => 'Europe/Paris');

	my $add_metadata_res = $SAMPLE->update_many({'sample_id'=>$sample_id}, {'$set'=>{'collection_date'=>$date_obj, 'selection_criterion'=>$selection_criterion, 'sex'=>$sex, 'age'=>$age, 'lab'=>$lab, 'mlu'=>$mlu, 'seqfacility'=>$sequencing_facility}});
    }
    
}
    

sub date {
    return (split ' ', $_[0])[0]
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

sub htmlencode {
    my $str = shift;

    $str =~ s/Å/&Aring;/g;
    $str =~ s/å/&aring;/g;
    $str =~ s/Ä/&Auml;/g;
    $str =~ s/ä/&auml;/g;
    $str =~ s/Ö/&Ouml;/g;
    $str =~ s/ö/&ouml;/g;

    return $str;
}
