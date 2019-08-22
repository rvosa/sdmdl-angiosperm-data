#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use JSON;
use File::Spec;
use Getopt::Long;
use LWP::UserAgent;
use Text::CSV 'csv';
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $infile;
my $outdir;
my $verbosity = WARN;
GetOptions(
	'infile=s' => \$infile,
	'outdir=s' => \$outdir,
	'verbose+' => \$verbosity,
);

# configure services
Bio::Phylo::Util::Logger->VERBOSE( '-level' => $verbosity, '-class' => 'main' );
my ( $v, $d, $f ) = File::Spec->splitpath( $infile );
my $dbh = DBI->connect( "dbi:CSV:f_dir=$d;csv_sep_char=\t" );

# query the taxon list
my $sth = $dbh->prepare("select allmb_name,gbif_species_key,family from $f");
$sth->execute();
while( my $row = $sth->fetchrow_hashref() ) {

	# fetch the occurrences
	INFO "Going to fetch occurrences for " . $row->{'allmb_name'};
	my $dir  = $outdir . '/' . $row->{'family'};	
	my $file =    $dir . '/' . $row->{'allmb_name'} . '.csv';	
	if ( not -e $file ) {
		my @occ = occurrences_for_species( $row );
		if ( @occ ) {
			INFO "Going to write ".scalar(@occ)." occurrences to file";
			mkdir $dir if not -d $dir;
			csv(
				'in'  => \@occ,
				'out' => $file,
			);
		}
		else {
			WARN "Fetched 0 occurrences";
		}
	}	
}

sub occurrences_for_species {
	my $row = shift;
	my %def = (
		'hasCoordinate' => 'true',
		'limit'         => 300,
		'basisOfRecord' => 'PRESERVED_SPECIMEN',
		'year'          => '1900,2019',
		'offset'        => '%s',
		'taxonKey'      => $row->{'gbif_species_key'},
	);
	my $api = 'http://api.gbif.org/v1/occurrence/search?';
	my $url = $api . join '&', map { $_ . '=' . $def{$_} } keys %def;		
	my $off = 0;
	my ( @result, %seen );
	FETCH: while(1) {
		DEBUG "Offset: $off";
		my $ua = LWP::UserAgent->new;
		my $response = $ua->get( sprintf( $url, $off ) );
		if ( $response->is_success ) {
			my $obj;
			eval { $obj = decode_json( $response->decoded_content ) };
			if ( $@ ) {
				ERROR "Problem decoding JSON: $@";
				last FETCH;
			}
			for my $r ( @{ $obj->{'results'} } ) {
				my %h = %$r;
				my ( $lat, $long, $key ) = @h{qw[decimalLatitude decimalLongitude key]};
				if ( $lat and $long and not $seen{"$lat/$long"}++ ) {			
					push @result, {
						'gbif_id'           => $key,
						'decimal_longitude' => $long,
						'decimal_latitude'  => $lat,
						'taxon_name'        => $row->{'allmb_name'},
					}
				}
			}
			last FETCH if $obj->{'endOfRecords'};
			$off += 300;
		}
		else {
			ERROR "Problem fetching $url: " . $response->status_line;
			last FETCH;
		}
	}
	return @result;
}