#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Text::CSV 'csv';
use Bio::DB::Taxonomy;
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $infile;
my $outfile;
my $verbosity = WARN;
GetOptions(
	'infile=s'  => \$infile,
	'outfile=s' => \$outfile,
	'verbose+'  => \$verbosity,
);

# configure services
Bio::Phylo::Util::Logger->VERBOSE( '-level' => $verbosity, '-class' => 'main' );
my $entrez = Bio::DB::Taxonomy->new( '-source' => 'entrez' );

# start reading the file
INFO "Going to read input file $infile";
my $data = csv(
	'in'       => $infile,
	'sep_char' => "\t",
	'headers'  => 'auto',
);

# store family for genus
my %lookup;
for my $r ( @$data ) {
	my $name = $r->{'allmb_name'};
	my ($genus) = split /_/, $name;
	
	# have not seen this name yet
	if ( not $lookup{$genus} ) {
		INFO "Going to query Entrez for $name";
		eval {
			if ( my $taxonid = $entrez->get_taxonid($name) ) {
				my $taxon = $entrez->get_taxon( '-taxonid' => $taxonid );
				my $name  = $taxon->scientific_name;
				my $rank  = $taxon->rank;
				while( $rank ne 'family' ) {
					$taxon = $taxon->ancestor;
					$name  = $taxon->scientific_name;
					$rank  = $taxon->rank;
					DEBUG "$name => $rank";
				}
				$lookup{$genus} = $name;
				$r->{$rank}     = $name;
			}	
		};
		if ( $@ ) {
			ERROR "Error querying Entrez: $@";
		}
	}		
	else {
		my $family = $lookup{$genus};
		$r->{'family'} = $family;
		DEBUG "Already seen $genus (genus), $family (family)";
	}
}

# write the output file
INFO "Going to write output file $outfile";
csv(
	'in'       => $data,
	'out'      => $outfile,
	'sep_char' => "\t",
);

