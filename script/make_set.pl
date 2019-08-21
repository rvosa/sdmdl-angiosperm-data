#!/usr/bin/perl
use strict;
use warnings;
use My::Taxon;
use Getopt::Long;
use Bio::Phylo::Forest::DBTree;
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $datafile;
my $treefile;
my $column;  # column in $datafile that contains the states of interest
my @states;  # state(s) that need to be present for a shift to be of interest
my $mindata = 100; # number of occurrences (on either side) for shifts to be interesting
my $verbosity = WARN;
GetOptions(
	'datafile=s' => \$datafile,
	'treefile=s' => \$treefile,
	'column=s'   => \$column,
	'states=s'   => \@states,
	'verbose+'   => \$verbosity,
);
	
# setup services
Bio::Phylo::Util::Logger->VERBOSE( '-level' => $verbosity, '-class' => 'main' );
my $data = My::Taxon->connect( $datafile )->resultset( 'Taxon' );
my $tree = Bio::Phylo::Forest::DBTree->connect( $treefile);
my $node = $tree->resultset('Node');
my $root = $tree->get_root;
my $taxa = $data->search( { $column => \'IS NOT NULL' } );
my $n = $taxa->count;

# we have to do two passes, first to label the nodes:
INFO "Attaching states to $n nodes";
my $i = 1;
my %anc;
TAXON: while( my $taxon = $taxa->next ) {

	# attach taxon to tip
	my $state = $taxon->$column;
	my $tip   = $node->find($taxon->allmb_id);
	my $tipid = $tip->id;
	$anc{$tipid} = $tip;
	$tip->set_generic( $column => { $state => $taxon } );
	DEBUG "Attached $column => $state to tip $i/$n";
	$i++;
	
	# carry over to ancestors
	my $ancrs = $tip->get_ancestors_rs;
	while( my $anc = $ancrs->next ) {
		my $cache    = $anc->get_generic($column) || {};
		my $exemplar = $cache->{$state} || $taxon;	
		my $ancid    = $anc->id;
		$anc{$ancid} = $anc;			
		if ( ( $taxon->num_occurrences || 0 ) >= ( $exemplar->num_occurrences || 0 ) ) {
			$cache->{$state} = $taxon;
			$anc->set_generic( $column => $cache );	
		}
		else {
			next TAXON;
		}		
	}
}

# second we do a pre-order traversal to fetch the earliest monotypic node
INFO "Going to find synapomorphy exemplars by pre-order traversal";
my @header = ( qw(allmb_name gbif_species_key num_occurrences), $column );
print join("\t",@header), "\n";
my %seen;

for my $node ( sort { $a->left <=> $b->left } values %anc ) {
	my $nid = $node->id;
	next if $seen{$nid};

	# if there is only one key the subtended clade is monotypic
	my %cache = %{ $node->get_generic($column) };
	if ( scalar(keys(%cache)) == 1 and scalar(grep {defined} map {$cache{$_}} @states) ) {
		DEBUG "Found monotypic node of interest";
				
		my $parent = $node->get_parent;			
		my ($n) = values %cache;
		my @values = map { $n->$_ } @header;
		print join("\t",@values), "\n";		
		
		INNER: for my $s ( values %{ $parent->get_generic($column) } ) {
			next INNER if $s->id == $n->id;
			my @values = map { $s->$_ } @header;
			print join( "\t", @values ), "\n";		
		}
		
		my ( $l, $r ) = ( $parent->left, $parent->right );
		for my $desc ( grep { $_->left > $l && $_->right < $r } values %anc ) {
			my $did = $desc->id;
			$seen{$did}++;
			DEBUG $did;
		}
	}
}
	



