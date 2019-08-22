#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use File::Spec;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $infile;
my $db;
my $state;
my $verbosity = WARN;
GetOptions(
	'infile=s' => \$infile,
	'db=s'     => \$db,
	'state=s'  => \$state,
	'verbose+' => \$verbosity,
);

# configure services
Bio::Phylo::Util::Logger->VERBOSE( '-level' => $verbosity, '-class' => 'main' );
my ( $v, $d, $f ) = File::Spec->splitpath( $infile );
my $dbh = DBI->connect( "dbi:CSV:f_dir=$d;csv_sep_char=\t" );
my $fac = Bio::Phylo::Factory->new;

# set up nexus scaffolding
my $proj   = $fac->create_project;
my $taxa   = $fac->create_taxa;
my $data   = $fac->create_matrix( '-type' => 'standard', '-taxa' => $taxa );
my $forest = $fac->create_forest( '-taxa' => $taxa );
$proj
	->insert($taxa)
	->insert($data)
	->insert($forest);

# build taxa block and character state matrix 
my $query = "select * from $f where $state is not null order by allmb_name";
INFO "Going to run query '$query'";
my $sth = $dbh->prepare($query);
$sth->execute();
my ( $counter, %label, %family ) = ( 0 );
while( my $row = $sth->fetchrow_hashref() ) {
	
	# create taxon and data record
	my $char  = $label{ $row->{$state} } // ( $label{ $row->{$state} } = $counter++ );
	my $taxon = $fac->create_taxon( 
		'-name'    => $row->{allmb_name},
		'-link'    => 'http://gbif.org/species/' . $row->{gbif_species_key},
	);
	my $datum = $fac->create_datum(
		'-type'  => 'standard',
		'-char'  => [ $char ],
		'-taxon' => $taxon,
		'-name'  => $row->{allmb_name},
	);
	$taxa->insert($taxon);
	$data->insert($datum);
	
	# add taxon to family
	my $famname = $row->{family};
	my $fam = $family{$famname} || [];
	push @$fam, $taxon;
	$family{$famname} = $fam;
	
	DEBUG sprintf("Created taxon %s with state %s", $taxon->get_name, $char);
}
$data->set_statelabels( 
	[ 
		[ 
			 map { "'$_'" } 
			sort { $label{$a} <=> $label{$b} } 
			keys %label 		
		] 	
	] 
);

# build tree
INFO "Going to extract subtree from $db";
my $names = join ',', map { $_->get_name } @{ $taxa->get_entities };
my $tree  = parse_tree(
	'-string' => `megatree-pruner -d $db -l $names`,
	'-format' => 'newick',
);
for my $tip ( @{ $tree->get_terminals } ) {
	$tip->set_taxon( $taxa->get_by_name( $tip->get_name ) );
}
$forest->insert($tree->set_name('https://doi.org/10.1002/ajb2.1019 - ALLMB.tre@v0.1'));

$_->set_name('') for @{ $tree->get_internals };

# label nodes by family
for my $famname ( keys %family ) {
	my @tips = map { $_->get_nodes->[0] } @{ $family{$famname} };
	no warnings 'recursion';
	my $mrca = $tree->get_mrca(\@tips);
	$mrca->set_name($famname) if $mrca->is_internal;
}

# write nexus
INFO "Going to write nexus";
print $proj->to_nexus( '-statelabels' => 1, '-nodelabels' => 1 );