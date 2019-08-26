#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# process command line arguments
my $stem = '^niche_traits';
my $indir = '../data/';
GetOptions(
    'stem=s'  => \$stem,
    'indir=s' => \$indir,
);

# keep global @header across CSV files. merge records into %table keyed on taxon
my ( @header, %table, %i );

# scan $dir for files that match $stem
opendir my $dh, $indir or die $!;
while( my $entry = readdir $dh ) {
    
    # file matches
    if ( $entry =~ /$stem/ ) {
        
        # keep the 1-based index embedded in the file name
        if ( $entry =~ /(\d+)/ ) {
            my $i = $1;
            $i{$i}++;
        }
        
        # keep a local header to see if it matches the global
        my @lh;
        
        # start reading the file
        open my $in, '<', "${indir}/${entry}" or die $!;
        LINE: while(<$in>) {
            chomp;
            my @line = split /,/, $_;
            if ( not @lh ) {
            
                # first column has no label, as in R's row names
                @lh = @line[ 1 .. $#line ];
                
                # process the header line, first time to create global header
                if ( not @header ) {
                    @header = @lh;
                }
                
                # subsequent times to match local against global
                else {
                    if ( grep { $header[$_] ne $lh[$_] } 0 .. $#header ) {
                        warn "Header mismatch in ${indir}/${entry}";
                        die "Global:\n'@header'\nLocal:\n'@lh'\n";
                    }
                }
                next LINE;
            }
            
            # process the data record(s), first cell is row name (i.e. taxon)
            my $taxon  = shift @line;
            my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
            die "Seeing $taxon again in ${indir}/${entry}" if $table{$taxon};
            $table{$taxon} = \%record;
        }
    }
}

# check which ones are missing
my ($max) = sort { $b <=> $a } keys %i;
my ($min) = sort { $a <=> $b } keys %i;
my @missing = grep { ! $i{$_} } $min .. $max;
warn "Missing trait files for taxa:\n@missing";

# print the merged table as tab separated
print join(",", 'taxon', @header), "\n";
for my $taxon ( sort { $a cmp $b } keys %table ) {
    my @values = @{ $table{$taxon} }{ @header };
    print join(",", $taxon, @values), "\n";
}