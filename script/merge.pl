#!/usr/bin/perl
use strict;
use warnings;

my @files = @ARGV;
my %records;
my %keys;
for my $f ( @files ) {
	my @header;
	open my $in, '<', $f or die $!;
	while(<$in>) {
		chomp;
		my @line = split /\t/, $_;
		if ( not @header ) {
			@header = @line;
		}
		else {
			my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
			my $id = $record{gbif_species_key};
			my $merged = $records{$id} || {};
			for my $k ( keys %record ) {
				$merged->{$k} = $record{$k};
				$keys{$k}++;
			}
			$records{$id} = $merged;
		}	
	}
}

my @keys = sort { $a cmp $b } keys %keys;

print join("\t",@keys), "\n";

for my $r ( sort { $a->{allmb_name} cmp $b->{allmb_name} } values %records ) {
	my %record = %$r;
	my @values = @record{@keys};
	no warnings 'uninitialized';
	print join("\t",@values), "\n";	
}