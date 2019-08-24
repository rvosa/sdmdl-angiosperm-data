#!/usr/bin/perl
use strict;
use warnings;

my @status = split /\n/, `git status -s`;
LINE: for my $line ( @status ) {
	chomp($line);
	if ( $line =~ /^ M (.+)/ ) {
		my $file = $1;
		open my $in, '<', "${file}.log" or die $!;
		my ( $test, @message );
		while(<$in>) {
			chomp;
			if ( /^Testing (.+)$/ ) {
				$test = $1;
			}
			if ( /^(Flagged|Removed) (\d+)/ ) {
				my ( $action, $result ) = ( $1, $2 );
				if ( $result > 0 ) {
					push @message, "${action} ${result} ${test}";
				}
			}
			if ( /^Error/ ) {
				warn $file;
				next LINE;
			}
		}
		my $m = join '. ', @message;		
		system( "git commit -m '$m' $file" );
		unlink( "${file}.log" );
	}
}