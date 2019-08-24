#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# process command line arguments
my $precision = 1;
my $infile;
GetOptions(
    'precision=i' => \$precision,
    'infile=s'    => \$infile,
);

# make sprintf template with floating precision
my $tmpl = '%.' . $precision . 'f';

# start reading the file
my ( $filtered, @header, @cleaned ) = ( 0 );
open my $in, '<', $infile or die $!;
LINE: while(<$in>) {
    chomp;
    my @line = split ',', $_;
    if ( not @header ) {
        @header = @line;
        next LINE;
    }
    my %record = map { $header[$_] => $line[$_] } 0 .. $#header;
    my ( $lat, $lon ) = @record{'decimal_latitude', 'decimal_longitude'};
    if ( sprintf($tmpl, $lat) != $lat and sprintf($tmpl, $lon) != $lon ) {
        push @cleaned, \%record;
    }
    else {
        $filtered++;
    }
}
close $in;

# replace, if we have to
if ( $filtered ) {
    open my $out, '>', $infile or die $!;
    print $out join(',',@header), "\n";
    for my $r ( @cleaned ) {
        my %h = %$r;
        my @v = @h{@header};
        print $out join(',',@v), "\n";
    }
    close $out;
    print "git commit -m 'Removed $filtered low-precision records' $infile\n";
}