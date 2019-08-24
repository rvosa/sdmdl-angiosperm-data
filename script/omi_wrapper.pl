#!/usr/bin/perl
use strict;
use warnings;
for ( my $i = 1; $i <= 1600; $i += 100 ) {
    my $j = $i + 99;
    system("./OMI.R -f $i -l $j -o ../data/niche_traits.csv -c ../config.yml");
}