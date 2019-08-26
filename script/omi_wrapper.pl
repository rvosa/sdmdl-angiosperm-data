#!/usr/bin/perl
use strict;
use warnings;
for my $i ( 309 .. 1585 ) {
    system("./OMI.R -f $i -l $i -o ../data/niche_traits$i.csv -c ../config.yml");
}
