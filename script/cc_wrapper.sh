#!/bin/bash
csvs=`find ../data/occurrences -name "*.csv"`
for csv in $csvs; do
	./coordinate_cleaner.R -i $csv -o $csv 2> $csv.log
done