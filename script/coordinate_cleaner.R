#!/usr/bin/env Rscript
library(CoordinateCleaner)
library(optparse)
library(tidyverse)

# Process command line arguments, Python style. The usage is thus:
# ./coordinate_cleaner.R -i <infile> -o <outfile>
opt <- parse_args( 
    OptionParser(
        option_list = list(
            make_option( c("-i", "--input"), type="character" ),
            make_option( c("-o", "--output"), type="character" )
        )
    ) 
)

# Read an input data file as CSV. These files meet the default expectations
# of read.csv in that they have a header and use the ',' as a field separator.
# DarwinCore uses CamelCase, which we had decided to change to underscore_case
# because that would allow us to import the data into SQLite if we ever needed
# to (which is case-insensitive for column names). However, CoordinateCleaner
# instead just flattens the CamelCase, so we have to do some column renaming
# here.
input <- select(
    read.csv(opt$input),
    gbif_id = gbif_id,
    decimallatitude = decimal_latitude,
    decimallongitude = decimal_longitude,
    species = taxon_name
)
input$dataset <- opt$input

# Run the CoordinateCleaner functions piped together. These functions have
# expectations for the default names of the columns in the data frame. We
# met these expectations by selecting and renaming the columns when we
# parsed the CSV.
clean <- input %>%
    cc_dupl()  %>%    # duplicates    
    cc_val()   %>%    # invalid lat/lon coordinates
    cc_equ()   %>%    # identical lat/lon
    cd_ddmm()  %>%    # degree conversion error   
    cd_round() %>%    # rasterized coordinates    
    cc_cap()   %>%    # coordinates near country capitals
    cc_cen()   %>%    # coordinates near country/province centroids
    cc_gbif()  %>%    # coordinates near GBIF headquarters
    cc_inst()  %>%    # coordinates near biodiversity institutions
    cc_sea()   %>%    # at sea
    cc_outl()         # statistical outliers

# Write the output:
write.csv(
    select(
        clean,
        gbif_id,
        taxon_name = species,
        decimal_longitude = decimallongitude,
        decimal_latitude = decimallatitude
    ), 
    file = opt$output, 
    quote = FALSE, 
    row.names = FALSE
)
