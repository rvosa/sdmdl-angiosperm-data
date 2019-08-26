#!/usr/bin/env Rscript
library(sp)
library(raster)
library(adehabitatMA)
library(optparse)
library(yaml)
library(dplyr)
library(ade4)
library(tibble)
library(logger)

# Process command line arguments, The usage is thus:
# ./OMI.R -c ../config.yml
opt <- optparse::parse_args( 
    OptionParser(
        option_list = list(
            optparse::make_option( 
                c("-c", "--config"), 
                type="character",
                default="config.yml"
            ),
            optparse::make_option(
                c("-o", "--outfile"),
                type="character",
                default="niche_traits.csv"
            ),
            optparse::make_option(
                c("-f", "--first"),
                type="integer",
                default=1
            ),
            optparse::make_option(
                c("-l", "--last"),
                type="integer",
                default=800
            )
        )
    ) 
)
logger::log_threshold(TRACE)
logger::log_info('Going to read config file {opt$config}')
config <- yaml::yaml.load_file(opt$config)

# Read and stack TIFF files
gis.layers  <- raster::stack()
layer.files <- config$gis_data$files
layer.proj  <- config$gis_data$proj
layer.datum <- config$gis_data$datum
logger::log_info('Going to stack {length(config$gis_data$files)} layers')
for ( layer.name in names(layer.files) ) {
    
    # Stack with previously read layers
    logger::log_debug('Adding layer {layer.name}')
    gis.layers <- raster::stack(
        gis.layers,
        
        # Read as raster
        raster::raster(layer.files[[layer.name]])
    )
}

# Apply user supplied layer names
names(gis.layers) <- names(layer.files)

# Make SpatialPixelsDataFrame with CRS string from GIS layers
logger::log_info('Going to coerce layers to SpatialPixelsDataFrame')
crs.obj <- CRS(
    sprintf("+proj=%s +datum=%s", layer.proj, layer.datum), 
    doCheckCRSArgs = T
)
gis.layers.spdf <- as(gis.layers, "SpatialPixelsDataFrame")
sp::proj4string(gis.layers.spdf) <- crs.obj

# Make SpatialPointsDataFrame from occurrences
# WARNING: the SpatialPointDataFrames are large files and therefore it might 
# not be possible to run all species. The datasets can be split by for example 
# changing 'i in length(taxa.names)' to 'i in 1:100' to calculate the 
# normalized mean values for the first 100 species. Splitting the datasets in 
# smaller subsets does not affect the outcomes because the normalization and 
# averaging is done per species.

# Create an empty SpatialPointDataFrame to populate in the following loop
occurrences.spdf <- new(
    "SpatialPointsDataFrame", 
    coords = structure(
        numeric(0), 
        .Dim = c(0L, 2L),
        .Dimnames = list( NULL, c("x", "y") )
    ),  
    bbox = structure(
        c(1,1,1,1), 
        .Dim = c(2L, 2L),                         
        .Dimnames = list( c("x","y"), c("min","max") )
    ),
    proj4string = crs.obj
) 

# Populate the empty dataframe with lat/lon values from the taxa list
# set for example i in 1:100 in the case of memory limits 
taxa.names <- names(config$occurrences)
logger::log_info('Going to read occurrences for taxa {opt$first}..{opt$last}')
for ( i in opt$first:opt$last ) {
    
    # Prevent out of bounds errors
    if ( i > length(taxa.names) ) {
        break
    }
    logger::log_debug('Going to read taxon {i}: {taxa.names[i]}')

	# Read occurrences
    csv.file <- config$occurrences[[taxa.names[i]]]
    points.df <- dplyr::select(
        read.csv(csv.file),
        decimal_longitude,
        decimal_latitude
    )
    sp::coordinates(points.df) <- ~ decimal_longitude + decimal_latitude
    
    # Populate SpatialPointsDataFrame for focal taxon
    focal.spdf <- SpatialPointsDataFrame(
    	points.df, 
    	data.frame(species = rep(taxa.names[i], NROW(points.df))), 
    	proj4string = crs.obj
    )
    proj4string(focal.spdf) <- crs.obj
    
    # Append to cumulative data frame
    occurrences.spdf <- rbind(occurrences.spdf, focal.spdf)
}

# data frame of values in the SpatialPixelsDataframe, i.e. the "traits"
logger::log_info('Going to make data frame of niche traits from GIS layers')
traitvalues.df <- slot(gis.layers.spdf, "data")
traitvalues.df <- tibble::rowid_to_column(traitvalues.df, "ID")
traitvalues.df <- na.omit(traitvalues.df)

# data frame of number of occurrences per pixel 
logger::log_info('Going to make data frame of occurrences per pixel')
npoints <- adehabitatMA::count.points(occurrences.spdf, gis.layers.spdf)
noccurrences.df <- slot(npoints, "data")
noccurrences.df <- tibble::rowid_to_column(noccurrences.df, "ID")
noccurrences.df <- subset(noccurrences.df, noccurrences.df$ID %in% traitvalues.df$ID)
noccurrences.df <- subset(noccurrences.df, select = -ID)

# PCA to collect and standardize environmental variables per occurrence
traitvalues.dudi <- ade4::dudi.pca(
    traitvalues.df[2:length(traitvalues.df)], 
    scannf = F
)

# niche analysis to average the standardized values fromt the PCA per species
niche.dudi <- ade4::niche(traitvalues.dudi, noccurrences.df, scannf = F)

# write averages per species for standardized environmental variables to CSV
rownames(niche.dudi$tab) <- gsub(
    x = rownames(niche.dudi$tab),
    pattern = "\\.",
    replacement = " "
)
write.csv( niche.dudi$tab, file = opt$outfile, append = T, quote = F )
