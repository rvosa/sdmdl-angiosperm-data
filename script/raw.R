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
                default="../config.yml"
            ),
            optparse::make_option(
                c("-o", "--outfile"),
                type="character",
                default="niche_traits.csv"
            )
        )
    ) 
)
logger::log_threshold(TRACE)
logger::log_info('Going to read config file {opt$config}')

# If the default fails it's because the wd is not set to the source file
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

# There will likely be a warning message here. This is because we reassign
# a new CRS (e.g. WGS84/longlat) without reprojecting. We can only get away with 
# that if we are sure that the GIS layers match what we say in the config.yml(!)
sp::proj4string(gis.layers.spdf) <- crs.obj

# List of taxon names as parsed out of the YAML
taxa.names <- names(config$occurrences)

# Empty data frame of taxa x GIS layers
rawmeans.df <- data.frame(taxa.names, row.names = 1)
for (n in names(gis.layers)) rawmeans.df[[n]] <- NA

# Start reading occurrences
logger::log_info('Going to read occurrences for taxa')
for ( i in 1:length(taxa.names) ) {
    logger::log_debug('Going to read taxon {i}: {taxa.names[i]}')
    
    # Read occurrences, note how we select long/lat
    csv.file <- config$occurrences[[taxa.names[i]]]
    points.df <- dplyr::select(
        read.csv(csv.file),
        decimal_longitude,
        decimal_latitude
    )
    sp::coordinates(points.df) <- ~ decimal_longitude + decimal_latitude
    
    # Iterate over GIS layers
    for ( gl in names(gis.layers) ) {
        val <- mean(extract(gis.layers[[gl]],points.df), na.rm=T)
        rawmeans.df[ taxa.names[i], gl ] <- val
    }
}

write.csv( rawmeans.df, file = opt$outfile, append = T, quote = F )
