#!/usr/bin/env Rscript
library(adehabitatHS, quietly = T)
library(raster, quietly = T)
library(SDMTools, quietly = T)
library(factoextra, quietly = T)
library(ecospat, quietly = T)
library(cluster, quietly = T)
library(ape, quietly=T)
library(adehabitatMA, quietly=T)
library(optparse)
library(yaml)

# process command line arguments
taxa.names.file # list of taxon names, one per line
REPO_HOME
tiff.dir # non-bioclim TIFF files, e.g. root/5_deg in Ungulates
bioclim.dir # where to download bioclim TIFF files
bioclim.res <- 5 # results in root/wc5/ in Ungulates
niche.file # output file with niche traits

###### GETOPTIONS #######

# Import the taxa list 
taxa.names <- scan(
    taxa.names.file, 
    sep = "\n", 
    what = character() 
)

# Download bioclim data
gis.layers <- raster::getData(
    "worldclim",
    var = "bio",
    res = 5,
    path = bioclim.dir,
    download = T
)

# Turn the file names into layer names: strip the prefix (which might include
# the resolution) and strip the file extension
gis.layers.names <- list.files(tiff.dir)
gis.layers.names <- gsub('current_5arcmin_','',gis.layers.names)
gis.layers.names <- gsub('.tif','',gis.layers.names)

# Combine the layer names with those we've already read from BIOCLIM
gis.layers.names <- c( names(gis.layers), gis.layers.names )

# Iterate over non-bioclim TIFF files
for (file.name in list.files(tiff.dir)) {
    
    # Stack with previously read layers
    gis.layers <- stack(
        gis.layers,
        
        # Read as raster
        raster(paste(tiff.dir, file.name, sep = ""))
    )
}

# Apply all names
names(gis.layers) <- gis.layers.names

# Make SpatialPixelsDataFrame with CRS string from GIS layers
crs.string <- "+proj=longlat +datum=WGS84"
gis.layers.spdf <- as(gis.layers, "SpatialPixelsDataFrame")
sp::proj4string(gis.layers.spdf) <- CRS(crs.string)

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
    proj4string = new( "CRS", projargs = crs.string )
) 

# Populate the empty dataframe with lat/lon values from the taxa list
# set for example i in 1:100 in the case of memory limits 
for ( i in 1:length(taxa.names)) {

	# Read occurrences as matrix
    csv.file <- sprintf('%s/data/filtered/%s.csv', REPO_HOME, taxa.names[i])
    occurrence.matrix <- as.matrix(read.csv(csv.file))

	# Convert matrix to coordinate data frame
    points.df <- as.data.frame(cbind(occurrence.matrix[, c("decimal_longitude", "decimal_latitude")]))
    points.df$decimal_longitude <- as.numeric(as.character(points.df$decimal_longitude))
    points.df$decimal_latitude <- as.numeric(as.character(points.df$decimal_latitude))
    coordinates(points.df) <- ~ decimal_longitude + decimal_latitude
    
    # Populate SpatialPointsDataFrame for focal taxon
    focal.spdf <- SpatialPointsDataFrame(
    	points.df, 
    	data.frame(species = rep(taxa.names[i], NROW(occurrence.matrix))), 
    	proj4string = CRS
    )
    proj4string(focal.spdf) <- CRS(crs.string)
    
    # Append to cumulative data frame
    occurrences.spdf <- rbind(occurrences.spdf, focal.spdf)
}

# data frame of values in the SpatialPixelsDataframe, i.e. the "traits"
traitvalues.df <- slot(gis.layers.spdf, "data")
traitvalues.df <- tibble::rowid_to_column(traitvalues.df, "ID")
traitvalues.df <- na.omit(traitvalues.df)

# data frame of number of occurrences per pixel 
npoints <- count.points(occurrences.spdf, gis.layers.spdf)
noccurrences.df <- slot(npoints, "data")
noccurrences.df <- tibble::rowid_to_column(noccurrences.df, "ID")
noccurrences.df <- subset(noccurrences.df, noccurrences.df$ID %in% traitvalues.df$ID)
noccurrences.df <- subset(noccurrences.df, select = -ID)

# PCA to collect and standardize environmental variables per occurrence
traitvalues.dudi <- ade4::dudi.pca(traitvalues.df[2:42], scannf = F)

# niche analysis to average the standardized values fromt the PCA per species
niche.dudi <- ade4::niche(traitvalues.dudi, noccurrences.df, scannf = F)

# write averages per species for standardized environmental variables to CSV
write.csv(niche.dudi$tab, niche.file)