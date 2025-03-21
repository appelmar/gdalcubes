% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdalcubes.R
\docType{package}
\name{gdalcubes}
\alias{gdalcubes-package}
\alias{gdalcubes}
\title{gdalcubes: Earth Observation Data Cubes from Satellite Image Collections}
\description{
Processing collections of Earth observation images as on-demand multispectral, multitemporal raster data cubes. Users
define cubes by spatiotemporal extent, resolution, and spatial reference system and let 'gdalcubes' automatically apply cropping, reprojection, and 
resampling using the 'Geospatial Data Abstraction Library' ('GDAL'). Implemented functions on data cubes include reduction over space and time, 
applying arithmetic expressions on pixel band values, moving window aggregates over time, filtering by space, time, bands, and predicates on pixel values, 
exporting data cubes as 'netCDF' or 'GeoTIFF' files, plotting, and extraction from spatial and or spatiotemporal features.  
All computational parts are implemented in C++, linking to the 'GDAL', 'netCDF', 'CURL', and 'SQLite' libraries. 
See Appel and Pebesma (2019) <doi:10.3390/data4030092> for further details.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/appelmar/gdalcubes}
  \item Report bugs at \url{https://github.com/appelmar/gdalcubes/issues/}
}

}
\author{
\strong{Maintainer}: Marius Appel \email{marius.appel@hs-bochum.de} (\href{https://orcid.org/0000-0001-5281-3896}{ORCID})

Other contributors:
\itemize{
  \item Edzer Pebesma (\href{https://orcid.org/0000-0001-8049-7069}{ORCID}) [contributor]
  \item Roger Bivand [contributor]
  \item Jeroen Ooms (\href{https://orcid.org/0000-0002-4035-0289}{ORCID}) [contributor]
  \item Lewis Van Winkle [copyright holder]
  \item Ole Christian Eidheim [copyright holder]
  \item Howard Hinnant [copyright holder]
  \item Adrian Colomitchi [copyright holder]
  \item Florian Dang [copyright holder]
  \item Paul Thompson [copyright holder]
  \item Tomasz Kamiński [copyright holder]
  \item Dropbox, Inc. [copyright holder]
}

}
