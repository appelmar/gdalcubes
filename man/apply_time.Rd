% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_time.R
\name{apply_time}
\alias{apply_time}
\title{Apply a function over (multi-band) pixel time series}
\usage{
apply_time(x, ...)
}
\arguments{
\item{x}{input data}

\item{...}{additional arguments passed to method implementations}
}
\value{
return value and type depend on the class of x
}
\description{
This generic function applies a function on pixel time series of a data cube, an R array, or other classes if implemented.
The resulting object is expected to have the same spatial and temporal shape as the input, i.e., no reduction is performed.
}
\examples{
# 1. input is data cube
# create image collection from example Landsat data only 
# if not already done in other examples
if (!file.exists(file.path(tempdir(), "L8.db"))) {
  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
                         ".TIF", recursive = TRUE, full.names = TRUE)
  create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
}

L8.col = image_collection(file.path(tempdir(), "L8.db"))
v = cube_view(extent=list(left=388941.2, right=766552.4, 
              bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
              srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
L8.cube = raster_cube(L8.col, v) 
L8.cube = select_bands(L8.cube, c("B04", "B05")) 
L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI")

# Apply a user defined R function
apply_time(L8.ndvi, names="NDVI_residuals", 
   FUN=function(x) {
      y = x["NDVI",]
      if (sum(is.finite(y)) < 3) {
         return(rep(NA,ncol(x)))
      }
      t = 1:ncol(x)
      return(predict(lm(y ~ t)) -  x["NDVI",])})

# 2. input is array
d <- c(4,16,32,32)
x <- array(rnorm(prod(d)), d)
z <- apply_time(x, function(v) {
  y = matrix(NA, ncol=ncol(v), nrow=2)
  y[1,] = (v[1,] + v[2,]) / 2
  y[2,] = (v[3,] + v[4,]) / 2
  y
})
dim(z)

}
\seealso{
\code{\link{apply_time.cube}}

\code{\link{apply_time.array}}
}
