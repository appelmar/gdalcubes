% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reduce.R
\name{reduce_space}
\alias{reduce_space}
\title{Reduce multidimensional data over space}
\usage{
reduce_space(x, ...)
}
\arguments{
\item{x}{object to be reduced}

\item{...}{further arguments passed to specific implementations}
}
\value{
return value and type depend on the class of x
}
\description{
This generic function applies a reducer function over a data cube, an R array, or other classes if implemented.
}
\examples{
# create image collection from example Landsat data only 
# if not already done in other examples
if (!file.exists(file.path(tempdir(), "L8.db"))) {
  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
                         ".TIF", recursive = TRUE, full.names = TRUE)
  create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
}

L8.col = image_collection(file.path(tempdir(), "L8.db"))
v = cube_view(extent=list(left=388941.2, right=766552.4, 
              bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
              srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
reduce_space(raster_cube(L8.col, v) , "median(B02)")  


d <- c(4,16,32,32)
x <- array(rnorm(prod(d)), d)
y <- reduce_space(x, function(v) {
  apply(v, 1, mean)
})
 
}
\seealso{
\code{\link{reduce_space.cube}}

\code{\link{reduce_space.array}}
}
