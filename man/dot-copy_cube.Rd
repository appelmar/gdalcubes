% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cube.R
\name{.copy_cube}
\alias{.copy_cube}
\title{Create a data cube proxy object copy}
\usage{
.copy_cube(cube)
}
\arguments{
\item{cube}{source data cube proxy object}
}
\value{
copied data cube proxy object
}
\description{
Copy a data cube proxy object without copying any data
}
\details{
This internal function copies the complete processing chain / graph of a data cube but does not copy any data
It is used internally to avoid in-place modification for operations with potential side effects on source data cubes.
}
