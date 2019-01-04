library(gdalcubes)
x = gcbs_read_stream_as_array()

out <- reduce_time(x, function(x) {
    ndvi <- (x[8,]-x[4,])/(x[8,]+x[4,])
    if (all(is.na(x))) return(NA)
    xx = max(ndvi,na.rm=TRUE) - min(ndvi, na.rm=T)
    return(xx)
})

#out = (x[8,,,,drop=FALSE] - x[4,,,,drop=FALSE]) / (x[8,,,,drop=FALSE] + x[4,,,,drop=FALSE])

gcbs_write_stream_from_array(out)


