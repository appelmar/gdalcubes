
# Mockup ideas for gdalcubes R package

```
library(gdalcubes)


gcbs_create_collection_format(name, bandnames, pattern_bands, pattern_datetime, pattern_image)

a <- gcbs_create_image_collection(name, paths, format )
gcbs_write_image_collection(collection, filename)
a <- gcbs_open_image_collection(filename)
summary(a)
print(a)

a$name
a$bands <- data frame(name, type, offset, scale, nodata)
a$images <- data.frame(proj, extent, datetime)
a$gdalrefs <- data.frame(id, bands... )

v <- gcbs_view(proj, nx, ny, dx, dy, l, r, u, b, dt, t0, t1, tile_z, tile_y, tile_x)
class(v) <- c("view", "list")
print(v)


c <- gcbs_cube(a, view, chunking) # if view missing, derive from image collection extent and size of graphics device

summary(a)
print(a)

gcbs_set_swarm(endpoints)
gcbs_set_threads(n)
gcbs_version()


x <- gcbs_stream(a, f)
print(x) 

y <- gcbs_reduce(a, reducer="median")

gcbs_eval(y)

plot(y)
image(y)



.Options(gdalcubes_gdalcache...)

gcbs_server(port=1111, endpoint="gdalcubes/api", whitelist=NULL)   # use processx package
```
