library(gdalcubes)

setwd("/media/marius/Samsung_T5/eodata/Sentinel2_L2A/Bialowieza/")

list.files(pattern="*.zip", full.names = TRUE)


gcbs_create_image_collection(list.files(pattern="*.zip", full.names = TRUE), 
                             format="Sentinel2_L2A", unroll_archives=TRUE, out_file="S2_timeseries.db") 

s2.col = gcbs_image_collection("S2_timeseries.db")

v = gcbs_view(gcbs_cube(s2.col))
v = gcbs_view(view=v, dt="P1M", nx=1000, ny=round(0.5477928*1000))

gcbs_set_threads(8)
system.time(gcbs_eval(gcbs_reduce(gcbs_select_bands(gcbs_cube(s2.col, v), c("B02", "B03", "B04")), "median"), "test.nc"))
system.time(plot(gcbs_reduce(gcbs_select_bands(gcbs_cube(s2.col, v), c("B02", "B03", "B04")), "median"), rgb=3:1, zlim=c(0,1800)))



grep

