
# Tutorial 1: Landsat 8 image collection for pre- and post- earthquake mapping

TODO: Docker demo image 

In this tutorial we use a small  (approx. 24 GB) Landsat 8 surface reflectance dataset including 20 images. The images captured the Minahasa Peninsula 
in Indonesia in September and October 2018, before and after the earthquake and tsunami.


The dataset can be downloaded [here](https://uni-muenster.sciebo.de/s/6OjnEyxzt4rk6px/download). After extraction of the ZIP file, you will see the atmospherically corrected images in subdirectories as shown below.

```
./LC081140612018091601T1-SC20181019060231
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_band7.tif
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_ANG.txt
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1.xml
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_band4.tif
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_band3.tif
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_MTL.txt
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_band5.tif
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_band2.tif
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_aerosol.tif
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_band6.tif
./LC081140612018091601T1-SC20181019060231/LC08_L1TP_114061_20180916_20180928_01_T1_sr_band1.tif
./LC081150612018100901RT-SC20181019053055
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_MTL.txt
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_band1.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_band5.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_band6.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT.xml
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_aerosol.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_band7.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_band3.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_band2.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_sr_band4.tif
./LC081150612018100901RT-SC20181019053055/LC08_L1TP_115061_20181009_20181009_01_RT_ANG.txt
./LC081150592018092301T1-SC20181019053124
...
```


At first we look at the file structure of the images. For Landsat, this is rather straightforward as each image is a directory that simply contains 
one GeoTIFF image per band.


## Create an image collection
To make this structure understandable for gdalcubes, we have to define a JSON collection format description, which afterwards will be used in `gdalcubes create_collection`.
In this description we must define how files relate to bands, how files are assigned to an image based on their path, and how to extract the
recording datetime.
 
**1.** Define the bands for the image collection and how filenames refer to specific bands

In this case, GeoTIFF filenames contain a string like "_band1.tif" to identify the band. We formulate this in our collection format as 


**2.** Define which files belong to the same image

In this case, this is straightforward, files of one image all have the same scene identifier in its path.


**3.** Define how to extract date / time from filenames

Looking at filenames such as `LC08_L1TP_113059_20181011_20181011_01_RT_sr_band1.tif`, we define a regular expression to extract
`LC08_L1TP_113059_**20181011**_20181011_01_RT_sr_band1.tif`. The expression basically matches all filenames starting with LC08_L1TP_ + 5 arbitrary characters + _ followed by the 
date that is extracted.


**4.** Add a global filename filter

Here, we ignore everything that is not a ".tif" file (see `pattern`).


Below, you can find the complete JSON format description. We will store this in a new file `L08SR.json`


```
{
  "pattern" : ".+\\.tif",
  "images" : {
    "pattern" : ".*/(.+)/.+\\.tif"
  },
  "datetime" : {
    "pattern" : ".*LC08_L1TP_.{6}_(.+)_.*\\.tif", 
    "format" : "%Y%m%d"
  },
  "bands": {
    "B01" : {
       "pattern" : ".+_sr_band1\\.tif"
    },
    "B02" : {
       "pattern" : ".+_sr_band2\\.tif"
    },
    "B03" : {
       "pattern" : ".+_sr_band3\\.tif"
    },
    "B04" : {
       "pattern" : ".+_sr_band4\\.tif"
    },
    "B05" : {
       "pattern" : ".+_sr_band5\\.tif"
    },
    "B06" : {
       "pattern" : ".+_sr_band6\\.tif"
    },
    "B07" : {
      "pattern" : ".+_sr_band7\\.tif"
    },
    "B_AEROSOL" : {
        "pattern" : ".+_sr_aerosol\\.tif"
    }
  }
}
```



```
gdalcubes create_collection  -R -f L08SR.json . L08.db
> IMAGE COLLECTION 'L08.db' has 20 images with 8 bands from 160 GDAL dataset references
```



```
gdalcubes info L08.db
> IMAGE COLLECTION 'L08.db' has 20 images with 8 bands from 160 GDAL dataset references
> DIMENSIONS: 
>   BANDS:       (B01) (B02) (B03) (B04) (B05) (B06) (B07) (B_AEROSOL)
>   DATETIME:    (2018-09-16T00:00:00 - 2018-10-18T00:00:00)
>   Y / LAT:     (-2.50226 - 2.4942)
>   X / LON:     (117.952 - 123.716)
```




## Derive a low resolution mosaic image


Define view:

```
{
  "aggregation" : "min",
  "resampling" : "bilinear",
  "space" :
  {
    "left" : 117.952,
    "right" : 123.716,
    "top" :  2.4942,
    "bottom" : -2.50226,
    "proj" : "EPSG:4326",
    "nx" : 500
  },
  "time" :
  {
    "t0" : "2018-09-16",
    "t1" : "2018-10-18",
    "dt" : "P7D"
  }
}
```


```
gdalcubes reduce -r median -t 2 -v view.json L08.db mosaic.tif
```


![](tutorial_mosaic.png)



```
{
  "aggregation" : "min",
  "resampling" : "bilinear",
  "space" :
  {
    "left" : 119.6,
    "right" : 120.2,
    "top" :  -0.6,
    "bottom" : -1.1,
    "proj" : "EPSG:4326",
    "nx" : 2500
  },
  "time" :
  {
    "t0" : "2018-09-16",
    "t1" : "2018-09-28",
    "dt" : "P7D"
  }
}
{
  "aggregation" : "min",
  "resampling" : "bilinear",
  "space" :
  {
    "left" : 119.6,
    "right" : 120.2,
    "top" :  -0.6,
    "bottom" : -1.1,
    "proj" : "EPSG:4326",
    "nx" : 2500
  },
  "time" :
  {
    "t0" : "2018-09-28",
    "t1" : "2018-10-18",
    "dt" : "P7D"
  }
}
```


```
gdalcubes reduce -r median -t 2 -v view_pre.json L08.db pre.tif
gdalcubes reduce -r median -t 2 -v view_post.json L08.db post.tif
```

![](tutorial_pre.png)
![](tutorial_post.png)


## Stream data into R





```
library(gdalcubes)
x = read_stream_as_array() # get input chunk as four dimensional array

out <- reduce_time_multiband(x, function(x) {
    clouds <- which(x[8,] == 8)
    x[5,clouds] <- NA
    x[4,clouds] <- NA
    ndvi <- (x[5,]-x[4,])/(x[5,]+x[4,])
    if (all(!is.finite(ndvi))) return(NA)
    if (sum(!is.na(ndvi)) == 1) return(NA)
    return(max(ndvi,na.rm=T) - min(ndvi, na.rm=T))
})

write_stream_from_array(out)
```
 
```
gdalcubes stream --exec="Rscript --vanilla stream.R" -r median -v "view.json" -t 2 -c "16 256 256" L08.db stream.tif
```



# Tutorial 2: Analysing the global CHIRPS precipitation dataset in R


```
gdalcubes create_collection -R -f collection_format.json . CHIRPS.db
> IMAGE COLLECTION 'CHIRPS.db' has 13787 images with 1 bands from 13787 GDAL dataset references
```

```
ls -lh CHIRPS.db
> -rw-r--r-- 1 root root 3.9M Oct 31 09:44 CHIRPS.db
```

```
gdalcubes info CHIRPS.db
> IMAGE COLLECTION 'CHIRPS.db' has 13787 images with 1 bands from 13787 GDAL dataset references
> DIMENSIONS: 
>   BANDS:       (precipitation)
>   DATETIME:    (1981-01-01T00:00:00 - 2018-09-30T00:00:00)
>   Y / LAT:     (-50 - 50)
>   X / LON:     (-180 - 180)
```


What is the maximum daily precipitation in 2015?

view_max_2015.json
```
{
  "aggregation" : "none",
  "resampling" : "bilinear",
  "space" :
  {
    "left" : -180,
    "right" : 180,
    "top" : 50,
    "bottom" : -50,
    "proj" : "EPSG:4326",
    "nx" : 1800,
    "ny" : 500
  },
  "time" :
  {
    "t0" : "2015-01-01",
    "t1" : "2015-01-31",
    "dt" : "P1D"
  }
}
```


```
gdalcubes reduce -r "max" -v view_max_2015.json CHIRPS.db max_2015.tif
```


https://uni-muenster.sciebo.de/s/akVSnlT2z6dFiZX/download

![](tutorial_chirps_all_max.png)




