#' Create an image collection from a STAC feature collection
#' 
#' This function creates an image collection from a STAC API collection response. It does not
#' need to read any image data. Additionally, bands can be filtered and asset links can be transformed to make them
#' readable for GDAL. 
#' 
#' @param s STAC feature collection
#' @param asset_names character vector with names of assets (e.g., bands) to be used, other assets will be ignored. By default (NULL), all asset names with "eo:bands" attributes will be used
#' @param asset_regex length 1 character defining a regular expression asset names must match to be considered
#' @param url_fun  optional function to modify URLs of assets, e.g, to add /vsicurl/ to URLS (the default)
#' @param out_file optional name of the output SQLite database file, defaults to a temporary file
#' @param property_filter optional function to filter STAC items (images) by their properties; see Details
#' @param skip_image_metadata logical, if TRUE per-image metadata (STAC item properties) will not be added to the image collection
#' @param srs character spatial reference system of images used either for images without corresponding STAC property ony or for all images
#' @param srs_overwrite logical, if FALSE, use srs only for images with unknown srs (missing STAC metadata)
#' @param duration character, if images represent time intervals, use either the"start" or "center" of time intervals
#' @note Currently, bbox results are expected to be WGS84 coordinates, even if bbox-crs is given in the STAC response.
#' @details 
#' 
#' The property_filter argument can be used to filter images by metadata such as cloud coverage. 
#' The functions receives all properties of a STAC item (image) as input list and is expected to produce a single logical value,
#' where an image will be ignored if the function returns FALSE.
#' 
#' Some STAC API endpoints may return items with duplicte IDs (image names), pointing to 
#' identical URLs. Such items are only added once during creation of the image collection.
#' 
#' @export
stac_image_collection <- function(s, out_file = tempfile(fileext = ".db"), 
                                  asset_names = NULL, asset_regex = NULL, 
                                  url_fun = .default_url_fun,
                                  property_filter = NULL, skip_image_metadata = FALSE,
                                  srs = NULL, srs_overwrite = FALSE, duration = c("center", "start")) {
  SUBBAND_SPLIT_CHAR = ":"
  if ("STACItemCollection" %in% class(s)) {
    s = s$features
  }
  if (!is.list(s)) {
    stop ("Input must be either a list or a STACItemCollection")
  }
  if (length(s) == 0) {
    stop ("Input does not include any STAC items")
  }
  
  if (file.exists(out_file)) {
    stop ("output file already exists")
  }
  
  
  # add band metadata if available
  bands = NULL
  asset_names_exist = rep(FALSE, length(asset_names))
  for (i in 1:length(s)) {
    for (j in 1:length(s[[i]]$assets)) {
      asset_name = names(s[[i]]$assets)[j]
      
      if (!is.null(asset_names)) { 
        if (asset_name %in% asset_names) {
          asset_names_exist[which(asset_names == asset_name)] = TRUE
        }
        else {
          next
        }
      }
      if (!is.null(asset_regex)) {
        if (!grepl(pattern = asset_regex, asset_name)) next
      }
      
      # assumption: use only assets with "eo:bands" metadata, unless its is part of asset_names
      if ("eo:bands" %in% names(s[[i]]$assets[[j]])) {
        nb = length(s[[i]]$assets[[j]][["eo:bands"]])
        if (nb == 1) {
          bands = union(bands, asset_name)
        }
        else {
          subnames = sapply(s[[i]]$assets[[j]][["eo:bands"]], function(x){x$name})
          joined_names = paste0(asset_name, SUBBAND_SPLIT_CHAR, subnames)
          bands = union(bands, joined_names)
        }
        next
      }

      # if no eo:bands metadata found, consider asset if type contains "image/*" and role equals "data"
      sroles = s[[i]]$assets[[j]]$roles
      stype = s[[i]]$assets[[j]]$type
      
      if (is.character(sroles) && is.character(stype)) {  
        if (grepl("image/", stype, fixed=TRUE) && "data" %in% sroles) {
          bands = union(bands, asset_name)
        }
      }
    }
  }
  
  # compare asset_names to found assets
  if (!is.null(asset_names)) { 
    a = which(!(asset_names %in% bands))
    if (length(a) > 0) {
      for (i in 1:length(a)) {
        if (asset_names_exist[a[i]]) {
          bands = union(bands, asset_names[a[i]])
          warning(paste0("STAC asset with name '", asset_names[a[i]] ,"' does not include eo:bands metadata and will be considered as a single band source"))
        }
        else {
          warning(paste0("STAC asset with name '", asset_names[a[i]] ,"' does not exist and will be ignored"))
        }    
      }
    }
  }
  
  
  bands_df = data.frame(id = 1:length(bands), name = bands, type = "", offset = NA,
                        scale = NA, unit = "", nodata = "", stringsAsFactors = FALSE)

  images_df = data.frame(id = integer(), name = character(), left = numeric(), 
                         top = numeric(), bottom = numeric(), right = numeric(), 
                         datetime = character(), proj = character())
  gdalrefs_df = data.frame(image_id = integer(), band_id = integer(), 
                           descriptor = character(),band_num = integer())
  image_md_df = data.frame(image_id = integer(),key = character(), value = character())
 
  
  for (i in 1:length(s)) {
    
    tryCatch(
      expr = {
          images_df_temp = data.frame(id = integer(), name = character(), left = numeric(), 
                                top = numeric(), bottom = numeric(), right = numeric(), 
                                datetime = character(), proj = character())
          gdalrefs_df_temp = data.frame(image_id = integer(), band_id = integer(), 
                                  descriptor = character(),band_num = integer())
          image_md_df_temp = data.frame(image_id = integer(),key = character(), value = character())
        if (!is.null(property_filter)) {
          if (!property_filter(s[[i]]$properties)) {
            next
          }
        }

        # bands
        img_has_bands = FALSE
        for (j in 1:nrow(bands_df)) {
          if (bands_df$name[j] %in% names(s[[i]]$assets)) {
            gdalrefs_df_temp = rbind(gdalrefs_df_temp, data.frame(image_id = i, band_id = j, 
                                    descriptor = url_fun(s[[i]]$assets[[bands_df$name[j]]]$href),
                                    band_num = 1, stringsAsFactors = FALSE))
            img_has_bands = TRUE
          }
          else {
            spl = strsplit(bands_df$name[j], SUBBAND_SPLIT_CHAR)[[1]]
            if (length(spl) == 2) {
              asset_name = spl[1]
              band_num = which(which(startsWith(bands_df$name, paste0(asset_name, SUBBAND_SPLIT_CHAR))) == j)
              
              if (length(band_num) == 1) {
                 gdalrefs_df_temp = rbind(gdalrefs_df_temp, data.frame(image_id = i, band_id = j, 
                                    descriptor = url_fun(s[[i]]$assets[[asset_name]]$href),
                                    band_num = band_num, stringsAsFactors = FALSE))
                img_has_bands = TRUE
              }
            }
          }
        }
        if (img_has_bands) {
          # fixes #60
          if (!is.null(srs) && srs_overwrite) {
            proj = srs
          }
          else {
            proj = s[[i]]$properties$"proj:epsg"
            if (!is.null(proj)) {
              if (!startsWith(toupper(proj), "EPSG:")) {
                proj = paste0("EPSG:", proj)
              }
            }
            if (is.null(proj)) {
              proj = s[[i]]$properties$"proj:wkt2"
            }
            if (is.null(proj)) {
              proj = s[[i]]$properties$"proj:projjson"
            }
            if (is.null(proj)) {
              if (!is.null(srs)) {
                proj = srs
              }
              else {
                stop(paste0("No projection info found in STAC item for image ", s[[i]]$id))
              }
            }
          }
          
          
          temp_datetime = NULL
          if (!is.null(s[[i]]$properties$datetime)) {
            temp_datetime = s[[i]]$properties$datetime
          }
          else if (!is.null(s[[i]]$properties$start_datetime)) {
            temp_starttime = s[[i]]$properties$start_datetime
            temp_endtime = s[[i]]$properties$end_datetime
            
            m = duration[1]
            if (m == "center") {
              if (!requireNamespace("lubridate", quietly = TRUE)) {
                stop("package lubridate required; please install first")
              }
              a = lubridate::ymd_hms(temp_starttime)
              b = lubridate::ymd_hms(temp_endtime)
              temp_datetime = format(a + (b-a)/2, "%Y-%m-%dT%H:%M:%S")
            }
            else {
              temp_datetime = temp_starttime
            }
          }
          else {
             stop(paste0("No datetime / start_datetime property found in STAC item for image ", s[[i]]$id))
          }
          
          # BBOX
          bbox = s[[i]]$bbox
          if (is.list(bbox)) {
            bbox = unlist(bbox)
          }
          
          # TO CHECK IN STAC SPEC, does BBOX always exist? Is it always WGS84? 
          # TODO: transform, if bbox-crs is given
          images_df_temp = data.frame(id = i, name = s[[i]]$id, left = bbox[1], top = bbox[4], bottom = bbox[2],
                                    right = bbox[3], datetime = temp_datetime, proj = proj, stringsAsFactors = FALSE)

          # Duplicate image name check
          if (s[[i]]$id %in% images_df$name) {
            if (.pkgenv$debug) {
              message(paste0("Skipping STAC item ", s[[i]]$id) , " due to duplicate id (image with identical name already exists)")
            }
            next
          }

          # image metadata
          if (!skip_image_metadata) {
            image_md_df_temp = data.frame(image_id = rep(i, length(s[[i]]$properties)),  
                                          key = names(s[[i]]$properties),
                                          value = as.character(s[[i]]$properties))
          }

          # add to result tables
          images_df = rbind(images_df, images_df_temp)
          gdalrefs_df = rbind(gdalrefs_df, gdalrefs_df_temp)
          if (!skip_image_metadata) {
            image_md_df = rbind(image_md_df, image_md_df_temp)
          }
          
        }
      },
      error = function(e){
        warning(paste0("Skipping STAC item ", i, " due to the following error: ", e))
      })
  }
  
  
  if (nrow(images_df) == 0) {
    stop("Collection does not contain any images")
  }
  
  
  if (skip_image_metadata) {
    gc_create_stac_collection(bands_df, images_df, gdalrefs_df, path.expand(out_file), data.frame())
  }
  else {
    gc_create_stac_collection(bands_df, images_df, gdalrefs_df, path.expand(out_file), image_md_df)
  }
  return(image_collection(out_file))
}






#' Default function to convert hrefs of STAC response aassets to GDAL dataset identifiers including VSI prefixes 
#' @param url a single URL
#' @examples 
#' .default_url_fun("s3://bucket/object/image.tif")
#' @noRd
.default_url_fun <- function(url) {
  
  if (startsWith(url, "s3://")) {
    return(paste0("/vsis3/", substr(url, 6, nchar(url))))
  }
  else if (startsWith(url, "gs://")) {
    return(paste0("/vsigs/", substr(url, 6, nchar(url))))
  }
  else if (startsWith(url, "http")) {
    return(paste0("/vsicurl/", url))
  }
  
  # default: try vsicurl
  return(paste0("/vsicurl/", url))
}


