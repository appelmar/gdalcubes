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
#' @note Currently, bbox results are expected to be WGS84 coordinates, even if bbox-crs is given in the STAC response.
#' @note This function is experimental.
#' @details 
#' 
#' The property_filter argument can be used to filter images by metadata such as cloud coverage. 
#' The functions receives all properties of a STAC item (image) as input list and is expected to produce a single logical value,
#' where an image will be ignored if the function returns FALSE.
#' 
#' 
#' @export
stac_image_collection <- function(s, out_file = tempfile(fileext = ".sqlite"), 
                                  asset_names = NULL, asset_regex = NULL, 
                                  url_fun = function(x) {paste0("/vsicurl/", x)},
                                  property_filter = NULL, skip_image_metadata = FALSE) {
  SUBBAND_SPLIT_CHAR = ":"
  if (!is.list(s)) {
    stop ("Input must be a list")
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
  
  
  bands_df = data.frame(id = 1:length(bands),
                        name = bands,
                        type = "",
                        offset = NA,
                        scale = NA,
                        unit = "",
                        nodata = "",
                        stringsAsFactors = FALSE)
  
  gdalrefs_image_id = NULL
  gdalrefs_band_id = NULL
  gdalrefs_descriptor = NULL
  gdalrefs_band_num = NULL
  
  images_id = NULL
  images_name = NULL
  images_left = NULL
  images_top = NULL
  images_bottom = NULL
  images_right = NULL
  images_datetime = NULL
  images_proj = NULL
  
  image_md_image_id = NULL
  image_md_key = NULL
  image_md_value = NULL
  
  for (i in 1:length(s)) {
    
    if (!is.null(property_filter)) {
      if (!property_filter(s[[i]]$properties)) {
        next
      }
    }
    
    
    # bands
    img_has_bands = FALSE
    for (j in 1:nrow(bands_df)) {
      if (bands_df$name[j] %in% names(s[[i]]$assets)) {
        gdalrefs_image_id = c(gdalrefs_image_id, i)
        gdalrefs_band_id = c(gdalrefs_band_id, j)
        gdalrefs_descriptor = c(gdalrefs_descriptor, url_fun(s[[i]]$assets[[bands_df$name[j]]]$href))
        gdalrefs_band_num = c(gdalrefs_band_num, 1)
        img_has_bands = TRUE
      }
      else {
        spl = strsplit(bands_df$name[j], SUBBAND_SPLIT_CHAR)[[1]]
        if (length(spl) == 2) {
          asset_name = spl[1]
          band_num = which(which(startsWith(bands_df$name, paste0(asset_name, SUBBAND_SPLIT_CHAR))) == j)
          
          if (length(band_num) == 1) {
            gdalrefs_image_id = c(gdalrefs_image_id, i)
            gdalrefs_band_id = c(gdalrefs_band_id, j)
            gdalrefs_descriptor = c(gdalrefs_descriptor, url_fun(s[[i]]$assets[[asset_name]]$href))
            gdalrefs_band_num = c(gdalrefs_band_num, band_num)
            img_has_bands = TRUE
          }
        }
      }
    }
    if (img_has_bands) {
      images_id = c(images_id, i)
      images_name = c(images_name, s[[i]]$id)
      images_proj =  c(images_proj, paste0("EPSG:",s[[i]]$properties$"proj:epsg"))
      images_datetime = c(images_datetime, s[[i]]$properties$datetime)
      
      # BBOX
      bbox = s[[i]]$bbox
      # TO CHECK IN STAC SPEC, does BBOX always exist? Is it always WGS84? 
      # TODO: transform, if bbox-crs is given
      
      images_left = c(images_left, bbox[1])
      images_right = c(images_right, bbox[3])
      images_top = c(images_top, bbox[4])
      images_bottom = c(images_bottom, bbox[2])
      
      # image metadata
      if (!skip_image_metadata) {
        image_md_image_id = c(image_md_image_id, rep(i, length(s[[i]]$properties)))
        image_md_key = c(image_md_key, names(s[[i]]$properties))
        image_md_value = c(image_md_value, as.character(s[[i]]$properties))
      }
      
    }
    
  }
  
  
  
  images_df = data.frame(id = images_id,
                         name = images_name,
                         left = images_left,
                         top = images_top,
                         bottom = images_bottom,
                         right = images_right,
                         datetime = images_datetime,
                         proj = images_proj, stringsAsFactors = FALSE)
  
  if (nrow(images_df) == 0) {
    stop("Collection does not contain any images")
  }
  
  
  gdalrefs_df = data.frame(
    image_id = gdalrefs_image_id,
    band_id = gdalrefs_band_id,
    descriptor = gdalrefs_descriptor,
    band_num = gdalrefs_band_num,
    stringsAsFactors = FALSE
  )
  
  if (skip_image_metadata) {
    libgdalcubes_create_stac_collection(bands_df, images_df, gdalrefs_df, path.expand(out_file), data.frame())
  }
  else {
    image_md_df = data.frame(
      image_id = image_md_image_id,
      key = image_md_key,
      value = image_md_value,
      stringsAsFactors = FALSE
    )
    libgdalcubes_create_stac_collection(bands_df, images_df, gdalrefs_df, path.expand(out_file), image_md_df)
  }
  return(image_collection(out_file))
}

