#' Plot a gdalcubes data cube
#'
#' @param x a data cube proxy object (class cube)
#' @param y __not used__
#' @param nbreaks number of breaks, should be one more than the number of colors given
#' @param breaks actual breaks used to assign colors to values; if missing, the function subsamples values and uses equally sized intervals between min and max or zlim[0] and zlim[1] if defined
#' @param col color definition, can be a character vector with nbreaks - 1 elements or a function such as \code{heat.colors}
#' @param key.pos position for the legend, 1 (bottom), 2 (left), 3 (top), or 4 (right). If NULL (the default), do not plot a legend.
#' @param bands integer vector with band numbers to plot (this must be band numbers, not band names)
#' @param t integer vector with time indexes to plot (this must be time indexes, not date / time)
#' @param rgb bands used to assign RGB color channels, vector of length 3 (this must be band numbers, not band names)
#' @param zlim vector of length 2, defining the minimum and maximum values to either derive breaks, or define black and white values in RGB plots
#' @param gamma gamma correction value, used for RGB plots only
#' @param periods.in.title logical value, if TRUE, the title of plots includes the datetime period length as ISO 8601 string
#' @param join.timeseries logical, for pure time-series plots, shall time series of multiple bands be plotted in a single plot (with different colors)?
#' @param axes logical, if TRUE, plots include axes
#' @param ncol number of columns for arranging plots with  \code{layout()}, see Details 
#' @param nrow number of rows for arranging plots with  \code{layout()}, see Details
#' @param na.color color used to plot NA pixels
#' @param downsample length-one integer or logical value used to select only every i-th pixel (in space only) for faster plots; by default (TRUE), downsampling will be determined automatically based on the resolution of the graphics device; set to FALSE to avoid downsampling.
#' @param ... further arguments passed to \code{image.default}
#' @note If caching is enabled for the package (see \code{\link{gdalcubes_options}}), repeated calls of plot
#' for the same data cube will not reevaluate the cube. Instead, the temporary result file will be reused, if possible.
#' @note Some parts of the function have been copied from the stars package (c) Edzer Pebesma
#' @details 
#' The style of the plot depends on provided parameters and on the shape of the cube, i.e., whether it is a pure time series and whether it contains multiple bands or not.
#' Multi-band, multi-temporal images will be arranged with \code{layout()} such that bands are represented by columns and time is represented by rows.
#' Time series plots can be combined to a single plot by setting \code{join.timeseries = TRUE}. The layout can be controlled with \code{ncol} and \code{nrow}, which define the number of rows and columns in the plot layout. Typically, only one of 
#' \code{ncol} and \code{nrow} is provided. For multi-band, multi-temporal plots, the actual number of rows or columns can be less if the input cube has less bands or time slices.
#' 
#' The \code{downsample} argument is used to speed-up plotting if the cube has much more pixels than the graphics device. 
#' If set to a scalar integer value > 1, the value is used to skip pixels in the spatial dimensions. For example, setting \code{downsample = 4} means
#' that every fourth pixel is used in the spatial dimensions. If TRUE (the default) \code{downsample} is derived automatically based on the 
#' sizes of the cube and the graphics device. If 1 or FALSE, no additional downsampling is performed. Notice that downsampling is only used for plotting.
#' The size of the data cube (and hence the computation time to process the data cube) is not modified.
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#'               
#' plot(select_bands(raster_cube(L8.col, v), c("B02", "B03", "B04")), rgb=3:1)
#'               
#' L8.cube = select_bands(raster_cube(L8.col, v), c("B04", "B05")) 
#' L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI") 
#' plot(reduce_time(L8.ndvi, "median(NDVI)"), key.pos=1, zlim=c(0,1))
#' 
#' @export
plot.cube  <- 
  function(x,
           y,
           ...,
           nbreaks = 11,
           breaks = NULL,
           col = grey(1:(nbreaks - 1) / nbreaks),
           key.pos = NULL,
           bands = NULL,
           t = NULL,
           rgb = NULL,
           zlim = NULL,
           gamma = 1,
           periods.in.title = TRUE,
           join.timeseries = FALSE,
           axes = TRUE,
           ncol = NULL,
           nrow = NULL,
           downsample = TRUE,
           na.color = "#AAAAAA") {
    stopifnot(is.cube(x))
    size = c(nbands(x), size(x))
    
    
    def.par <-
      par(no.readonly = TRUE) 
    on.exit({
      layout(matrix(1))
      par(def.par)  
    })
    
    if (size[3] == 1 && size[4] == 1) {
      # cube is a (potentially multi-band) time-series
      if (!is.null(breaks))
        warning("data cube is a pure time-series, ignoring breaks")
      if (!is.null(key.pos))
        warning("data cube is a pure time-series, ignoring key.pos")
      if (!is.null(t))
        warning("data cube is a pure time-series, ignoring t")
      if (!is.null(rgb))
        warning("data cube is a pure time-series, ignoring rgb")
      
      dtvalues = gc_datetime_values(x)
      #if(periods.in.title) dtvalues = paste(dtvalues, cube_view(x)$time$dt)
      
      if ("ncdf_cube" %in% class(x)) {
        fn = jsonlite::parse_json(as_json(x))$file
        if (is.null(fn)) {
          stop("Invalid ncdf cube; missing file reference")
        }
      }
      else if (.pkgenv$use_cube_cache) {
        j = gc_simple_hash(as_json(x))
        if (!is.null(.pkgenv$cube_cache[[j]])
            && file.exists(.pkgenv$cube_cache[[j]])) {
          fn = .pkgenv$cube_cache[[j]]
        }
        else {
          fn = tempfile(fileext = ".nc")
          write_ncdf(x, fn)
          .pkgenv$cube_cache[[j]] = fn
        }
      }
      else {
        fn = tempfile(fileext = ".nc")
        write_ncdf(x, fn)
      }
      
      
      f <- ncdf4::nc_open(fn)
      # derive name of variables but ignore non three-dimensional variables (e.g. crs)
      vars <- names(which(sapply(f$var, function(v) {
        if (v$ndims == 3)
          return(v$name)
        return("")
      }) != ""))
      
      if (!is.null(bands)) {
        if (is.character(bands)) {
          stopifnot(all(bands %in% vars))
          vars = bands
        }
        else {
          vars = vars[bands]
        }
      }
      
      
      if (join.timeseries) {
        # if not zlim, sample data to derive ylim
        if (!is.null(ncol)) {
          warning("producing a single time series plot, ignoring ncol")
        }
        if (!is.null(nrow)) {
          warning("producing a single time series plot, ignoring nrow")
        }
        val <- NULL
        if (is.null(zlim)) {
          for (b in vars) {
            dat <- ncdf4::ncvar_get(f, b, raw_datavals = TRUE)
            val = c(val, as.vector(dat)[seq(1, size[2], length.out = min(10000 %/% size[1], size[2]))])
          }
          #zlim <- quantile(val, c(0.05, 0.95),na.rm = TRUE)
          zlim <- range(val, na.rm = TRUE)
        }
        
        #dat <- ncdf4::ncvar_get(f, b)
        if (!is.function(col)) {
          col <- rainbow(n = size[1], v = 0.5, s = 0.9)
        }
        else {
          col = col(size[1])
        }
        
        # dummy plot
        plot(
          range(1:size[2]),
          zlim,
          t = "n",
          ylim = zlim,
          xlim = range(1:size[2]),
          axes = FALSE,
          xlab = "",
          ylab = "",
          ...
        ) # add x axis labels
        
        if (length(vars) > 1) {
          for (bi in 1:length(vars)) {
            dat <- ncdf4::ncvar_get(f, vars[bi], raw_datavals = TRUE)
            lines(dat, col = col[bi], type = "b", ...)
          }
        }
        if (axes) {
          axis(2)
          axis(1, at = 1:size[2], labels = dtvalues)
        }
        box()
        legend(
          x = "topright",
          legend = vars,
          col = col,
          lty = par("lty"),
          pch = par("pch"),
          ...
        )
      }
      else {
        
        if (!is.null(ncol) && !is.null(nrow)) {
          icol = ncol
          irow = nrow
          if (size[1] > icol*irow) {
            warning("cube has more than ncol * nrow bands, some of the bands will not be plotted")
            size[1] = icol*irow
            vars=vars[1:(icol*irow)]
          }
        }
        else if (!is.null(ncol)) {
          # derive nrow
          icol = ncol
          irow = ceiling(size[1] / icol)
        }
        else if (!is.null(nrow)) {
          # derive ncol
          irow <- nrow
          icol <- ceiling(size[1] / irow)
        }
        else {
          # find a good layout
          irow <- round(sqrt(size[1]))
          icol <- ceiling(size[1] / irow)
        }

        layout(matrix(c((1:size[1]), rep(
          0, irow * icol - size[1]
        )), irow, icol, byrow = T), respect = FALSE)
        for (b in vars) {
          dat <- ncdf4::ncvar_get(f, b, raw_datavals = TRUE)
          if (!is.null(zlim)) {
            plot(
              dat,
              ylim = zlim,
              axes = FALSE,
              type = "b",
              xlab = "",
              ylab = "",
              ...
            )
          }
          else {
            plot(
              dat,
              axes = FALSE,
              type = "b",
              xlab = "",
              ylab = "",
              ...
            )
          }
          title(main = b)
          if (axes) {
            axis(2)
            axis(1, at = 1:size[2], labels = dtvalues)
          }
          box()
        }
      }
      ncdf4::nc_close(f)
      
      layout(matrix(1))
      par(def.par)  # reset to default
      
    }
    else {
      # not a time-series-only plot
      stopifnot(!(!is.null(rgb) &&
                    !is.null(bands))) # RGB and bands parameters are mutually exclusive
      if (!is.null(rgb) && !is.null(key.pos)) {
        warning("Ignoring key.pos for RGB plots")
        key.pos <- NULL
      }
      is.wholenumber <-
        function(x, tol = .Machine$double.eps ^ 0.5)
          abs(x - round(x)) < tol
      
      if (!is.null(bands)) {
        if (is.numeric(bands)) {
          stopifnot(all(is.wholenumber(bands)))
          bands = as.integer(bands)
        }
        stopifnot(is.integer(bands) || is.character(bands))
        stopifnot(anyDuplicated(bands) == 0)
        if (!is.character(bands)) {
          if (is.integer(bands)) {
            stopifnot(min(bands) >= 1 && max(bands) <= size[1])
          }
        }
      }
      else {
        bands <- 1:size[1]
      }
    
      if (!is.null(rgb)) {
        stopifnot(length(rgb) == 3)
        stopifnot(size[1] >= 3)
        if (is.numeric(rgb)) {
          stopifnot(all(is.wholenumber(rgb)))
          rgb = as.integer(rgb)
        }
        stopifnot(is.integer(rgb) || is.character(rgb))
        stopifnot(anyDuplicated(rgb) == 0)
        if (!is.character(rgb)) {
          if (is.integer(rgb)) {
            stopifnot(min(rgb) >= 1 && max(rgb) <= size[1])
          }
        }
        bands <- rgb
      }
      dtvalues = gc_datetime_values(x)
      if (periods.in.title)
        dtvalues = paste(dtvalues, dimensions(x)$t$pixel_size)
      if (!is.null(t)) {
        if (is.numeric(t)) {
          stopifnot(all(is.wholenumber(t)))
          t = as.integer(t)
        }
        stopifnot(is.integer(t))
        stopifnot(min(t) >= 1 && max(t) <= size[2])
        stopifnot(anyDuplicated(t) == 0)
        size[2] = length(t)
        dtvalues = dtvalues[t]
      }
      else {
        t <- 1:size[2]
      }
      
      stopifnot(is.numeric(gamma))
      stopifnot(gamma > 0)
      
      
      if ("ncdf_cube" %in% class(x)) {
        fn = jsonlite::parse_json(as_json(x))$file
        if (is.null(fn)) {
          stop("Invalid ncdf cube; missing file reference")
        }
      }
      else if (.pkgenv$use_cube_cache) {
        j = gc_simple_hash(as_json(x))
        if (!is.null(.pkgenv$cube_cache[[j]])
            && file.exists(.pkgenv$cube_cache[[j]])) {
          fn = .pkgenv$cube_cache[[j]]
        }
        else {
          fn = tempfile(fileext = ".nc")
          write_ncdf(x, fn)
          .pkgenv$cube_cache[[j]] = fn
        }
      }
      else {
        fn = tempfile(fileext = ".nc")
        write_ncdf(x, fn)
      }
      
  
      
      # read nc and plot individual slices as table, x = band, y = t
      def.par <-
        par(no.readonly = TRUE) # save default, for resetting...
      par(mar = c(2, 2, 2, 2))
      if (!is.null(rgb)) {
        if (!is.null(ncol) && !is.null(nrow)) {
          icol = ncol
          irow = nrow
        }
        else if (!is.null(ncol)) {
          icol = ncol
          irow = ceiling(size[2] / ncol)
        }
        else if (!is.null(nrow)) {
          irow = nrow
          icol = ceiling(size[2] / nrow)
        }
        else {
          # find a good (close to square) layout
          irow <- round(sqrt(size[2]))
          icol <- ceiling(size[2] / irow)
        }
        layout(matrix(c(1:size[2], rep(
          0, irow * icol - size[2]
        )), irow, icol, byrow = T), respect = FALSE)
      }
      else {
        if (is.null(key.pos)) {
          if (size[1] == 1 || size[2] == 1) {
            if (!is.null(ncol) && !is.null(nrow)) {
              icol = ncol
              irow = nrow
            }
            else if (!is.null(ncol)) {
              icol <- ncol
              irow <- ceiling(size[1] * size[2] / icol)
            }
            else if (!is.null(nrow)) {
              irow <-nrow
              icol <- ceiling(size[1] * size[2] / irow)
            }
            else {
              # find a good layout
              irow <- round(sqrt(size[1] * size[2]))
              icol <- ceiling(size[1] * size[2] / irow)
            }

            layout(matrix(c(
              1:(size[1] * size[2]), rep(0, irow * icol - size[1] * size[2])
            ), irow, icol, byrow = TRUE), respect = FALSE)
          }
          # general case, time dimension is vertical, bands horizontal
          else {
            icol = size[1]
            irow = size[2]
            layout(matrix(1:(size[1] * size[2]), size[2], size[1], byrow =
                            FALSE), respect = FALSE)
          }
        }
        else {
          s <- lcm(1.8)
          ww <- 1
          wh <- 1
          
          if (size[1] == 1 || size[2] == 1) {
            #icol = size[1]
            #irow = size[2]
            if (!is.null(ncol) && !is.null(nrow)) {
              icol = ncol
              irow = nrow
            }
            else if (!is.null(ncol)) {
              icol <- ncol
              irow <- ceiling(size[1] * size[2] / icol)
            }
            else if (!is.null(nrow)) {
              irow <-nrow
              icol <- ceiling(size[1] * size[2] / irow)
            }
            else {
              # find a good (close to square) layout
              irow <- round(sqrt(size[1] * size[2]))
              icol <- ceiling(size[1] * size[2] / irow)
            }
            switch (
              key.pos,
              layout(
                rbind(matrix(c(
                  1:(size[1] * size[2]), rep(0, irow * icol - size[1] * size[2])
                ), irow, icol, byrow = TRUE), size[1] * size[2] + 1),
                widths = rep.int(ww, icol),
                heights = c(rep.int(wh, irow), s),
                respect = FALSE
              ),
              layout(
                cbind(size[1] * size[2] + 1, matrix(c(
                  1:(size[1] * size[2]), rep(0, irow * icol - size[1] * size[2])
                ), irow, icol, byrow = TRUE)),
                widths = c(s, rep.int(ww, icol)),
                heights = rep.int(wh, irow),
                respect = FALSE
              ),
              layout(
                rbind(size[1] * size[2] + 1, matrix(c(
                  1:(size[1] * size[2]), rep(0, irow * icol - size[1] * size[2])
                ), irow, icol, byrow = TRUE)),
                widths = rep.int(ww, icol),
                heights = c(s, rep.int(wh, irow)),
                respect = FALSE
              ),
              layout(
                cbind(matrix(c(
                  1:(size[1] * size[2]), rep(0, irow * icol - size[1] * size[2])
                ), irow, icol, byrow = TRUE), size[1] * size[2] + 1),
                widths = c(rep.int(ww, icol), s),
                heights = c(rep.int(wh, irow)),
                respect = FALSE
              )
            )
          }
          else {
            # general case, time dimension is vertical, bands horizontal
            icol = size[1]
            irow = size[2]
            switch (
              key.pos,
              layout(
                rbind(matrix(
                  1:(size[1] * size[2]), size[2], size[1], byrow = FALSE
                ), size[1] * size[2] + 1),
                widths = rep.int(ww, size[1]),
                heights = c(rep.int(wh, size[2]), s),
                respect = FALSE
              ),
              layout(
                cbind(size[1] * size[2] + 1, matrix(
                  1:(size[1] * size[2]), size[2], size[1], byrow = FALSE
                )),
                widths = c(s, rep.int(ww, size[1])),
                heights = rep.int(wh, size[2]),
                respect = FALSE
              ),
              layout(
                rbind(size[1] * size[2] + 1, matrix(
                  1:(size[1] * size[2]), size[2], size[1], byrow = FALSE
                )),
                widths = rep.int(ww, size[1]),
                heights = c(s, rep.int(wh, size[2])),
                respect = FALSE
              ),
              layout(
                cbind(matrix(
                  1:(size[1] * size[2]), size[2], size[1], byrow = FALSE
                ), size[1] * size[2] + 1),
                widths = c(rep.int(ww, size[1]), s),
                heights = c(rep.int(wh, size[2])),
                respect = FALSE
              )
            )
            
          }
          
        }
      }
      
      dev_size = dev.size("px")
      if (is.null(downsample) || isTRUE(downsample)) {
        ds_x = max(floor(size[4] / (dev_size[1] / icol)),1)
        ds_y = max(floor(size[3] / (dev_size[2] / irow)),1)
        downsample = min(ds_x, ds_y)
      }
      else if (isFALSE(downsample)) {
        downsample = 1
      }
      else {
        if (!is.wholenumber(downsample)) {
          ds_x = max(floor(size[4] / (dev_size[1] / icol)),1)
          ds_y = max(floor(size[3] / (dev_size[2] / irow)),1)
          downsample = min(ds_x, ds_y)
          warning(paste0("Provided downsampling factor is not an integer number; using ", downsample, " instead"))
        }
        if (downsample < 1) {
          ds_x = max(floor(size[4] / (dev_size[1] / icol)),1)
          ds_y = max(floor(size[3] / (dev_size[2] / irow)),1)
          downsample = min(ds_x, ds_y)
          warning(paste0("Provided downsampling factor is < 1; using ", downsample, " instead"))
        }
      }
  
      dims <- dimensions(x)
      
      f <- ncdf4::nc_open(fn)
      # dtunit = strsplit(f$dim$time$units, " since ")[[1]]
      # dtunit_str = switch(dtunit[1],
      #        years = paste("years since", format(strptime(dtunit[2], format="%Y-%m-%dT%H:%M:%S")),"%Y"),
      #        months = paste("months since", format(strptime(dtunit[2], format="%Y-%m-%dT%H:%M:%S")),"%Y-%m"),
      #        days = paste("days since", format(strptime(dtunit[2], format="%Y-%m-%dT%H:%M:%S")),"%Y-%m-%d"),
      #        hours = paste("hours since", format(strptime(dtunit[2], format="%Y-%m-%dT%H:%M:%S")),"%Y-%m-%dT%H:%M:%S"),
      #        minutes = paste("minutes since", format(strptime(dtunit[2], format="%Y-%m-%dT%H:%M:%S")),"%Y-%m-%dT%H:%M:%S"),
      #        seconds = paste("seconds since", format(strptime(dtunit[2], format="%Y-%m-%dT%H:%M:%S")),"%Y-%m-%dT%H:%M:%S"))
      #
      
      # derive name of variables but ignore non three-dimensional variables (e.g. crs)
      vars <- names(which(sapply(f$var, function(v) {
        if (v$ndims == 3)
          return(v$name)
        return("")
      }) != ""))
      
      if (!is.null(bands)) {
        if (is.character(bands)) {
          stopifnot(all(bands %in% vars))
          vars = bands
        }
        else {
          vars = vars[bands]
        }
      }
      
      
      
      dimsx = seq(dims$x$low, dims$x$high, length.out = size[4])
      dimsy = seq(dims$y$low, dims$y$high, length.out = size[3])
      asp = ((dims$y$high - dims$y$low) / dims$y$count) / ((dims$x$high - dims$x$low) / dims$x$count)
      ylim = c(dims$y$low, dims$y$high)
      xlim = c(dims$x$low, dims$x$high)
      
      #dimst = seq(from = dims$low[1], by = (dims$high[1] - dims$low[1] + 1) %/% size[2], length.out = size[2])
      
      
      # if breaks will be computed from the data,
      # read data from all bands to get a good sample of pixels
      # TODO avoid reading the data twice
      if (is.null(rgb)) {
        if (is.null(breaks)) {
          if (is.null(zlim)) {
            val <- NULL
            for (b in vars) {
              if (!is.null(bands)) {
                if (is.character(bands)) {
                  if (b %in% bands)
                    next
                }
              }
              dat <- ncdf4::ncvar_get(f, b, raw_datavals = TRUE)
              if (length(dim(dat)) == 2) {
                val = c(val, as.vector(dat)[seq(1,
                                                prod(size[2:4]),
                                                length.out = min(10000 %/% size[1], prod(size[2:4])))])
              }
              else {
                val = c(val, as.vector(dat[, , t])[seq(1,
                                                       prod(size[2:4]),
                                                       length.out = min(10000 %/% size[1], prod(size[2:4])))])
              }
            }
            zlim[1] = min(dat, na.rm = TRUE)
            zlim[2] = max(dat, na.rm = TRUE)
            breaks = seq(zlim[1], zlim[2], length.out = nbreaks)
          }
          else {
            breaks = seq(zlim[1], zlim[2], length.out = nbreaks)
          }
        }
        stopifnot(length(breaks) ==  nbreaks) # TODO: clean up graphics state
        
        if (is.function(col)) {
          col = col(n = nbreaks - 1, ...)
        }
        col = c(col, na.color)
        nbreaks = nbreaks + 1
        breaks=c(breaks, breaks[length(breaks)] + diff(range(breaks))/length(breaks))
        stopifnot(length(col) == nbreaks - 1)
      }
      
      
      if (!is.null(rgb)) {
        dat_R <- ncdf4::ncvar_get(f, vars[1], raw_datavals = TRUE)
        dat_G <- ncdf4::ncvar_get(f, vars[2], raw_datavals = TRUE)
        dat_B <- ncdf4::ncvar_get(f, vars[3], raw_datavals = TRUE)
        
        
        rng_R <- range(dat_R, na.rm = T, finite = T)
        rng_G <- range(dat_G, na.rm = T, finite = T)
        rng_B <- range(dat_B, na.rm = T, finite = T)
        
        if (is.null(zlim)) {
          zlim <- quantile(c(dat_R, dat_G, dat_B), c(0.1, 0.9), na.rm = T)
        }
        
        #rng <- range(c(rng_R, rng_B, rng_G), na.rm = T, finite=T)
        scale  = (zlim[2] - zlim[1])
        offset = zlim[1]
        
        dat_R <- (dat_R - offset) / scale
        dat_G <- (dat_G - offset) / scale
        dat_B <- (dat_B - offset) / scale
        
        if (gamma != 1) {
          dat_R = dat_R^gamma
          dat_G = dat_G^gamma
          dat_B = dat_B^gamma
        }
        
        dat_R[which(is.na(dat_R) , arr.ind = T)] <- col2rgb(na.color)[1] / 255
        dat_G[which(is.na(dat_G) , arr.ind = T)] <- col2rgb(na.color)[2] / 255
        dat_B[which(is.na(dat_B) , arr.ind = T)] <- col2rgb(na.color)[3] / 255
        
        
        dat_R[which(dat_R < 0 , arr.ind = T)] <- 0
        dat_G[which(dat_G < 0 , arr.ind = T)] <- 0
        dat_B[which(dat_B < 0 , arr.ind = T)] <- 0
        
        dat_R[which(dat_R > 1 , arr.ind = T)] <- 1
        dat_G[which(dat_G > 1 , arr.ind = T)] <- 1
        dat_B[which(dat_B > 1 , arr.ind = T)] <- 1
        
        #asp = diff(ylim) / diff(xlim)
        
        xlim_new = xlim
        ylim_new = ylim
        
        for (ti in 1:size[2]) {
          plot(
            1,
            1,
            t = "n",
            ylim = c(ylim_new[1], ylim_new[2]),
            xlim = c(xlim_new[1], xlim_new[2]),
            asp = asp,
            axes = axes,
            xlab = "",
            ylab = "",
            xaxs = "i",
            yaxs = "i"
          )
          
          if (length(dim(dat_R)) == 2) {
            #ar <- array(c(dat_R[size[4]:1,], dat_G[size[4]:1,], dat_B[size[4]:1,]),dim=c(dim(dat_R)[1],dim(dat_R)[2], 3))
            ar <- array(c(dat_R[, ], dat_G[, ], dat_B[, ]), dim = c(dim(dat_R)[1], dim(dat_R)[2], 3))
            
            if (downsample > 1) {
              ds_1_idx = seq(1,dim(ar)[1], by = downsample) 
              ds_2_idx = seq(1,dim(ar)[2], by = downsample) 
              ar = ar[ds_1_idx,ds_2_idx,]
            }
            
            rasterImage(
              aperm(ar, c(2, 1, 3)),
              xleft = xlim[1],
              xright = xlim[2],
              ybottom = ylim[1],
              ytop = ylim[2]
            ) # TODO: add interpolate argument
          }
          else {
            #ar <- array(c(dat_R[size[4]:1,,t[ti]], dat_G[size[4]:1,,t[ti]], dat_B[size[4]:1,,t[ti]]),dim=c(dim(dat_R)[1],dim(dat_R)[2], 3))
            ar <-
              array(c(dat_R[, , t[ti]], dat_G[, , t[ti]], dat_B[, , t[ti]]), dim = c(dim(dat_R)[1], dim(dat_R)[2], 3))
            if (downsample > 1) {
              ds_1_idx = seq(1,dim(ar)[1], by = downsample) 
              ds_2_idx = seq(1,dim(ar)[2], by = downsample) 
              ar = ar[ds_1_idx,ds_2_idx,]
            }
            rasterImage(
              aperm(ar, c(2, 1, 3)) ,
              xleft = xlim[1],
              xright = xlim[2],
              ybottom = ylim[1],
              ytop = ylim[2]
            ) # TODO: add interpolate argument
          }
          title(dtvalues[ti])
          box()
        }
        
      }
      else {
        # non-RGB plot
        for (b in vars) {
          #ncdf4::ncvar_get(f, b, start=c(t,1,1), count = c(1, f$dim[[2]]$len,f$dim[[3]]$len))
          dat <- ncdf4::ncvar_get(f, b, raw_datavals = TRUE)
          
          dat[which(dat < breaks[1] | dat > breaks[length(breaks)-1], arr.ind = TRUE)] <- breaks[1] - 1 # do not plot values outside breaks / zlim
          dat[which(is.na(dat),arr.ind = TRUE)] <- breaks[length(breaks)] # plot NA with special color
          
          xaxt = "s"
          yaxt = "s"
          if (!axes) {
            xaxt = "n"
            yaxt = "n"
          }
          
          for (ti in 1:size[2]) {
            #image(aperm(dat[,,t], 2:1))
            if (length(dim(dat)) == 2) {
              ar = dat[, size[3]:1]
              dimsx_new = dimsx
              dimsy_new = dimsy
              if (downsample > 1) {
                ds_1_idx = seq(1,dim(ar)[1], by = downsample) 
                ds_2_idx = seq(1,dim(ar)[2], by = downsample) 
                dimsx_new = dimsx[ds_1_idx]
                dimsy_new = dimsy[ds_2_idx]
                ar = ar[ds_1_idx,ds_2_idx]
              }
              image.default(
                dimsx_new,
                dimsy_new,
                ar,
                col = col,
                asp = asp,
                xaxt = xaxt,
                yaxt = yaxt,
                breaks = breaks,
                xlim = xlim,
                ylim = ylim,
                useRaster = dev.capabilities("rasterImage")$rasterImage == "yes",
                ...
              )
              title(paste(b, " | ",  dtvalues[ti], sep = ""))
            }
            else {
              ar =  dat[, size[3]:1, t[ti]]
              dimsx_new = dimsx
              dimsy_new = dimsy
              if (downsample > 1) {
                ds_1_idx = seq(1,dim(ar)[1], by = downsample) 
                ds_2_idx = seq(1,dim(ar)[2], by = downsample) 
                dimsx_new = dimsx[ds_1_idx]
                dimsy_new = dimsy[ds_2_idx]
                ar = ar[ds_1_idx,ds_2_idx]
              }
              image.default(
                dimsx_new,
                dimsy_new,
                ar,
                col = col,
                asp = asp,
                xaxt = xaxt,
                yaxt = yaxt,
                breaks = breaks,
                xlim = xlim,
                ylim = ylim,
                useRaster = dev.capabilities("rasterImage")$rasterImage == "yes",
                ...
              )
              title(paste(b, " | ",  dtvalues[ti], sep = ""))
            }
          }
        }
      }
      
      ncdf4::nc_close(f)
      
      
      if (!is.null(key.pos)) {
        #plot.new()
        
        if (key.pos %in% c(1, 3)) {
          if (key.pos == 1)
            par(mar = c(2.1, 2, 0.5, 2))
          else
            par(mar = c(0.5, 2, 2.1, 2))
          plot(
            range(breaks)[1],
            0,
            t = "n",
            ylim = c(0, 1),
            xlim = range(breaks[1:(length(breaks)-1)]),
            axes = FALSE,
            xlab = "",
            ylab = "",
            xaxs = "i",
            yaxs = "i"
          )
        }
        
        if (key.pos %in% c(2, 4)) {
          if (key.pos == 2)
            par(mar = c(2, 2.1, 2, 0.5))
          else
            par(mar = c(2, 0.5, 2, 2.1))
          
          plot(
            0,
            range(breaks)[1],
            t = "n",
            ylim = range(breaks[1:(length(breaks)-1)]),
            xlim = c(0, 1),
            axes = FALSE,
            xlab = "",
            ylab = "",
            xaxs = "i",
            yaxs = "i"
          )
        }
        
        
        
        
        for (i in 1:length(col)) {
          switch(
            key.pos,
            rect(
              xleft = breaks[1:(length(breaks) - 2)],
              ybottom = 0 ,
              xright = breaks[2:(length(breaks) - 1)],
              ytop = 1,
              col = col[1:(length(col)-1)],
              border = NA
            ),
            rect(
              ybottom = breaks[1:(length(breaks) - 2)],
              xleft = 0 ,
              ytop = breaks[2:(length(breaks) - 1)],
              xright = 1,
              col = col[1:(length(col)-1)],
              border = NA
            ),
            rect(
              xleft = breaks[1:(length(breaks) - 2)],
              ybottom = 0 ,
              xright = breaks[2:(length(breaks) - 1)],
              ytop = 1,
              col = col[1:(length(col)-1)],
              border = NA
            ),
            rect(
              ybottom = breaks[1:(length(breaks) - 2)],
              xleft = 0 ,
              ytop = breaks[2:(length(breaks) - 1)],
              xright = 1,
              col = col[1:(length(col)-1)],
              border = NA
            )
          )
        }
        box()
        axis(key.pos, at = pretty(range(breaks)))
      }
    }
  }
