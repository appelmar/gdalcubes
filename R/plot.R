#' Plot a gdalcubes data cube
#'
#' @param x a data cube proxy object (class gcbs_cube)
#' @param y _not used_
#' @param nbreaks number of breaks, should be one more than the number of colors given
#' @param breaks actual breaks used to assign colors to values, if missing, the function subsamples values and uses equally sized intervals between min and max or zlim[0] and zlim[1] if provided
#' @param col color definition, can be a character vector with nbreaks - 1 elements or a function such as heat.colors
#' @param key.pos position for the legend, 1 (bottom), 2 (left), 3 (top), or 4 (right). If NULL, do not plot a legend.
#' @param bands integer vector with band numbers to plot (currently this must be band numbers, not band names)
#' @param t integer vector time indexes to plot (currently this must be time index, not date / time)
#' @param rgb bands used to assign RGB color channels, vector of length 3  (currently this must be band numbers, not band names)
#' @param zlim vector of length 2, defining the maximum and minimum values to either derive breaks, or define black and white values in RGB plots
#' @param ... further arguments passed to \code{image.default}
#' @note There is currently no way to plot the result of \çode{gcbs_eval} without reevaluating the cube
#' @export
plot.gcbs_cube  <- function(x, y, ..., nbreaks=11, breaks=NULL,col=grey(1:(nbreaks-1)/nbreaks), key.pos=NULL, bands=NULL, t=NULL, rgb=NULL, zlim=NULL) {
  stopifnot(is.gcbs_cube(x))
  size = c(gcbs_nbands(x), gcbs_size(x))
  stopifnot(!(!is.null(rgb) && !is.null(bands))) # RGB and bands parameters ase mutually exclusive
  if (!is.null(rgb) && !is.null(key.pos)) {
    warning("Ignoring key.pos for RGB plots")
    key.pos <- NULL
  }
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.null(bands)) {
    if (is.numeric(bands)) {
      stopifnot(all(is.wholenumber(bands)))
      bands = as.integer(bands)
    }
    stopifnot(is.integer(bands) || is.character(bands))
    stopifnot(anyDuplicated(bands) == 0)
    if (!is.character(bands)) {
      if(is.integer(bands)) {
        stopifnot(min(bands) >= 1 && max(bands) <= size[1])
      }
    }
    x = gcbs_select_bands(x, as.character(gcbs_bands(x)$name[bands])) # optimization to only store needed bands
    bands = 1:length(bands)
    size[1] = length(bands)
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
      if(is.integer(rgb)) {
        stopifnot(min(rgb) >= 1 && max(rgb) <= size[1])
      }
    }
    bands <- rgb
    x = gcbs_select_bands(x, as.character(gcbs_bands(x)$name[bands])) # optimization to only store needed bands
    rgb = 1:3
    bands = 1:3
    size[1] <- 3
  }
  dtvalues = libgdalcubes_datetime_values(x)
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
  if (is.null(rgb)) {
    if (size[1] > 5)
    {
      warning("too many bands, plotting only bands 1 to 5")
      bands = bands[1:5]
      size[1] = 5
    }
    if (size[2] > 5) {
      warning("too many time instances, plotting only times 1 to 5")
      t = t[1:5]
      size[2] = 5
    }
  }
  else {
    if (size[2] > 25) {
      warning("too many time instances, plotting only times 1 to 25")
      t = t[1:25]
      size[2] = 25
    }
  }


  fn = tempfile(fileext=".nc")
  gcbs_eval(x, fn)

  # read nc and plot individual slices as table, x = band, y = t
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  par(mar=c(2,2,2,2))
  if(!is.null(rgb)) {
    ncol = min(5, size[2])
    nrow = ceiling(size[2] / ncol)
    layout(matrix(c(1:size[2], rep(0,nrow*ncol-size[2])), nrow, ncol, byrow=T), respect = TRUE)
  }
  else {
    if (is.null(key.pos))
      layout(matrix(1:(size[1]*size[2]), size[2], size[1], byrow=T), respect = TRUE)
    else {
      s <- lcm(1.8)
      ww <- 1
      wh <- 1
      switch (key.pos,
              layout(rbind(matrix(1:(size[1]*size[2]), size[2], size[1], byrow=T),size[1]*size[2]+1), widths = rep.int(ww,size[1]), heights = c(rep.int(wh,size[2]), s), respect = TRUE),
              layout(cbind(size[1]*size[2]+1, matrix(1:(size[1]*size[2]), size[2], size[1], byrow=T)), widths = c(s, rep.int(ww,size[1])), heights = rep.int(wh,size[2]), respect = TRUE),
              layout(rbind(size[1]*size[2]+1, matrix(1:(size[1]*size[2]), size[2], size[1], byrow=T)), widths = rep.int(ww,size[1]), heights = c(s, rep.int(wh,size[2])), respect = TRUE),
              layout(cbind(matrix(1:(size[1]*size[2]), size[2], size[1], byrow=T),size[1]*size[2]+1), widths = c(rep.int(ww,size[1]), s), heights = c(rep.int(wh,size[2])), respect = TRUE))

    }
  }

  dims <- gcbs_dimensions(x)

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
      stopifnot( all(bands %in% vars))
      vars = bands
    }
    else {
      vars = vars[bands]
    }
  }



  dimsx = seq(dims$low[3], dims$high[3], length.out = size[4])
  dimsy = seq(dims$low[2], dims$high[2], length.out = size[3])
  asp = ((dims$high[2] - dims$low[2])/dims$size[2])/((dims$high[3] - dims$low[3])/dims$size[3])
  ylim = c(dims$low[2], dims$high[2])
  xlim =c(dims$low[3], dims$high[3])

  #dimst = seq(from = dims$low[1], by = (dims$high[1] - dims$low[1] + 1) %/% size[2], length.out = size[2])


  # if breaks will be computed from the data,
  # read data from all bands to get a good sample of pixels
  # TODO avoid reading the data twice
  if (is.null(rgb)) {
    if(is.null(breaks)) {
      if (is.null(zlim)) {
        val <- NULL
        for (b in vars) {
          if (!is.null(bands)) {
            if (is.character(bands)) {
              if (b %in% bands)
                next
            }
          }
          dat <- ncdf4::ncvar_get(f, b)
          if (length(dim(dat)) == 2) {
            val = c(val, as.vector(dat)[seq(1,prod(size[2:4]), length.out = min(10000 %/% size[1], prod(size[2:4])))])
          }
          else {
            val = c(val, as.vector(dat[,,t])[seq(1,prod(size[2:4]), length.out = min(10000 %/% size[1], prod(size[2:4])))])
          }
        }
        breaks = seq(min(dat,na.rm = TRUE), max(dat,na.rm = TRUE), length.out = nbreaks)
      }
      else {
        breaks = seq(zlim[1], zlim[2], length.out = nbreaks)
      }
      #breaks = quantile(dat, seq(0,1,length.out=nbreaks+2)[2:(nbreaks+1)], na.rm=TRUE)
    }
    stopifnot(length(breaks) ==  nbreaks) # TODO: clean up graphics state

    if (is.function(col)) {
      col = col(n=nbreaks-1, ...)
    }
    stopifnot(length(col) == nbreaks - 1)
  }


  if (!is.null(rgb)) {
    dat_R <- ncdf4::ncvar_get(f, vars[1])
    dat_G <- ncdf4::ncvar_get(f, vars[2])
    dat_B <- ncdf4::ncvar_get(f, vars[3])


    rng_R <- range(dat_R, na.rm = T, finite=T)
    rng_G <- range(dat_G, na.rm = T, finite=T)
    rng_B <- range(dat_B, na.rm = T, finite=T)

    if (is.null(zlim)) {
      zlim <- quantile(c(dat_R, dat_G, dat_B), c(0.1, 0.9), na.rm = T)
    }

    #rng <- range(c(rng_R, rng_B, rng_G), na.rm = T, finite=T)
    scale  = (zlim[2]-zlim[1])
    offset = zlim[1]

    dat_R <- (dat_R - offset)/scale
    dat_G <- (dat_G - offset)/scale
    dat_B <- (dat_B - offset)/scale

    dat_R[which(is.nan(dat_R) ,arr.ind = T)] <- 1
    dat_G[which(is.nan(dat_G) ,arr.ind = T)] <- 1
    dat_B[which(is.nan(dat_B) ,arr.ind = T)] <- 1


    dat_R[which(dat_R < 0 ,arr.ind = T)] <- 0
    dat_G[which(dat_G < 0 ,arr.ind = T)] <- 0
    dat_B[which(dat_B < 0 ,arr.ind = T)] <- 0

    dat_R[which(dat_R > 1 ,arr.ind = T)] <- 1
    dat_G[which(dat_G > 1 ,arr.ind = T)] <- 1
    dat_B[which(dat_B > 1 ,arr.ind = T)] <- 1

    for (ti in 1:size[2]) {
      plot(1, 1, t = "n", ylim = c(0,1), xlim = c(0,1), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")

      if (length(dim(dat_R)) == 2) {
        #ar <- array(c(dat_R[size[4]:1,], dat_G[size[4]:1,], dat_B[size[4]:1,]),dim=c(dim(dat_R)[1],dim(dat_R)[2], 3))
        ar <- array(c(dat_R[,], dat_G[,], dat_B[,]),dim=c(dim(dat_R)[1],dim(dat_R)[2], 3))


        rasterImage(aperm(ar, c(2,1,3)), xleft = 0, xright = 1, ybottom = 0, ytop = 1) # TODO: add interpolate argument
      }
      else {
        #ar <- array(c(dat_R[size[4]:1,,t[ti]], dat_G[size[4]:1,,t[ti]], dat_B[size[4]:1,,t[ti]]),dim=c(dim(dat_R)[1],dim(dat_R)[2], 3))
        ar <- array(c(dat_R[,,t[ti]], dat_G[,,t[ti]], dat_B[,,t[ti]]),dim=c(dim(dat_R)[1],dim(dat_R)[2], 3))
        rasterImage(aperm(ar, c(2,1,3)) ,xleft = 0, xright = 1, ybottom = 0, ytop = 1) # TODO: add interpolate argument
      }
      title(dtvalues[ti])
      box()
      #TODO axis??
    }

  }
  else {
    for (b in vars) {
      #ncdf4::ncvar_get(f, b, start=c(t,1,1), count = c(1, f$dim[[2]]$len,f$dim[[3]]$len))
      dat <- ncdf4::ncvar_get(f, b)

      for (ti in 1:size[2]) {
        #image(aperm(dat[,,t], 2:1))
        if (length(dim(dat)) == 2) {
          # add asp?
          image.default(dimsx, dimsy, dat[,size[3]:1], col=col, asp=asp,  breaks=breaks, xlim=xlim, ylim=ylim, ...)
          title(paste(b, " | ",  dtvalues[ti], sep="")) # TODO: replace t with string
        }
        else {
          # add asp?
          image.default(dimsx, dimsy, dat[,size[3]:1,t[ti]], col=col, asp=asp, breaks=breaks, xlim=xlim, ylim=ylim, ...)
          #title(paste(b, " | ",  t[ti], "[", dtunit_str, "]", sep="")) # TODO: replace t with string
          title(paste(b, " | ",  dtvalues[ti], sep="")) # TODO: replace t with string
          }
      }
    }
  }

  ncdf4::nc_close(f)


  if (!is.null(key.pos)) {
    #plot.new()

    if (key.pos %in% c(1,3)) {
      if (key.pos == 1)
        par(mar=c(2.1,2,0.5,2))
      else
        par(mar=c(0.5,2,2.1,2))
      plot(range(breaks)[1], 0, t = "n", ylim = c(0,1), xlim = range(breaks), axes = FALSE,  xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    }

    if (key.pos %in% c(2,4)) {
      if (key.pos == 2)
        par(mar=c(2,2.1,2,0.5))
      else
        par(mar=c(2,0.5,2,2.1))

      plot(0, range(breaks)[1], t = "n", ylim = range(breaks), xlim = c(0,1), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    }




    for (i in 1:length(col)) {
      switch(key.pos,
             rect(xleft = breaks[1:(length(breaks)-1)], ybottom = 0 ,xright = breaks[2:(length(breaks))], ytop = 1, col=col, border=NA),
             rect(ybottom = breaks[1:(length(breaks)-1)], xleft = 0 ,ytop = breaks[2:(length(breaks))], xright = 1, col=col, border=NA),
             rect(xleft = breaks[1:(length(breaks)-1)], ybottom = 0 ,xright = breaks[2:(length(breaks))], ytop = 1, col=col, border=NA),
             rect(ybottom = breaks[1:(length(breaks)-1)], xleft = 0 ,ytop = breaks[2:(length(breaks))], xright = 1, col=col, border=NA))
    }
    box()
    axis(key.pos, at = pretty(range(breaks)))
  }


  layout(matrix(1))
  par(def.par)  # reset to default

}




