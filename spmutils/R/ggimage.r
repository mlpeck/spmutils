ggimage <- function(zmat, x=NULL, y=NULL, col=viridis::viridis(256),
                    xlab=NULL, ylab=NULL, legend=NULL, title=NULL,
                    xrev=FALSE, yrev=FALSE, asp=1,
                    addcontour=FALSE, binwidth=NULL
                   ) {
    require(ggplot2)
    if (is.null(x)) x <- 1:nrow(zmat)
    if (is.null(y)) y <- 1:ncol(zmat)
    if (!is.null(dim(x))) {
      x <- colMeans(x)
      y <- rowMeans(y)
    }
    xy <- expand.grid(x,y)
    df <- data.frame(cbind(xy, as.vector(zmat)))
    names(df) <- c("x","y","z")
    g1 <- ggplot(df, aes(x=x, y=y, z=z)) + 
        geom_raster(aes(fill=z)) + 
        scale_fill_gradientn(colors=col, na.value="#FFFFFF00") +
        coord_fixed(ratio = asp)
    if (addContour) {
        if (!is.null(binwidth)) {
            g1 <- g1 + geom_contour(binwidth=binwidth, na.rm=TRUE)
        } else {
            g1 <- g1 + geom_contour(na.rm=TRUE)
        }
    }
    if (xrev) g1 <- g1 + scale_x_reverse()
    if (yrev) g1 <- g1 + scale_y_reverse()
    if (!is.null(xlab)) {
        g1 <- g1 + xlab(xlab)
    }
    if (!is.null(ylab)) {
        g1 <- g1 + ylab(ylab)
    }
    if (!is.null(legend)) {
        g1 <- g1 + labs(fill=legend)
    }
    if (!is.null(title)) {
        g1 <- g1 + ggtitle(title)
    }
    invisible(list(df=df, graph=g1))
}

## simple 2D spatial interpolation

fillpoly <- function(ra, dec, zvals, dxy=0.5, min_ny=100, usefields=TRUE) {
    ny <- max(round((max(dec)-min(dec))*3600/dxy), min_ny)
    nx <- ny
    if (usefields) {
      xy <- cbind(ra, dec)
      zxy <- fields::Tps(x=xy, Y=zvals, lon.lat=TRUE, miles=FALSE, lambda=0)
      zsurf <- fields::predictSurface(zxy, nx=nx, ny=ny)
    } else {
      allok <- complete.cases(ra, dec, zvals)
      zsurf <- akima::interp(x=ra[allok], y=dec[allok], z=zvals[allok], nx=nx, ny=ny)
    }
    list(ra=zsurf$x, dec=zsurf$y, z=zsurf$z)
}

ggbinplot <- function (gdat, z, zlab = NULL, addfiberpos = TRUE, addcentroid = FALSE, 
                addborders = TRUE, addfp = FALSE, addcontour = FALSE,
                show.legend = TRUE, na.value = "grey95", 
                palette = "Set1", colors = viridis::viridis(256), 
                contourcolor="black", cbinwidth=10, fpcolor = "red") {
  require(ggplot2)
  df <- data.frame(x = gdat$xpos, y = gdat$ypos, z = z)
  if (exists("x.orig", gdat)) {
    df2 <- data.frame(x = gdat$x.orig, y = gdat$y.orig)
  }
  else {
    df2 <- df[ , -3]
  }
  df3 <- df2[grDevices::chull(df2), ]
  df3[df3 < 0] <- df3[df3 < 0] - 0.7
  df3[df3 >= 0] <- df3[df3 >= 0] + 0.7
  prange <- c(range(df2$x), range(df2$y)) + c(-1, 1, -1, 1)
  g1 <- ggplot(df) + ggforce::geom_voronoi_tile(aes(x = x, 
    y = y, fill = z, group=-1L), na.rm = TRUE, bound = df3, show.legend = show.legend)
  if (addborders) {
    g1 <- g1 + ggforce::geom_voronoi_segment(aes(x = x, 
      y = y), bound = prange)
  }
  if (addcentroid) {
    g1 <- g1 + geom_point(aes(x = x, y = y), shape = 21, 
      size = 2)
  }
  if (addfiberpos) {
    g1 <- g1 + geom_point(aes(x = x, y = y), data = df2)
  }
  if (addcontour) {
    allok <- complete.cases(df)
    zsurf <- akima::interp(x=df$x[allok], y=df$y[allok], z=df$z[allok])
    xy <- expand.grid(zsurf$x, zsurf$y)
    df4 <- data.frame(x=xy[,1], y=xy[,2], z=as.vector(zsurf$z))
    g1 <- g1 + geom_contour(aes(x=x, y=y, z=z), data=df4, 
                            color=contourcolor, binwidth=cbinwidth, na.rm=TRUE)
  }
  if (addfp) {
    g1 <- g1 + ggalt::geom_encircle(aes(x = x, y = y), data = df2, 
      expand = 0.02, color = fpcolor, size = 2)
  }
  if (!is.null(zlab) & show.legend) {
    g1 <- g1 + labs(fill = zlab)
  }
  if (is.factor(z)) {
    g1 <- g1 + scale_fill_brewer(type = "qual", palette = palette, 
      na.value = na.value, drop = FALSE)
  }
  else {
    g1 <- g1 + scale_fill_gradientn(colors = colors, na.value = na.value)
  }
  plot(g1)
  invisible(list(df=df, graph=g1))
}


plotci <- function(df, names,
                   col=NULL, shape=NULL, palette="Set1", 
                   legend=NULL, legend2=NULL,
                   xlab=NULL, ylab=NULL, title=NULL) {
  require(ggplot2)
  if (length(names) == 4) {
    if (!is.na(names[2])) {
      xmin <- df[, names[1]] - df[, names[2]]
      xmax <- df[, names[1]] + df[, names[2]]
    }
    ymin <- df[, names[3]] - df[, names[4]]
    ymax <- df[, names[3]] + df[, names[4]]
    yname <- as.name(names[3])
  } else if (length(names) == 6) {
    if (!is.na(names[2])) {
      xmin <- df[, names[2]]
      xmax <- df[, names[3]]
    }
    ymin <- df[, names[5]]
    ymax <- df[, names[6]]
    yname <- as.name(names[4])
  } else stop("Need 4 or 6 names to plot")

  g1 <- ggplot(df, aes_string(x=as.name(names[1]), y=yname), na.rm=TRUE)
  if (is.null(shape)) {
        g1 <- g1 + geom_point(na.rm=TRUE)
  } else {
      g1 <- g1 + geom_point(aes_string(shape=as.name(shape)), na.rm=TRUE)
  }
  g1 <- g1 + geom_errorbar(aes(ymin = ymin, ymax = ymax), na.rm=TRUE)
  if (!is.null(col)) {
    g1 <- g1 + geom_point(aes_string(color = as.name(col)), na.rm=TRUE)
    g1 <- g1 + geom_errorbar(aes_string(ymin=ymin, ymax=ymax, color = as.name(col)), na.rm=TRUE)
  }
  if (!is.na(names[2])) {
    g1 <- g1 + geom_errorbarh(aes(xmin = xmin, xmax = xmax), na.rm=TRUE)
    if (!is.null(col)) {
      g1 <- g1 + geom_errorbarh(aes_string(xmin = xmin, xmax = xmax, color = as.name(col)), na.rm=TRUE)
    }
  }
  if (is.factor(col)) {
    g1 <- g1 + scale_colour_brewer(drop=FALSE, type="qual", na.value=NA,
                                       palette=palette)
  }
  if (!is.null(legend)) {
    g1 <- g1 + labs(color=legend)
  }
  if (!is.null(legend2)) {
    g1 <- g1 + labs(shape=legend2)
  }
  if (!is.null(xlab)) {
    g1 <- g1 + xlab(xlab)
  }
  if (!is.null(ylab)) {
    g1 <- g1 + ylab(ylab)
  }
  if (!is.null(title)) {
    g1 <- g1 + ggtitle(title)
  }
  g1
}

ggmapfac <- function(fac, x, y, palette="Set1",
                     xlab=NULL, ylab=NULL, legend=NULL, title=NULL,
                     xrev=FALSE, yrev=FALSE, asp=1) {
    require(ggplot2)
    options("warn" = -1)
    if (is.null(x)) x <- 1:nrow(zmat)
    if (is.null(y)) y <- 1:ncol(zmat)
    if (!is.null(dim(x))) {
      x <- colMeans(x)
      y <- rowMeans(y)
    }
    xy <- expand.grid(x,y)
    xyz <- data.frame(cbind(xy, fac))
    names(xyz) <- c("x","y","fac")
    g1 <- ggplot(xyz, aes(x=x,y=y))+geom_tile(aes(fill=fac)) +
        scale_fill_brewer(drop=FALSE, type="qual",na.value="#FFFFFF00",palette=palette) +
        coord_fixed(ratio = asp)
    if (xrev) g1 <- g1 + scale_x_reverse()
    if (yrev) g1 <- g1 + scale_y_reverse()
    if (!is.null(xlab)) {
        g1 <- g1 + xlab(xlab)
    }
    if (!is.null(ylab)) {
        g1 <- g1 + ylab(ylab)
    }
    if (!is.null(legend)) {
        g1 <- g1 + labs(fill=legend)
    }
    if (!is.null(title)) {
        g1 <- g1 + ggtitle(title)
    }
    invisible(list(df=xyz, graph=g1))
}
                     
