ggimage <- function(zmat, x=NULL, y=NULL, col=viridis::viridis(256),
                    xlab=NULL, ylab=NULL, legend=NULL, title=NULL,
                    xrev=FALSE, yrev=FALSE, asp=1,
                    addcontour=FALSE, binwidth=NULL
                   ) {
    require(ggplot2)
    if (is.list(zmat)) {
      zmat <- zmat$zmat
      x <- zmat$ra
      y <- zmat$dec
      xrev <- TRUE
    }
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
    if (addcontour) {
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
    plot(g1)
    invisible(list(df=df, graph=g1))
}

## simple 2D spatial interpolation

fillpoly <- function(ra, dec, zvals, dxy=0.5, min_ny=0, usefields=TRUE) {
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
    list(zmat=zsurf$z, ra=zsurf$x, dec=zsurf$y)
}

ggbinplot <- function (gdat, z, zlab = NULL, addfiberpos = TRUE, addcentroid = FALSE, 
                addborders = addfiberpos, addfp = FALSE, addcontour = !addfiberpos,
                addbinno = !addfiberpos,
                show.legend = TRUE, na.value = "grey95", 
                palette = "Set1", colors = viridis::viridis(256), 
                contourcolor="black", cbinwidth=25, fpcolor = "red") {
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
  if (addbinno) {
    g1 <- g1 + annotate(geom="text", x=df$x, y=df$y, label=1:length(df$x))
  }
  if (addcontour) {
    allok <- complete.cases(df)
    zsurf <- akima::interp(x=df$x[allok], y=df$y[allok], z=df$z[allok])
    xy <- expand.grid(zsurf$x, zsurf$y)
    df4 <- data.frame(x=xy[,1], y=xy[,2], z=as.vector(zsurf$z))
    g1 <- g1 + metR::geom_contour2(aes(x=x, y=y, z=z, label = after_stat(level)), data=df4,
                            color=contourcolor, binwidth=cbinwidth, label_size=5, na.rm=TRUE)
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
                     
