sum_spec <- function(flux, ivar, indexes, wl.ind, cfn) {
  if (length(indexes)==1) {
    flux <- flux[indexes,]
    ivar <- ivar[indexes,]
  } else {
    flux <- apply(flux[indexes,], 2, sum)
    ivar <- 1/(apply(1/ivar[indexes,], 2, sum))
  }
  snr <- cfn(sqrt(pmax(flux[wl.ind[1]:wl.ind[2]], 0)^2 * ivar[wl.ind[1]:wl.ind[2]]), na.rm=TRUE)
  list(flux=flux, ivar=ivar, snr=snr)
}


bin2d <- function(gdat, snrthresh=6.25, f.fluxok=0.8, wl.limits = c(3700, 7500), cfn=mean, WVT=FALSE, PLOT=TRUE) {
  voronoi_tessellation <- function(x, y, xnode, ynode, scale) {
    if (scale[1] == 1) {
      classe <- nabor::knn(cbind(xnode, ynode), query=cbind(x, y), k=1)
      classe <- as.vector(classe$nn.idx)
    } else {
      classe <- apply((outer(x, xnode, '-')^2 + outer(y, ynode, '-')^2)/scale^2, 1, which.min)
    }
    classe
  }
  
  
  accretion <- function(flux, ivar, x, y, snr, snrthresh=5, wl.ind, cfn) {
    
    n <- length(x)
    classe <- numeric(n)
    good <- numeric(n)
    
    cbin <- which.max(snr)
    SN <- snr[cbin]
    
    
    for (ind in 1:n) {
      classe[cbin] <- ind
      xbar <- x[cbin]
      ybar <- y[cbin]
      
      repeat {
        if (all(classe>0)) break
          unbinned <- which(classe == 0)
          k <- which.min((x[unbinned]-xbar)^2+(y[unbinned]-ybar)^2)
          mindist <- sqrt(min((x[cbin] - x[unbinned[k]])^2 + (y[cbin]-y[unbinned[k]])^2))
          
          nextbin <- c(cbin, unbinned[k])
          
          SN.old <- SN
          binspec <- sum_spec(flux, ivar, nextbin, wl.ind=wl.ind, cfn=cfn)
          SN <- binspec$snr
          if (abs(SN-snrthresh) > abs(SN.old-snrthresh) || SN.old > SN) {
            if (SN.old > 0.8*snrthresh) {
              good[cbin] <- 1
              break
            }
          }
          classe[unbinned[k]] <- ind
          cbin <- nextbin
          good[cbin] <- 1
          
          xbar <- mean(x[cbin])
          ybar <- mean(y[cbin])
      }
      binned <- which(classe > 0)
      if (all(classe > 0)) break
        unbinned <- setdiff(1:n, binned)
        xbar <- mean(x[binned])
        ybar <- mean(y[binned])
        k <- which.min((x[unbinned]-xbar)^2+(y[unbinned]-ybar)^2)
        cbin <- unbinned[k]
        SN <- snr[unbinned[k]]
    }
    
    classe <- classe*good
    classe
  }
  
  reassign_bad_bins <- function(classe, x, y) {
    good <- unique(classe[classe>0])
    xnode <- numeric(length(good))
    ynode <- numeric(length(good))
    
    for (i in seq_along(good)) {
      xnode[i] <- mean(x[which(classe==good[i])])
      ynode[i] <- mean(y[which(classe==good[i])])
    }
    bad <- which(classe == 0)
    if (length(bad) > 0) {
      index <- voronoi_tessellation(x[bad], y[bad], xnode, ynode, scale=1)
      classe[bad] <- good[index]
      good <- unique(classe)
      for (i in seq_along(good)) {
        xnode[i] <- mean(x[which(classe==good[i])])
        ynode[i] <- mean(y[which(classe==good[i])])
      }
    }
    list(classe=classe, xnode=xnode, ynode=ynode)
  }
  
  cvt_equal_mass <- function(flux, ivar, x, y, snr, xnode, ynode, wl.ind, cfn, WVT) {
    dens2 <- snr^4
    scale <- rep(1, length(xnode))
    
    for (it in seq_along(xnode)) {
      xnode_old <- xnode
      ynode_old <- ynode
      classe <- voronoi_tessellation(x, y, xnode, ynode, scale)
      bins <- unique(classe)
      for (k in bins) {
        ind <- which(classe == k)
        if (WVT) {
          xnode[k] <- mean(x[ind])
          ynode[k] <- mean(y[ind])
          sn <- sum_spec(flux, ivar, ind, wl.ind=wl.ind, cfn=cfn)$snr
          scale[k] <- sqrt(length(ind)/sn)
        } else {
          mass <- sum(dens2[ind])
          xnode[k] <- sum(x[ind]*dens2[ind])/mass
          ynode[k] <- sum(y[ind]*dens2[ind])/mass
        }
      }
      diff <- sqrt(sum((xnode-xnode_old)^2 + (ynode-ynode_old)^2))
      if (diff < 0.1) break
    }
    if (diff > 0) {
      classe <- voronoi_tessellation(x, y, xnode, ynode, scale)
      bins <- unique(classe)
    }
    list(classe=classe, xnode = xnode[bins], ynode=ynode[bins], scale=scale[bins])
  }
 
 ## start of main routine
  
  if (is.null(wl.limits)) {
    wl.ind <- c(1, length(gdat$lambda))
  }
  else {
    wl.ind <- findInterval(wl.limits, gdat$lambda/(1+gdat$meta$z), all.inside=TRUE)
  }
  
  fluxok <- !is.na(gdat$flux[,,wl.ind[1]:wl.ind[2]])
  p.fluxok <- apply(fluxok, 1, sum)/(wl.ind[2]-wl.ind[1]+1)
  ok2bin <- (p.fluxok >= f.fluxok)
  
  x <- gdat$xpos[ok2bin]
  y <- gdat$ypos[ok2bin]
  flux <- gdat$flux[ok2bin,,]
  ivar <- gdat$ivar[ok2bin,,]
  ra.f <- gdat$ra.f[ok2bin]
  dec.f <- gdat$dec.f[ok2bin]
  
  snr <- numeric(length(x))
  for (i in seq_along(x)) {
    snr[i] <- cfn(sqrt(pmax(flux[i,wl.ind[1]:wl.ind[2]], 0)^2 * ivar[i,wl.ind[1]:wl.ind[2]]), na.rm=TRUE)
  }
  
  classe <- accretion(flux, ivar, x, y, snr, snrthresh=snrthresh, wl.ind=wl.ind, cfn=cfn)
  starts <- reassign_bad_bins(classe, x, y)
  final <- cvt_equal_mass(flux, ivar, x, y, snr, starts$xnode, starts$ynode, wl.ind=wl.ind, cfn=cfn, WVT=WVT)
  classe <- final$classe
  bins <- unique(classe)
  nbins <- length(bins)
  fibersinbin <- numeric(nbins)
  flux.b <- array(0, dim=c(nbins, 1, length(gdat$lambda)))
  ivar.b <- array(0, dim=c(nbins, 1, length(gdat$lambda)))
  snr <- matrix(0, nbins, 1)
  
  xpos <- final$xnode
  ypos <- final$ynode
  dec.f <- gdat$meta$dec + ypos/3600
  ra.f <- gdat$meta$ra - xpos/3600/cos(dec.f*pi/180)

  for (i in 1:nbins) {
    spec <- sum_spec(flux, ivar, which(classe==bins[i]), wl.ind=wl.ind, cfn=cfn)
    fibersinbin[i] <- length(which(classe==bins[i]))
    flux.b[i,1,] <- spec$flux
    ivar.b[i,1,] <- spec$ivar
    snr[i, 1] <- spec$snr
  }
  bin.fiber <- rep(NA, length(ok2bin))
  bin.fiber[ok2bin] <- classe
  gdat.bin <- list(meta=gdat$meta, lambda=gdat$lambda, flux=flux.b, ivar=ivar.b, snr=snr, 
                   xpos=xpos, ypos=ypos, ra.f=ra.f, dec.f=dec.f,
                   bin.fiber=bin.fiber, fibersinbin=fibersinbin, 
                   x.orig = gdat$xpos, y.orig = gdat$ypos)
  if (PLOT) {
    ggbinplot(gdat.bin, unique(gdat.bin$bin.fiber[!is.na(gdat.bin$bin.fiber)]), addfiberpos=TRUE, addcentroid=TRUE,
              addborders=TRUE, show.legend=FALSE)
  }
  gdat.bin
}
