## utilities for manipulating positions

radec2hms <- function(...) {
  radec <- unlist(list(...))
  dec <- radec[2]
  ra <- radec[1]/15
  ra.h <- trunc(ra)
  ra <- 60*(ra-ra.h)
  ra.m <- trunc(ra)
  ra.s <- 60*(ra-ra.m)
  ds <- sign(dec)
  dec.d <- trunc(dec)
  dec <- 60*abs(dec-dec.d)
  dec.m <- trunc(dec)
  dec.s <- 60*(dec-dec.m)
  cat(paste(ra.h, ra.m, format(ra.s, digits=2),sep=":"))
  if (ds >= 0) cat(" +") else cat(" -")
  cat(paste(dec.d, dec.m, format(dec.s, digits=1), sep=":"), "\n")
}

hms2radec <- function(...) {
  radec <- unlist(list(...))
  nr <- length(radec)/2
  radec <- matrix(radec,nr,2)
  ra <- radec[,1]
  dec <- radec[,2]
  ra.h <- trunc(ra)
  ra <- 100*(ra-ra.h)
  ra.m <- trunc(ra)
  ra.s <- 100*(ra-ra.m)
  ra.d <- 15*(ra.h + ra.m/60 + ra.s/3600)
  dec.sign <- sign(dec)
  dec <- abs(dec)
  dec.d <- trunc(dec)
  dec <- 100*abs(dec-dec.d)
  dec.m <- trunc(dec)
  dec.s <- 100*(dec-dec.m)
  dec.d <- dec.sign*(dec.d+dec.m/60+dec.s/3600)
  list(ra=ra.d, dec=dec.d)
}


lb2radec <- function(...) {
  rotm <- matrix(c(-.054876, -.873437, -.483835,
        .494109, -0.444830, .746982,
        -0.867666, -0.198076, 0.455984), 3, 3)
  lb <- unlist(list(...))
  nr <- length(lb)/2
  lb <- matrix(lb, nr, 2)
  l <- pi*lb[,1]/180
  b <- pi*lb[,2]/180
  x.lb <- rbind(cos(b)*cos(l), cos(b)*sin(l), sin(b))
  x.ad <- rotm %*% x.lb
  a <- atan2(x.ad[2,], x.ad[1,])
  a[a<0] <- 2*pi+a[a<0]
  d <- asin(x.ad[3,])
  ad <- data.frame(ra=a*180/pi, dec=d*180/pi)
  ad
}

radec2lb <- function(...) {
    rotm <- t(matrix(c(-.054876, -.873437, -.483835,
        .494109, -0.444830, .746982,
        -0.867666, -0.198076, 0.455984), 3, 3))
    radec <- unlist(list(...))
    nr <- length(radec)/2
    radec <- matrix(radec,nr,2)
    ra <- pi*radec[,1]/180
    dec <- pi*radec[,2]/180
    x.ad <- rbind(cos(ra)*cos(dec),sin(ra)*cos(dec),sin(dec))
    x.lb <- rotm %*% x.ad
    l <- atan2(x.lb[2,], x.lb[1,])
    l[l<0] <- 2*pi+l[l<0]
    b <- asin(x.lb[3,])
    lb <- data.frame(l=180*l/pi, b=180*b/pi)
    lb
}

