prep_data_mod_sdss <- function (gdat, nnfits, which.spax=1) {
  i <- which.spax
  z <- gdat$meta$z
  flux <- gdat$flux
  ivar <- gdat$ivar
  lambda.rest <- gdat$lambda/(1 + z)
  logl <- log10(lambda.rest)
  T.gyr <- 10^(ages - 9)
  dT <- diff(c(0, T.gyr))
  lib.ssp$lambda <- airtovac(lib.ssp$lambda)
  lib.st <- regrid(lambda.rest, lib.ssp)
  x.st <- blur.lib(lib.st, nnfits$vdisp.st[i])
  n.st <- ncol(x.st)
  nz <- length(Z)
  nt <- n.st/nz
  ind.young <- (0:(nz - 1)) * nt + 1
  allok <- complete.cases(flux, ivar, x.st)
  lambda <- lambda.rest[allok]
  gflux <- flux[allok]
  x.st <- x.st[allok,]
  wl_ind <- findInterval(c(5375,5625), lambda)
  nl <- length(which(allok))
  norm_g <- mean(pmax(gflux[wl_ind[1]:wl_ind[2]], 0))
  norm_st <- colMeans(x.st[wl_ind[1]:wl_ind[2], ])
  x.st <- scale(x.st, center = FALSE, scale = norm_st)
  x.em <- make_emlib(emlines, nnfits$vdisp.em[i], logl, allok)
  in.em <- x.em$in_em
  x.em <- x.em$x_em
  n.em <- ncol(x.em)
  norm_em <- 1
  list(nt = nt, nz = nz, nl = nl, n_em = n.em, ind_young = ind.young, 
    lambda = lambda, gflux = gflux/norm_g, 
    g_std = 1/sqrt(ivar[allok])/norm_g, norm_g = norm_g, 
    norm_st = norm_st, norm_em = norm_em, sp_st = x.st, 
    dT = rep(dT, nz), sp_em = x.em[allok, ], in_em = in.em)
}

