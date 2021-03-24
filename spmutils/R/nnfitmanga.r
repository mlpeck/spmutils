## estimate redshift offsets

getdz <- function(gdat, lib, snrthresh=5, nlthresh=2000, dzlim=0.003, searchinterval=1e-4) {
  nr <- nrow(gdat$flux)
  dz <- numeric(nr)
  dz.err <- numeric(nr)
  z <- gdat$meta$z
  lambda <- gdat$lambda
  fmla <- as.formula(paste("f ~", paste(names(lib)[-1], collapse="+")))
  fn <- function(x) {
    lambda.out <- lambda/(1+z+x)
    lib.out <- regrid(lambda.out, lib)
    lfit <- lm(fmla, data=lib.out, weights=iv)
    deviance(lfit)
  }
  fn_v <- Vectorize(fn)
  z_search <- seq(-dzlim, dzlim, by=searchinterval)
  
  for (i in 1:nr) {
    if (is.na(gdat$snr[i]) || gdat$snr[i]<=snrthresh) {
      dz[i] <- NA
      dz.err[i] <- NA
      next
    }
    f <- gdat$flux[i, ]
    if (length(which(!is.na(f))) < nlthresh) {
      dz[i] <- NA
      dz.err[i] <- NA
      next
    }
    iv <- gdat$ivar[i, ]
    dev_grid <- fn_v(z_search)
    z0 <- z_search[which.min(dev_grid)]
    bestz <- Rsolnp::solnp(pars=z0, fun=fn, 
                           LB=z0-searchinterval, UB=z0+searchinterval, 
                           control=list(trace=0))
    dz[i] <- bestz$pars
    dz.err[i] <- sqrt(2/bestz$hessian)
  }
  list(dz=dz, dz.err=dz.err)
}

## Moore-Penrose pseudo inverse (needed for uncertainty estimates in fit.nn)

mpinv <- function(X) {
    S <- svd(X)
    eps <- .Machine$double.eps * max(dim(X)) * S$d[1]
    dinv <- numeric(length(S$d))
    dinv[S$d >= eps] <- 1/S$d[S$d >= eps]
    tcrossprod(S$v %*% diag(dinv), S$u)
}


## nnls fits to manga data cube or rss file

nnfitmanga <- function(gdat, dz,
                       nz=length(Z), nt=length(ages),
                       snrthresh=5, tsf=0.1, rlaw=calzetti, dlogl=1.e-4, 
                       starts = c(0.25, 1., 1.), lb=c(0, 0.7, 0), ub=c(3., 5., 5.),
                       flux.em.bad=1.e5,
                       which.lick=c(1, 13:15, 20),
                       PLOT=TRUE) {
  require(nnls)
  require(cosmo)
  if (PLOT) {
    require(ggplot2)
    require(reshape2)
    options("warn" = -1)
  }
  
  ## variables and function definition for nonlinear fit to get tau, vdisp.st, vdisp.em
  
  fit.nn <- NULL
  x.st <- NULL
  x.em <- NULL
  in.em <- NULL
  allok <- NULL
  fn <- function(pars) {
    tauv <- pars[1]
    vdisp.em <- 100*pars[2]
    vdisp.st <- 100*pars[3]
    x.st <<- blur.lib(lib.ssp, vdisp.st)
    allok <<- complete.cases(flux, ivar, x.st)
    lib.em <- make_emlib(lambda.em, vdisp.em, logl, allok)
    x.em <<- lib.em$x_em
    in.em <<- lib.em$in_em
    att <- as.vector(rlaw(lambda.rest[allok], tauv))
    if (is.null(in.em)) {
      fit.nn <<- nnls(x.st[allok,]*att*sqrt(ivar[allok]), flux[allok]*sqrt(ivar[allok]))
    } else {
      fit.nn <<- nnls(cbind(x.st[allok,]*att, x.em[allok, ])*sqrt(ivar[allok]),
                      flux[allok]*sqrt(ivar[allok]))
    }
    fit.nn$deviance/2
  }
  
  dz <- dz$dz
  nr <- length(dz)
  lib.ssp$lambda <- airtovac(lib.ssp$lambda)
  olib.ssp <- lib.ssp
  tauv_err <- rep(NA, nr)
  vdisp.em_err <- rep(NA, nr)
  vdisp.st_err <- rep(NA, nr)
  tauv <- rep(NA, nr)
  vdisp.em <- rep(NA, nr)
  vdisp.st <- rep(NA, nr)
  
  T.gyr <- 10^(ages-9)
  isf <- which.min(abs(tsf-T.gyr))
  lambda.em <- lambda_em
  n.em <- length(lambda.em)
  ##needed to correct for log lambda grid
  em.mult <- lambda.em * log(10)/10000
  n.st <- ncol(lib.ssp)-1
  
  nnfits <- matrix(NA, nrow=nr, ncol=n.em+n.st)
  log_lik <- rep(NA, nr)
  tbar <- rep(NA, nr)
  tbar.lum <- rep(NA, nr)
  zbar <- rep(NA, nr)
  Mstar <- rep(NA, nr)
  gri <- matrix(NA, nrow=nr, ncol=nrow(gri.ssp))
  d4000_n <- rep(NA, nr)
  d4000_n_err <- rep(NA, nr)
  lick <- matrix(NA, nrow=nr, ncol=length(which.lick))
  lick.err <- matrix(NA, nrow=nr, ncol=length(which.lick))
  flux.em <- matrix(NA, nrow=nr, ncol=n.em)
  flux.em.err <- matrix(NA, nrow=nr, ncol=n.em)
  
  for (i in 1:nr) {
    if (is.na(gdat$snr[i]) || gdat$snr[i]<=snrthresh || is.na(dz[i])) {
      next
    }
    lambda <- gdat$lambda
    flux <- gdat$flux[i,]
    ivar <- gdat$ivar[i,]
    z <- gdat$meta$z+dz[i]
    lambda.rest <- lambda/(1+z)
    logl <- log10(lambda.rest)
    lib.ssp <- regrid(lambda.rest, olib.ssp)
    fitij <- Rsolnp::solnp(pars=starts, fn, LB=lb, UB=ub)
    tauv[i] <- fitij$pars[1]
    vdisp.em[i] <- 100*fitij$pars[2]
    vdisp.st[i] <- 100*fitij$pars[3]
    errs <- tryCatch(sqrt(diag(solve(fitij$hessian))), error=function(e) rep(NA, 3))
    tauv_err[i] <- errs[1]
    vdisp.em_err[i] <- 100*errs[2]
    vdisp.st_err[i] <- 100*errs[3]
    
    b <- fit.nn$x
    b.st <- b[1:n.st]
    b.em <- b[(n.st+1):length(b)]
    nnfits[i, 1:length(b)] <- b
    ni.em <- length(b)-n.st
    nl <- length(lambda.rest[allok])
    fitted <- cbind(x.st*as.vector(rlaw(lambda.rest,tauv[i])), x.em) %*% b
    fitted.em <- x.em %*% b.em
    gflux.net <- flux-fitted.em
    residual <- (flux-fitted)*sqrt(ivar)
    
    ## basic measures of goodness of fit            
    log_lik[i] <- fit.nn$deviance/2
    
    ## plot spectrum and fit
    tdat <- data.frame(lambda=lambda.rest, gflux=flux, fitted=fitted,
                       residual=residual)
    tlong <- melt(tdat, id.vars="lambda")
    val <- c(rep(" obs",length(lambda.rest)), rep("fitted",length(lambda.rest)),
             rep("residual",length(lambda.rest)))
    tlong <- cbind(tlong, val = val)
    tlong$variable[tlong$variable=="gflux"] <- "fitted"
    base <- qplot(lambda, value, data=tlong, geom="line", xlab=expression(lambda), 
                  ylab="", col=val, 
                  main=paste("(",i,") log_lik= ", 
                             format(log_lik[i], digits=0), sep=""))
    add.resid <- facet_grid(variable ~ ., scale="free_y")
    plot(base+add.resid)
    
    ## emission line fluxes and errors
    flux.em[i,in.em] <- b.em*em.mult[in.em]
    nzeros <- sort(fit.nn$passive)
    X <- cbind(x.st*as.vector(rlaw(lambda.rest,tauv[i])), x.em)[allok,nzeros]
    V <- max(log_lik[i]*2/nl, 1)*mpinv(crossprod(X, diag(ivar[allok])) %*% X)
    sd.b <- rep(NA, n.st+ni.em)
    sd.b[nzeros] <- sqrt(diag(V))
    flux.em.err[i,in.em] <- sd.b[(n.st+1):(n.st+ni.em)]*em.mult[in.em]
    
    ## gri magnitudes
    gri[i,] <- -2.5*log10(gri.ssp %*% b.st) + 20.092
    
    ## summaries from estimated sfh
    tbar[i] <- log10(sum(rep(T.gyr, nz)*b.st)/sum(b.st))+9
    tbar.lum[i] <- log10(sum(rep(T.gyr, nz) * b.st * gri.ssp["r",])/
    sum(b.st * gri.ssp["r",]))+9
    m.st <- cosmo::lum.sol(1, z)*b.st
    Mstar[i] <- log10(sum(m.st*mstar))
    
    ##d4000 and lick indices
    d4 <- d4000n(lambda.rest, gflux.net, ivar)
    d4000_n[i] <- d4$d4000_n
    d4000_n_err[i] <- d4$d4000_n_err
    
    d4 <- lickew(lambda.rest, gflux.net, ivar, which.index=which.lick)
    lick[i,] <- d4[1:length(which.lick)]
    lick.err[i,] <- d4[(length(which.lick)+1):length(d4)]            
  }    
  
  flux.em[flux.em > flux.em.bad] <- NA
  flux.em.err[!is.finite(flux.em)] <- NA
  dimnames(lick)[[2]] <- as.list(names(d4)[1:length(which.lick)])
  dimnames(lick.err)[[2]] <- as.list(names(d4)[(length(which.lick)+1):length(d4)])
  dimnames(gri)[[2]] <- dimnames(gri.ssp)[[1]]
  dimnames(flux.em)[[2]] <- as.list(names(lambda.em))
  dimnames(flux.em.err)[[2]] <- as.list(paste(names(lambda.em), "_err", sep=""))
  options("warn"=0)
  returns <- list(tauv=tauv, vdisp.em=vdisp.em, vdisp.st=vdisp.st, 
                  tauv_err=tauv_err, vdisp.em_err=vdisp.em_err, vdisp.st_err=vdisp.st_err,
                  d4000_n=d4000_n, d4000_n_err=d4000_n_err, 
                  lick=lick, lick.err=lick.err, 
                  flux.em=flux.em, flux.em.err=flux.em.err, gri=gri, 
                  tbar=tbar, tbar.lum=tbar.lum, Mstar=Mstar,
                  log_lik=log_lik, nnfits=nnfits)
  returns
}
                       
            
replot <- function(gdat, dz, nnfit, lib.ssp, 
                   which.spax=gdat$meta$cpix, rlaw=calzetti,
                   title=NULL) {
    require(ggplot2)
    require(reshape2)
    options("warn" = -1)
    
    i <- which.spax[1]
    j <- which.spax[2]
    flux <- gdat$flux[i, j, ]
    ivar <- gdat$ivar[i, j, ]
    lambda <- gdat$lambda
    lambda.rest <- lambda/(1+gdat$meta$z+dz[i,j])
    logl <- log10(lambda.rest)
    lib.ssp$lambda <- airtovac(lib.ssp$lambda)
    lib.ssp <- regrid(lambda.rest, lib.ssp)
    x.st <- blur.lib(lib.ssp, nnfit$vdisp.st[i,j])
    n.st <- ncol(x.st)
    allok <- complete.cases(flux, ivar, x.st)
    x.st <- x.st[allok, ]
    lambda.rest <- lambda.rest[allok]
    lambda.em <- lambda_em
    x.em <- make_emlib(lambda.em, nnfit$vdisp.em[i,j], allok, logl)
    in.em <- x.em$in.em
    x.em <- x.em$x.em
    b <- nnfit$nnfits[i, j, ]
    b <- b[!is.na(b)]
    b.st <- b[1:n.st]
    b.em <- b[(n.st+1):length(b)]
    fitted <- cbind(x.st*as.vector(rlaw(lambda.rest,nnfit$tauv[i,j])), x.em[allok,]) %*% b
    fitted.em <- x.em[allok,] %*% b.em
    gflux.net <- flux[allok]-fitted.em
    residual <- (flux[allok]-fitted)*sqrt(ivar[allok])
    tdat <- data.frame(lambda=lambda.rest, flux=flux[allok], fitted=fitted,
                       fitted.em=fitted.em, residual=residual)
    tlong <- melt(tdat, id.vars="lambda")
    val <- c(rep(" obs",length(lambda.rest)), rep("fitted",length(lambda.rest)),
             rep("em", length(lambda.rest)), rep("residual",length(lambda.rest)))
    tlong <- cbind(tlong, val = val)
    tlong$variable[tlong$variable !="residual"] <- "flux"
    base <- qplot(lambda, value, data=tlong, geom="line", xlab=expression(lambda), 
                  ylab="", col=val)
    if (!is.null(title)) {
        base <- base + ggtitle(title)
    }
    add.resid <- facet_grid(variable ~ ., scale="free_y")
    g1 <- base+add.resid
    plot(g1)
    options("warn" = 0)
    g1
}

getbpt <- function(nnfits, snthresh=3, PLOT=TRUE) {
    flux.em <- nnfits$flux.em
    flux.em.err <- nnfits$flux.em.err
    dims <- dim(flux.em)
    nr <- dims[1]
    nc <- dims[2]
    
    n2halpha <- as.vector(log10(flux.em[,,'nii_6584']/flux.em[,,'h_alpha']))
    sn.halpha <- as.vector(flux.em[,,'h_alpha']/flux.em.err[,,'h_alpha_err'])
    sn.nii <- as.vector(flux.em[,,'nii_6584']/flux.em.err[,,'nii_6584_err'])
    n2halpha[!is.finite(n2halpha) | sn.halpha<snthresh | sn.nii<snthresh] <- NA
    
    oiiihbeta <- as.vector(log10(flux.em[,,'oiii_5007']/flux.em[,,'h_beta']))
    sn.hbeta <- as.vector(flux.em[,,'h_beta']/flux.em.err[,,'h_beta_err'])
    sn.oiii <- as.vector(flux.em[,,'oiii_5007']/flux.em.err[,,'oiii_5007_err'])
    oiiihbeta[!is.finite(oiiihbeta) | sn.hbeta<snthresh | sn.oiii<snthresh] <- NA
    
    df <- data.frame(n2halpha=n2halpha, o3hbeta=oiiihbeta,
                     oii_flux=as.vector(flux.em[,,'oii_3729']), 
                     oii_flux_err=as.vector(flux.em.err[,,'oii_3729_err']),
                     h_beta_flux=as.vector(flux.em[,,'h_beta']),
                     h_beta_flux_err=as.vector(flux.em.err[,,'h_beta_err']),
                     oiii_flux=as.vector(flux.em[,,'oiii_5007']),
                     oiii_flux_err=as.vector(flux.em.err[,,'oiii_5007_err']),
                     h_alpha_flux=as.vector(flux.em[,,'h_alpha']),
                     h_alpha_flux_err=as.vector(flux.em.err[,,'h_alpha_err']),
                     nii_6584_flux=as.vector(flux.em[,,'nii_6584']),
                     nii_6584_flux_err=as.vector(flux.em.err[,,'nii_6584_err']))
    bpt <- emclass(df, snthresh=snthresh, PLOT=PLOT)
    bpt[is.na(nnfits$log_lik)] <- NA
    n2halpha <- matrix(n2halpha, nr, nc)
    o3hbeta <- matrix(oiiihbeta, nr, nc)
    list(bpt=bpt, n2halpha=n2halpha, o3hbeta=o3hbeta)
}

