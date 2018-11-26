nnfitsdss <- function(files, lib.ssp, age, Z, gri.ssp, mstar, 
                      dname="spectra",
                      lambda.em = lambda_em,
                      nz=length(Z), nt=length(age),
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
  
  nr <- length(files)
  lib.ssp$lambda <- airtovac(lib.ssp$lambda)
  olib.ssp <- lib.ssp
  tauv_err <- rep(NA, nr)
  vdisp.em_err <- rep(NA, nr)
  vdisp.st_err <- rep(NA, nr)
  tauv <- rep(NA, nr)
  vdisp.em <- rep(NA, nr)
  vdisp.st <- rep(NA, nr)
  
  T.gyr <- 10^(age-9)
  isf <- which.min(abs(tsf-T.gyr))
  n.em <- length(lambda.em)
  ##needed to correct for log lambda grid
  em.mult <- lambda.em * log(10)
  n.st <- ncol(lib.ssp)-1
  
  nnfits <- matrix(NA, nr, n.em+n.st)
  xi2 <- rep(NA, nr)
  tbar <- rep(NA, nr)
  tbar.lum <- rep(NA, nr)
  zbar <- rep(NA, nr)
  Mstar <- rep(NA, nr)
  gri <- matrix(NA, nr, nrow(gri.ssp))
  d4000_n <- rep(NA, nr)
  d4000_n_err <- rep(NA, nr)
  lick <- matrix(NA, nr, length(which.lick))
  lick.err <- matrix(NA, nr, length(which.lick))
  flux.em <- matrix(NA, nr, n.em)
  flux.em.err <- matrix(NA, nr, n.em)
  
  for (i in 1:nr) {
    gdat <- readsdssfit(files[i], dname=dname, extcor=TRUE, rest=TRUE)
    lambda.rest <- gdat$lambda
    flux <- gdat$flux
    ivar <- gdat$ivar
    z <- gdat$meta$z
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
    xi2[i] <- fit.nn$deviance/(nl-fit.nn$nsetp-1)
    
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
                  main=paste("(",i,") xi2= ", 
                             format(xi2[i], digits=3), sep=""))
    add.resid <- facet_grid(variable ~ ., scale="free_y")
    plot(base+add.resid)
    
    ## emission line fluxes and errors
    flux.em[i,in.em] <- b.em*em.mult[in.em]
    nzeros <- sort(fit.nn$passive)
    X <- cbind(x.st*as.vector(rlaw(lambda.rest,tauv[i])), x.em)[allok,nzeros]
    V <- max(xi2[i], 1)*zernike::mpinv(crossprod(X, diag(ivar[allok])) %*% X)
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
                  xi2=xi2, nnfits=nnfits)
  returns
}
                      
