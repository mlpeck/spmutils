## estimate redshift offsets

getdz <- function(gdat, lib, snrthresh=5, nlthresh=2000, dzlim=0.003, searchinterval=1e-4) {
  nr <- nrow(gdat$flux)
  dz <- numeric(nr)
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
      next
    }
    f <- gdat$flux[i, ]
    if (length(which(!is.na(f))) < nlthresh) {
      dz[i] <- NA
      next
    }
    iv <- gdat$ivar[i, ]
    dev_grid <- fn_v(z_search)
    z0 <- z_search[which.min(dev_grid)]
    bestz <- Rsolnp::solnp(pars=z0, fun=fn, 
                           LB=z0-searchinterval, UB=z0+searchinterval, 
                           control=list(trace=0))
    dz[i] <- bestz$pars
  }
  dz
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

nnfitmanga <- function(gdat, dz, lib.mod,
                       nz=length(Z), nt=length(ages),
                       snrthresh=5, tsf=0.1, rlaw=calzetti, dlogl=1.e-4, 
                       starts = c(0.25, 1., 1.), lb=c(0, 0.7, 0), ub=c(3., 5., 5.),
                       flux.em.bad=1.e5,
                       which.lick=c(1, 13:15, 20),
                       PLOT=TRUE) {
  attach(lib.mod)
  on.exit(detach(lib.mod))
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
    fit.nn$deviance
  }
  
  nr <- length(dz)
#  lib.ssp$lambda <- airtovac(lib.ssp$lambda)
  olib.ssp <- lib.ssp
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
  d4000_n <- rep(NA, nr)
  d4000_n_err <- rep(NA, nr)
  lick <- matrix(NA, nrow=nr, ncol=length(which.lick))
  lick.err <- matrix(NA, nrow=nr, ncol=length(which.lick))
  flux.em <- matrix(NA, nrow=nr, ncol=n.em)
  
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
    log_lik[i] <- sum(dnorm(flux-fitted, sd=1/sqrt(ivar), log=TRUE), na.rm=TRUE)
    
    ## plot spectrum and fit
    tdat <- data.frame(lambda=lambda.rest, gflux=flux, fitted=fitted,
                       residual=residual)
    tlong <- melt(tdat, id.vars="lambda")
    val <- c(rep(" obs",length(lambda.rest)), rep("fitted",length(lambda.rest)),
             rep("residual",length(lambda.rest)))
    tlong <- cbind(tlong, val = val)
    tlong$variable[tlong$variable=="gflux"] <- "fitted"
    base <- ggplot(tlong, aes(x=lambda, y=value, color=val)) + geom_line() + xlab(expression(lambda)) +
                    ylab(NULL) + ggtitle(paste("(",i,") log_lik= ", 
                                        format(log_lik[i], digits=1), sep=""))
    add.resid <- facet_grid(variable ~ ., scale="free_y")
    plot(base+add.resid)
    
    ## emission line fluxes and errors
    flux.em[i,in.em] <- b.em*em.mult[in.em]

    ##d4000 and lick indices
    d4 <- d4000n(lambda.rest, gflux.net, ivar)
    d4000_n[i] <- d4$d4000_n
    d4000_n_err[i] <- d4$d4000_n_err
    
    d4 <- lickew(lambda.rest, gflux.net, ivar, which.index=which.lick)
    lick[i,] <- d4[1:length(which.lick)]
    lick.err[i,] <- d4[(length(which.lick)+1):length(d4)]            
  }    
  
  flux.em[flux.em > flux.em.bad] <- NA
  dimnames(lick)[[2]] <- as.list(names(d4)[1:length(which.lick)])
  dimnames(lick.err)[[2]] <- as.list(names(d4)[(length(which.lick)+1):length(d4)])
  dimnames(flux.em)[[2]] <- as.list(names(lambda.em))
  options("warn"=0)
  returns <- list(tauv=tauv, vdisp.em=vdisp.em, vdisp.st=vdisp.st, 
                  d4000_n=d4000_n, d4000_n_err=d4000_n_err, 
                  lick=lick, lick.err=lick.err, 
                  flux.em=flux.em,
                  log_lik=log_lik, nnfits=nnfits)
  returns
}
                       
            
