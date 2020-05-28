## stan model specific utilities to wrangle data and provide initial values to Stan

## normalize galaxy and stellar library fluxes by mean over all good data points

prep_data_avg <- function(gdat, dz, nnfits, which.spax) {
  i <- which.spax
  z <- gdat$meta$z+dz[i]
  flux <- gdat$flux[i, ]
  ivar <- gdat$ivar[i, ]
  lambda.rest <- gdat$lambda/(1+z)
  logl <- log10(lambda.rest)
  T.gyr <- 10^(ages-9)
  dT <- diff(c(0, T.gyr))
  
  lib.ssp$lambda <- airtovac(lib.ssp$lambda)
  lib.st <- regrid(lambda.rest, lib.ssp)
  x.st <- blur.lib(lib.st, nnfits$vdisp.st[i])
  n.st <- ncol(x.st)
  nz <- length(Z)
  nt <- n.st/nz
  ind.young <- (0:(nz-1))*nt+1
  allok <- complete.cases(flux, ivar, x.st)
  nl <- length(which(allok))
  norm_g <- mean(pmax(flux[allok], 0))
  norm_st <- colMeans(x.st[allok,])
  x.st <- scale(x.st[allok,], center=FALSE, scale=norm_st)
  
  x.em <- make_emlib(emlines, nnfits$vdisp.em[i], logl, allok)
  in.em <- x.em$in_em
  x.em <- x.em$x_em
  n.em <- ncol(x.em)
  norm_em <- 1
  
  list(nt=nt, nz=nz, nl=nl, n_em=n.em,
       ind_young=ind.young,
       lambda=lambda.rest[allok],
       gflux=flux[allok]/norm_g,
       g_std=1/sqrt(ivar[allok])/norm_g,
       norm_g=norm_g,
       norm_st=norm_st,
       norm_em=norm_em,
       sp_st=x.st,
       dT=rep(dT, nz),
       sp_em=x.em[allok,],
       in_em=in.em
      )
}

## normalize by mean over wavelength interval (5375, 5625)

prep_data_mod <- function (gdat, dz, nnfits, which.spax) {
  i <- which.spax
  z <- gdat$meta$z + dz[i]
  flux <- gdat$flux[i, ]
  ivar <- gdat$ivar[i, ]
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

## initialize Stan's optimizer to estimate MAP solution

init_opt_mod <- function(spm_data, nnfits, which.spax, jv) {
  n_st <- ncol(spm_data$sp_st)
  n_em <- spm_data$n_em
  norm_st <- spm_data$norm_st
  norm_em <- spm_data$norm_em
  norm_g <- spm_data$norm_g
  b <- nnfits$nnfits[which.spax,]
  b_st_s <- b[1:n_st] * norm_st/norm_g + runif(n_st, min=jv/10, max=jv)
  a <- sum(b_st_s)
  b_st_s <- b_st_s/a
  b_em <- b[(n_st+1):(n_st+n_em)] * norm_em/norm_g + runif(n_em, min=jv/10, max=jv)
  tauv <- nnfits$tauv[which.spax]
  if (tauv == 0) tauv=runif(1, min=jv/10, max=jv)
  delta <- rnorm(1, 0, jv)
  list(a=a, b_st_s=b_st_s, b_em=b_em, tauv=tauv, delta=delta)
}

## initialize Stan's sampler

init_sampler_mod <- function(X, stan_opt, jv) {
  a <- as.numeric(stan_opt$a + rnorm(1, sd=jv))
  b_st_s <- stan_opt$b_st_s + runif(length(stan_opt$b_st_s), min=jv/10, max=jv)
  b_st_s <- b_st_s/sum(b_st_s)
  b_em <- stan_opt$b_em + runif(length(stan_opt$b_em), min=jv/10, max=jv)
  tauv <- as.numeric(stan_opt$tauv + runif(1, min=jv/10, max=jv))
  delta <- as.numeric(stan_opt$delta + rnorm(1, sd=jv))
  list(a=a, b_st_s=b_st_s, b_em=b_em, tauv=tauv, delta=delta)
}

## Stan model specific data to track in stanfit_batch

init_tracked_mod <- function(nsim, n_st, n_em, nr) {
  assign("b_st", array(NA, dim=c(nsim, n_st, nr)), parent.frame())
  assign("b_em", array(NA, dim=c(nsim, n_em, nr)), parent.frame())
  assign("a", matrix(NA, nrow=nsim, ncol=nr), parent.frame())
  assign("tauv", matrix(NA, nrow=nsim, ncol=nr), parent.frame())
  assign("delta", matrix(NA, nrow=nsim, ncol=nr), parent.frame())
  assign("in_em", matrix(NA, nrow=n_em, ncol=nr), parent.frame())
  assign("ll", matrix(NA, nrow=nsim, ncol=nr), parent.frame())
  assign("walltime", rep(NA, nr), parent.frame())
  assign("divergences", rep(NA, nr), parent.frame())
  assign("max_treedepth", rep(NA, nr), parent.frame())
  assign("norm_g", rep(NA, nr), parent.frame())
  assign("norm_st", matrix(NA, nrow=n_st, ncol=nr), parent.frame())
  assign("norm_em", rep(NA, nr), parent.frame())
}

update_tracked_mod <- function(i, sfit, fpart) {
  post <- rstan::extract(sfit$stanfit)
  env <- parent.frame()
  env$b_st[,,i] <- post$b_st
  env$b_em[,sfit$in_em,i] <- post$b_em
  env$a[,i] <- post$a
  env$tauv[,i] <- post$tauv
  env$delta[,i] <- post$delta
  env$in_em[sfit$in_em,i] <- sfit$in_em
  env$ll[,i] <- post$ll
  env$walltime[i] <- max(rowSums(rstan::get_elapsed_time(sfit$stanfit)))
  sp <- rstan::get_sampler_params(sfit$stanfit, inc_warmup=FALSE)
  env$divergences[i] <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  env$max_treedepth[i] <- max(sapply(sp, function(x) max(x[, "treedepth__"])))
  env$norm_g[i] <- sfit$norm_g
  env$norm_st[,i] <- sfit$norm_st
  env$norm_em[i] <- sfit$norm_em
  save(b_st, b_em, tauv, delta, a, in_em, ll,
       walltime, divergences, max_treedepth,
       norm_g, norm_st, norm_em, envir=env, file = fpart)
}

return_tracked_mod <- function() {
  env <- parent.frame()
  retval <- list(b_st= env$b_st, b_em= env$b_em, a= env$a, tauv= env$tauv, delta= env$delta, in_em= env$in_em, ll= env$ll,
       walltime= env$walltime, divergences= env$divergences, max_treedepth= env$max_treedepth,
       norm_g= env$norm_g, norm_st= env$norm_st, norm_em= env$norm_em)
  retval
}

## replace a single stan run in batch fits

replace_sfit_mod <- function(sfit.all, sfit.one, which.spax) {
  post <- rstan::extract(sfit.one$stanfit)
  sp <- rstan::get_sampler_params(sfit.one$stanfit, inc_warmup=FALSE)
  sfit.all$b_st[,,which.spax] <- post$b_st
  sfit.all$b_em[,sfit.one$in.em,which.spax] <- post$b_em
  sfit.all$a[,which.spax] <- post$a
  sfit.all$tauv[,which.spax] <- post$tauv
  sfit.all$delta[,which.spax] <- post$delta
  sfit.all$in_em[sfit.one$in.em,which.spax] <- sfit.one$in.em
  sfit.all$ll[,which.spax] <- post$ll
  sfit.all$walltime[which.spax] <- max(rowSums(rstan::get_elapsed_time(sfit.one$stanfit)))
  sfit.all$divergences[which.spax] <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  sfit.all$max_treedepth[which.spax] <- max(sapply(sp, function(x) max(x[, "treedepth__"])))
  sfit.all$norm_g[which.spax] <- sfit.one$norm_g
  sfit.all$norm_st[,which.spax] <- sfit.one$norm_st
  sfit.all$norm_em[which.spax] <- sfit.one$norm_em
  sfit.all
}


