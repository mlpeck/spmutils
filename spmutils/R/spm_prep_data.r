## stan model specific utilities to wrangle data and provide initial values to Stan

prep_data_mod <- function(nnfits, which.spax) {
    lib.ssp$lambda <- airtovac(lib.ssp$lambda)
    lib.st <- regrid(lambda.rest, lib.ssp)
    x.st <- blur.lib(lib.st, nnfits$vdisp.st[which.spax])
    n.st <- ncol(x.st)
    nz <- length(Z)
    nt <- n.st/nz
    ind.young <- (0:(nz-1))*nt+1
    allok <- complete.cases(flux, ivar, x.st)
    nl <- length(which(allok))
    norm_g <- mean(flux[allok])
    norm_st <- colMeans(x.st[allok,])
    x.st <- scale(x.st[allok,], center=FALSE, scale=norm_st)
    
    x.em <- make_emlib(emlines, nnfits$vdisp.em[which.spax], logl, allok)
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
         sp_em=x.em[allok,])
}

init_opt_mod <- function(nnfits, which.spax, jv) {
    b <- nnfits$nnfits[which.spax,]
    b_st_s <- b[1:n.st] * norm_st/norm_g + runif(n.st, min=jv/10, max=jv)
    a <- sum(b_st_s)
    b_st_s <- b_st_s/a
    b_em <- b[(n.st+1):(n.st+n.em)] * norm_em/norm_g + runif(n.em, min=jv/10, max=jv)
    tauv <- nnfits$tauv[which.spax]
    if (tauv == 0) tauv=runif(1, min=jv/10, max=jv)
    delta <- rnorm(1, 0, jv)
    inits <- list(a=a, b_st_s=b_st_s, b_em=b_em, tauv=tauv, delta=delta)
  
}

init_sampler_mod <- function(X, stan_opt, jv) {
  a <- stan_opt$a + rnorm(1, sd=jv)
  b_st_s <- stan_opt$b_st_s + runif(length(stan_opt$b_st_s), min=jv/10, max=jv)
  b_st_s <- b_st_s/sum(b_st_s)
  b_em <- stan_opt$b_em + runif(length(stan_opt$b_em), min=jv/10, max=jv)
  tauv <- stan_opt$tauv + runif(1, min=jv/10, max=jv)
  delta <- stan_opt$delat + rnorm(1, sd=jv)
  list(a=a, b_st_s=b_st_s, b_em=b_em, tauv=tauv, delta=delta)
}

## Stan model specific data to track in stanfit_batch

init_tracked_mod <- function() {
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

update_tracked_mod <- function() {
  post <- extract(sfit$stanfit)
  b_st[,,i] <<- post$b_st
  b_em[,sfit$in.em,i] <<- post$b_em
  a[,i] <<- post$a
  tauv[,i] <<- post$tauv
  delta[,i] <<- post$delta
  in_em[sfit$in.em,i] <<- sfit$in.em
  ll[,i] <<- post$ll
  walltime[i] <<- max(rowSums(get_elapsed_time(sfit$stanfit)))
  sp <<- get_sampler_params(sfit$stanfit, inc_warmup=FALSE)
  divergences[i] <<- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  max_treedepth[i] <<- max(sapply(sp, function(x) max(x[, "treedepth__"])))
  norm_g[i] <<- sfit$norm_g
  norm_st[,i] <<- sfit$norm_st
  norm_em[i] <<- sfit$norm_em
  save(b_st, b_em, tauv, a, in_em, ll,
       walltime, divergences, max_treedepth,
       norm_g, norm_st, norm_em, file = fpart)
}

return_tracked_mod <- function() {
  list(b_st=b_st, b_em=b_em, a=a, tauv=tauv, delta=delta, in_em=in_em, ll=ll,
       walltime=walltime, divergences=divergences, max_treedepth=max_treedepth,
       norm_g=norm_g, norm_st=norm_st, norm_em=norm_em)
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


