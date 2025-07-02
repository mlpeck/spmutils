stanfit_one <- function(gdat, nnfits, dz, which.spax,
                        prep_data = prep_data_mod,
                        init_opt = init_opt_mod,
                        init_sampler = init_sampler_mod,
                        stan_model=NULL,
                        stan_file="spm_dust_mod_psum.stan", stan_filedir="~/spmcode/",
                        iter_opt=5000, jv=1.e-4,
                        iter=1000, warmup=250, thin=1, chains=4,
                        threads_per_chain=4,
                        max_treedepth=10, adapt_delta=0.9,
                        open_progress=FALSE, ...) {
    
    require(rstan)
    if (is.null(stan_model)) {
      stan_model <- rstan::stan_model(file=file.path(stan_filedir, stan_file))
    }
    rstan_options(threads_per_chain=threads_per_chain)

    spm_data <- prep_data(gdat, nnfits, dz, which.spax)
    inits <- init_opt(spm_data, nnfits, which.spax, jv)
    spm_opt <- optimizing(stan_model, data=spm_data, init=inits, as_vector=FALSE, verbose=TRUE, iter=iter_opt)
    
    init_pars <- lapply(X=1:chains, init_sampler, stan_opt=spm_opt$par, jv=jv)
    
    stanfit <- sampling(stan_model, data=spm_data,
                     chains=chains, iter=iter, warmup=warmup, thin=thin,
                     cores=getOption("mc.cores"),
                     init=init_pars, open_progress=open_progress,
                     control=list(max_treedepth=max_treedepth, adapt_delta=adapt_delta), ...)
    
    list(spm_data=spm_data, stanfit=stanfit, 
         norm_g=spm_data$norm_g, norm_st=spm_data$norm_st, norm_em=spm_data$norm_em, in_em=spm_data$in_em)
}


                        
stanfit_batch <- function(gdat, nnfits, dz, lib.mod,
                        init_tracked = init_tracked_mod,
                        update_tracked = update_tracked_mod,
                        return_tracked = return_tracked_mod,
                        prep_data = prep_data_mod,
                        init_opt = init_opt_mod,
                        init_sampler = init_sampler_mod,
                        stan_file="spm_dust_mod_psum.stan", stan_filedir="~/spmcode/",
                        iter_opt=5000, jv=1.e-4,
                        iter=1000, warmup=250, chains=4,
                        threads_per_chain=4,
                        max_treedepth=10, adapt_delta=0.9,
                        open_progress=FALSE,
                        start=NULL, end=NULL, fpart="bfits.rda", ...) {
  attach(lib.mod)
  on.exit(detach(lib.mod))
  dims <- dim(gdat$flux)
  nsim <- (iter-warmup)*chains
  nt <- length(lib.mod$ages)
  nr <- dims[1]
  n_st <- ncol(lib.mod$lib.ssp)-1
  n_em <- length(emlines)
  nl <- length(gdat$lambda)
  smodel <- rstan::stan_model(file.path(stan_filedir, stan_file))
  rstan::rstan_options(threads_per_chain=threads_per_chain)
  if (is.null(start) || !file.exists(fpart)) {
    init_tracked(nsim, n_st, n_em, nr)
    start <- 1
  } else {
    load(fpart)
  }
  if (is.null(end)) {
    end <- nr
  }
  for (i in start:end) {
    if (is.na(dz[i]) || is.na(nnfits$tauv[i])) next
    sfit <- stanfit_one(gdat, lib.mod, nnfits, dz, which.spax=i,
                          prep_data,
                            init_opt,
                            init_sampler,
                            stan_model=smodel,
                            iter_opt=iter_opt, jv=jv,
                            iter = iter, warmup = warmup, chains = chains,
                            max_treedepth = max_treedepth, adapt_delta = adapt_delta,
                            open_progress=open_progress, ...)
    plotpp(sfit, title=paste("fiber", i))
    update_tracked(i, sfit, fpart)
    rm(sfit)
    }
    return_tracked()
}

