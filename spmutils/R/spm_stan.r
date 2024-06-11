stanfit_one <- function(gdat, dz, nnfits, which.spax,
                        prep_data = prep_data_mod,
                        init_opt = init_opt_mod,
                        init_sampler = init_sampler_mod,
                        stan_model=NULL,
                        stan_file="spm_dust_mod_simpl.stan", stan_filedir="~/spmcode/",
                        iter_opt=5000, 
                        jv=1.e-4,
                        iter=1000, warmup=250, thin=1, chains=4, 
                        open_progress=TRUE, ...) {
    
    require(rstan)
    if (is.null(stan_model)) {
      stan_model <- rstan::stan_model(file=file.path(stan_filedir, stan_file))
    }
    
    spm_data <- prep_data(gdat, dz, nnfits, which.spax)
    inits <- init_opt(spm_data, nnfits, which.spax, jv)
    spm_opt <- optimizing(stan_model, data=spm_data, init=inits, as_vector=FALSE, verbose=TRUE, iter=iter_opt)
    
    init_pars <- lapply(X=1:chains, init_sampler, stan_opt=spm_opt$par, jv=jv)
    
    stanfit <- sampling(stan_model, data=spm_data,
                     chains=chains, iter=iter, warmup=warmup, thin=thin,
                     cores=getOption("mc.cores"),
                     init=init_pars, open_progress=open_progress, ...)
    
    list(spm_data=spm_data, stanfit=stanfit, 
         norm_g=spm_data$norm_g, norm_st=spm_data$norm_st, norm_em=spm_data$norm_em, in_em=spm_data$in_em)
}


                        
stanfit_batch <- function(gdat, dz, nnfits,
                        init_tracked = init_tracked_mod,
                        update_tracked = update_tracked_mod,
                        return_tracked = return_tracked_mod,
                        prep_data = prep_data_mod,
                        init_opt = init_opt_mod,
                        init_sampler = init_sampler_mod,
                        stan_file="spm_dust_mod_psum.stan", stan_filedir="~/spmcode/",
                        iter_opt=5000, 
                        jv=1.e-4,
                        iter=1000, warmup=250, chains=4,
                        open_progress=TRUE,
                        start=NULL, end=NULL, fpart="bfits.rda", ...) {
    dims <- dim(gdat$flux)
    nsim <- (iter-warmup)*chains
    nt <- length(ages)
    nr <- dims[1]
    n_st <- ncol(lib.ssp)-1
    n_em <- length(emlines)
    nl <- length(gdat$lambda)
    smodel <- rstan::stan_model(file.path(stan_filedir, stan_file))
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
        sfit <- stanfit_one(gdat, dz, nnfits, which.spax=i,
                            prep_data,
                            init_opt,
                            init_sampler,
                            stan_model=smodel,
                            iter_opt=iter_opt, jv=jv,
                            iter = iter, warmup = warmup, chains = chains, 
                            open_progress=open_progress, ...)
        plot(plotpp(sfit)+ggtitle(paste("fiber =", i)))
        update_tracked(i, sfit, fpart)
        rm(sfit)
    }
    return_tracked()
}

