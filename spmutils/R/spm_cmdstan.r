cmdstanfit_one <- function(gdat, dz, nnfits, which.spax,
                           prep_data = prep_data_mod,
                           init_opt = init_opt_mod,
                           init_sampler = init_sampler_cmd,
                           stan_model=NULL,
                           stan_file="spm_dust_mod_psum.stan", stan_filedir="~/spmcode/",
                           iter_opt=5000, 
                           jv=1.e-5,
                           iter=750, warmup=250, thin=1, 
                           chains=4, cores=4, threads=4,
                           maxtree=11, adapt_delta=0.9
) {
  
  require(cmdstanr)
  if (is.null(stan_model)) {
    stan_model <- cmdstan_model(stan_file=file.path(stan_filedir, stan_file), cpp_options = list(stan_threads = TRUE))
  }
  spm_data <- prep_data(gdat, dz, nnfits, which.spax)
  inits <- init_opt(spm_data, nnfits, which.spax, jv)
  spm_opt <- stan_model$optimize(data=spm_data, init=list(inits), iter=iter_opt, threads=1)
 
  init_files <- init_sampler(stan_opt=spm_opt$mle(), jv=jv, chains=chains)
  
#  init_pars <- lapply(X=1:chains, init_sampler, stan_opt=spm_opt$par, jv=jv)
  
  
  spm_sample <- stan_model$sample(data=spm_data, init=init_files,
                                  chains=chains, parallel_chains=cores, threads_per_chain=threads,
                                  iter_warmup=warmup, iter_sampling=iter, thin=thin,
                                  max_treedepth=maxtree, adapt_delta=adapt_delta
  )
  stanfit <- rstan::read_stan_csv(spm_sample$output_files())
  
  list(spm_data=spm_data, stanfit=stanfit, 
       norm_g=spm_data$norm_g, norm_st=spm_data$norm_st, norm_em=spm_data$norm_em, in_em=spm_data$in_em)
}



cmdstanfit_batch <- function(gdat, dz, nnfits,
                        init_tracked = init_tracked_mod,
                        update_tracked = update_tracked_mod,
                        return_tracked = return_tracked_mod,
                        prep_data = prep_data_mod,
                        init_opt = init_opt_mod,
                        init_sampler = init_sampler_cmd,
                        stan_file="spm_dust_mod_psum.stan", stan_filedir="~/spmcode/",
                           iter_opt=5000, 
                           jv=1.e-5,
                           iter=750, warmup=250, thin=1, 
                           chains=4, cores=4, threads=4,
                           maxtree=11, adapt_delta=0.9,
                           start=NULL, end=NULL, fpart="bfits.rda"
) {
  
    dims <- dim(gdat$flux)
    dz <- dz$dz
    nsim <- iter*chains
    nt <- length(ages)
    nr <- dims[1]
    n_st <- ncol(lib.ssp)-1
    n_em <- length(emlines)
    nl <- length(gdat$lambda)
    smodel <- cmdstanr::cmdstan_model(stan_file=file.path(stan_filedir, stan_file), cpp_options = list(stan_threads = TRUE))
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
        cat(paste("fiber", i, "\n"))
        if (is.na(dz[i]) || is.na(nnfits$tauv[i])) next
        sfit <- cmdstanfit_one(gdat, dz, nnfits, which.spax=i,
                            prep_data = prep_data,
                            init_opt = init_opt,
                            init_sampler = init_sampler,
                            stan_model=smodel,
                            iter_opt=iter_opt, jv=jv,
                            iter = iter, warmup = warmup, thin=thin, 
                            chains = chains, cores=cores, threads=threads,
                            maxtree=maxtree, adapt_delta=adapt_delta
                            )
        plotpp(sfit, title=paste("fiber", i))
        update_tracked(i, sfit, fpart)
        rm(sfit)
    }
    return_tracked()
}

init_opt_cmd <- function(spm_data, nnfits, which.spax, jv) {
  init_vals <- init_opt_mod(spm_data, nnfits, which.spax, jv)
  fname <- tempfile(fileext = ".json")
  cmdstanr::write_stan_json(init_vals, fname)
  fname
}

init_sampler_cmd <- function(stan_opt, jv, chains) {
  a_ind <- grep("\\ba\\b", names(stan_opt))
  b_st_s_ind <- grep("b_st_s", names(stan_opt))
  b_em_ind <- grep("b_em", names(stan_opt))
  tauv_ind <- grep("tauv", names(stan_opt), fixed=TRUE)
  delta_ind <- grep("delta", names(stan_opt), fixed=TRUE)
  files <- character(chains)
  for (i in 1:chains) {
    files[i] <- tempfile(fileext = ".json")
    a <- as.numeric(stan_opt[a_ind] + rnorm(1, sd=jv))
    b_st_s <- stan_opt[b_st_s_ind] + runif(length(b_st_s_ind), min=jv/10, max=jv)
    b_st_s <- as.vector(b_st_s/sum(b_st_s))
    b_em <- as.vector(stan_opt[b_em_ind] + runif(length(b_em_ind), min=jv/10, max=jv))
    tauv <- as.numeric(stan_opt[tauv_ind] + runif(1, min=jv/10, max=jv))
    delta <- as.numeric(stan_opt[delta_ind] + rnorm(1, sd=jv))
    init_vals <- list(a=a, b_st_s=b_st_s, b_em=b_em, tauv=tauv, delta=delta)
    cmdstanr::write_stan_json(init_vals, files[i])
  }
  files
}
    
