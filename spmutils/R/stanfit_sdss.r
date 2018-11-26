stanfit_sdss <- function(gdat, nnfits, 
                        lib.ssp, age, Z,
                        emlines=spmutils::lambda_em,
                        stan_model=NULL,
                        stan_file="spm_dust_norm.stan", stan_filedir="~/spmcode/",
                        iter=500, warmup=250, thin=1, chains=4, 
                        OP=FALSE, ...) {
    
    require(rstan)
    z <- gdat$meta$z
    flux <- gdat$flux
    ivar <- gdat$ivar
    lambda.rest <- gdat$lambda
    logl <- log10(lambda.rest)
    nsim <- (iter-warmup)*chains
    T.gyr <- 10^(age-9)
        
    lib.ssp$lambda <- spmutils::airtovac(lib.ssp$lambda)
    lib.st <- regrid(lambda.rest, lib.ssp)
    x.st <- blur.lib(lib.st, nnfits$vdisp.st)
    n.st <- ncol(x.st)
    nz <- length(Z)
    nt <- n.st/nz
    ind.young <- (0:(nz-1))*nt+1
    allok <- complete.cases(flux, ivar, x.st)
    nl <- length(which(allok))
    
    x.em <- make_emlib(emlines, nnfits$vdisp.em, logl, allok)
    in.em <- x.em$in_em
    x.em <- x.em$x_em
    n.em <- ncol(x.em)
    
    tauv <- nnfits$tauv
    
    b.st <- nnfits$nnfits[1:n.st]
    b.em <- nnfits$nnfits[(n.st+1):(n.st+n.em)]
    norm.st <- max(b.st)
    norm.em <- max(b.em)
    x.st <- x.st*norm.st
    x.em <- x.em*norm.em
    dT <- rep(diff(c(0,T.gyr)), nz)
    
    spm_data <- list(nt=nt, nz=nz, nl=nl, n_em=n.em,
                     ind_young=ind.young,
                     lambda=lambda.rest[allok],
                     gflux=flux[allok],
                     g_std=1/sqrt(ivar[allok]),
                     sp_st=x.st[allok,],
                     dT=dT,
                     sp_em=x.em[allok,])
    deps <- sqrt(.Machine$double.eps)
    b_st <- b.st/norm.st
    b_em <- b.em/norm.em
    b_em[b_em==0] <- deps
    if (tauv == 0) tauv <- deps
    
    init_pars <- function(chain_id) {
        b_st0 <- b_st
        atzero <- which(b_st0 == 0)
        b_st0[atzero] <- runif(length(atzero), 0, .001)
        list(b_st=b_st0, b_em=b_em, tauv=tauv)
    }
    
    if (!is.null(stan_model)) {
      stanfit <- sampling(stan_model, data=spm_data,
                     chains=chains, iter=iter, warmup=warmup, thin=thin,
                     cores=min(chains, getOption("mc.cores")),
                     init=init_pars, open_progress=OP, ...)
    
    } else {
      stanfit <- stan(file=file.path(stan_filedir, stan_file), data=spm_data,
                     chains=chains, iter=iter, warmup=warmup, thin=thin,
                     cores=min(chains, getOption("mc.cores")),
                     init=init_pars, open_progress=OP, ...)
    }
    list(spm_data=spm_data, stanfit=stanfit, norm.st=norm.st, norm.em=norm.em, in.em=in.em)
}
