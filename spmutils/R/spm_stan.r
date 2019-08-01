stanfit_one <- function(gdat, dz, nnfits, which.spax, 
                        stan_model=NULL,
                        stan_file="spm_dust_simpl.stan", stan_filedir="~/spmcode/",
                        iter_opt=5000, 
                        init_opt = list(npars=4, is_simpl=c(FALSE, TRUE, FALSE, FALSE), jitter=TRUE, jv=1.e-5),
                        iter=500, warmup=250, thin=1, chains=4, 
                        OP=FALSE, ...) {
    
    require(rstan)
    i <- which.spax
    z <- gdat$meta$z+dz[i]
    flux <- gdat$flux[i, ]
    ivar <- gdat$ivar[i, ]
    lambda.rest <- gdat$lambda/(1+z)
    logl <- log10(lambda.rest)
    nsim <- (iter-warmup)*chains
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
    norm_g <- mean(flux[allok])
    norm_st <- colMeans(x.st[allok,])
    x.st <- scale(x.st[allok,], center=FALSE, scale=norm_st)
    
    x.em <- make_emlib(emlines, nnfits$vdisp.em[i], logl, allok)
    in.em <- x.em$in_em
    x.em <- x.em$x_em
    n.em <- ncol(x.em)
    norm_em <- 1
    
    spm_data <- list(nt=nt, nz=nz, nl=nl, n_em=n.em,
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
    
    if (is.null(stan_model)) {
      stan_model <- rstan::stan_model(file=file.path(stan_filedir, stan_file))
    }
    
    spm_opt <- optimizing(stan_model, data=spm_data, as_vector=FALSE, verbose=TRUE, iter=iter_opt)
    names_pars <- names(spm_opt$par)[1:init_opt$npars]
    
    init_pars <- function(chain_id) {
      retval <- list()
      for (i in 1:init_opt$npars) {
        assign(names_pars[i], spm_opt$par[[names_pars[i]]])
        if (init_opt$jitter) {
          assign(names_pars[i], get(names_pars[i]) +
                  runif(length(get(names_pars[i])), max= init_opt$jv))
        }
        if (init_opt$is_simpl[i]) {
          assign(names_pars[i], get(names_pars[i])/sum(get(names_pars[i])))
        }
        retval[[names_pars[i]]] <- get(names_pars[i])
      }
      retval
    }
    
    stanfit <- sampling(stan_model, data=spm_data,
                     chains=chains, iter=iter, warmup=warmup, thin=thin,
                     cores=min(chains, getOption("mc.cores")),
                     init=init_pars, open_progress=OP, ...)
    
    list(spm_data=spm_data, stanfit=stanfit, norm_g=norm_g, norm_st=norm_st, norm_em=norm_em, in.em=in.em)
}


                        
stanfit_batch <- function(gdat, dz, nnfits,
                        stan_file="spm_dust_norm.stan", stan_filedir="~/spmcode/",
                        iter_opt=5000, 
                        init_opt = list(npars=4, is_simpl=c(FALSE, TRUE, FALSE, FALSE), jitter=TRUE, jv=1.e-5),
                        iter=500, warmup=250, chains=4,
                        OP=FALSE,
                        start=NULL, end=NULL, fpart="bfits.rda", ...) {
    dims <- dim(gdat$flux)
    dz <- dz$dz
    nsim <- (iter-warmup)*chains
    nt <- length(ages)
    nr <- dims[1]
    n_st <- ncol(lib.ssp)-1
    n_em <- length(emlines)
    nl <- length(gdat$lambda)
    smodel <- rstan::stan_model(file.path(stan_filedir, stan_file))
    if (is.null(start) || !file.exists(fpart)) {
      b_st <- array(NA, dim=c(nsim, n_st, nr))
      b_em <- array(NA, dim=c(nsim, n_em, nr))
      a <- matrix(NA, nrow=nsim, ncol=nr)
      tauv <- matrix(NA, nrow=nsim, ncol=nr)
      in_em <- matrix(NA, nrow=n_em, ncol=nr)
      walltime <- rep(NA, nr)
      divergences <- rep(NA, nr)
      max_treedepth <- rep(NA, nr)
      norm_g <- rep(NA, nr)
      norm_st <- matrix(NA, nrow=n_st, ncol=nr)
      norm_em <- rep(NA, nr)
      start <- 1
    } else {
        load(fpart)
    }
    if (is.null(end)) {
      end <- nr
    }
    for (i in start:end) {
        if (is.na(dz[i]) || is.na(nnfits$Mstar[i])) next
        sfit <- stanfit_one(gdat, dz, nnfits, which.spax=i,
                            stan_model=smodel,
                            iter_opt=iter_opt, init_opt=init_opt,
                            iter = iter, warmup = warmup, chains = chains, 
                            OP=OP, ...)
        plot(plotpp(sfit)+ggtitle(paste("fiber =", i)))
        post <- extract(sfit$stanfit)
        b_st[,,i] <- post$b_st
        b_em[,sfit$in.em,i] <- post$b_em
        a[,i] <- post$a
        tauv[,i] <- post$tauv
        in_em[sfit$in.em,i] <- sfit$in.em
        walltime[i] <- max(rowSums(get_elapsed_time(sfit$stanfit)))
        sp <- get_sampler_params(sfit$stanfit, inc_warmup=FALSE)
        divergences[i] <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
        max_treedepth[i] <- max(sapply(sp, function(x) max(x[, "treedepth__"])))
        norm_g[i] <- sfit$norm_g
        norm_st[,i] <- sfit$norm_st
        norm_em[i] <- sfit$norm_em
        save(b_st, b_em, tauv, a, in_em,
                  walltime, divergences, max_treedepth,
                  norm_g, norm_st, norm_em, file = fpart)
        rm(sfit, post, sp)
    }
    list(b_st=b_st, b_em=b_em, a=a, tauv=tauv, in_em=in_em,
                  walltime=walltime, divergences=divergences, max_treedepth=max_treedepth,
                  norm_g=norm_g, norm_st=norm_st, norm_em=norm_em)
}

## star formation history, mass growth history, etc.

get_sfh <- function(..., z, fibersinbin=1, density=TRUE, tsf=0.1) {
  ins <- list(...)
  if (is.list(ins[[1]])) {
    ins <- ins[[1]]
    post <- rstan::extract(ins$stanfit)
    b_st <- post$b_st
    norm_st <- ins$norm_st
    if (exists("norm_g", ins)) {
      b_st <- b_st*ins$norm_g
      norm_st <- 1/norm_st
    }
  } else {
    b_st <- ins$b_st
    norm_st <- ins$norm_st
  }
  nsim <- nrow(b_st)
  nt <- length(ages)
  nz <- ncol(b_st)/nt
  T.gyr <- 10^(ages-9)
  isf <- which.min(abs(tsf-T.gyr))
  if (density) {
    binarea <- log10(pi*fibersinbin*cosmo::ascale(z)^2)
  } else {
    binarea <- 0
  }
  b_st <- t(t(b_st)*norm_st)*cosmo::lum.sol(1, z)
  rmass <- t(t(b_st) * mstar)
  sfh_post <- matrix(0, nsim, nt)
  mgh_post <- matrix(0, nsim, nt)
  for (i in 1:nz) {
    sfh_post <- sfh_post + b_st[,((i-1)*nt+1):(i*nt)]
    mgh_post <- mgh_post + rmass[,((i-1)*nt+1):(i*nt)]
  }
  totalmg_post <- cbind(rowSums(mgh_post), rowSums(mgh_post) - t(apply(mgh_post, 1, cumsum)))
  mgh_post <- cbind(rep(1, nsim), 1 - (t(apply(mgh_post, 1, cumsum))/rowSums(mgh_post)))
  sfr <- log10(rowSums(sfh_post[,1:isf])/T.gyr[isf])-9.
  relsfr <- sfr - log10(rowSums(sfh_post)/(cosmo::dcos(Inf)$dT-cosmo::dcos(z)$dT)) + 6
  sigma_sfr <- sfr - binarea
  sigma_mstar <- log10(b_st %*% mstar) - binarea
  ssfr <- sigma_sfr - sigma_mstar
  list(sfh_post=sfh_post, mgh_post=mgh_post, totalmg_post=totalmg_post,
       sigma_mstar=sigma_mstar, sigma_sfr=sigma_sfr, ssfr=ssfr, relsfr=relsfr)
}

batch_sfh <- function(gdat, sfits, tsf=0.1) {
  b_st <- sfits$b_st
  norm_st <- sfits$norm_st
  dims <- dim(b_st)
  nsim <- dims[1]
  nt <- length(ages)
  nz <- dims[2]/nt
  nf <- dims[3]
  if (exists("norm_g", sfits)) {
    norm_g <- sfits$norm_g
    norm_st <- 1/norm_st
  } else {
    norm_g <- rep(1, nf)
  }
  z <- gdat$meta$z
  if (!exists("fibersinbin", gdat)) {
    fibersinbin <- rep(1, nf)
  } else {
    fibersinbin <- gdat$fibersinbin
  }
  
  sfh_post <- array(NA, dim=c(nsim, nt, nf))
  mgh_post <- array(0, dim=c(nsim, nt+1, nf))
  totalmg_post <- matrix(0, nsim, nt+1)
  sigma_mstar <- matrix(NA, nsim, nf)
  sigma_sfr <- matrix(NA, nsim, nf)
  ssfr <- matrix(NA, nsim, nf)
  relsfr <- matrix(NA, nsim, nf)
  
  for (i in 1:nf) {
    if (is.na(b_st[1, 1, i])) next
    sfi <- get_sfh(b_st=b_st[,,i]*norm_g[i], norm_st=norm_st[,i], z=z, fibersinbin=fibersinbin[i])
    sfh_post[,,i] <- sfi$sfh_post
    mgh_post[,,i] <- sfi$mgh_post
    totalmg_post <- totalmg_post + sfi$totalmg_post
    sigma_mstar[, i] <- sfi$sigma_mstar
    sigma_sfr[, i] <- sfi$sigma_sfr
    ssfr[, i] <- sfi$ssfr
    relsfr[,i] <- sfi$relsfr
  }
  list(sfh_post=sfh_post, mgh_post=mgh_post, totalmg_post=totalmg_post,
       sigma_mstar=sigma_mstar, sigma_sfr=sigma_sfr, ssfr=ssfr, relsfr=relsfr)
}
  
  
## some sorta useful summary measures

get_proxies <- function(...) {
  ins <- list(...)
  if (is.list(ins[[1]])) {
    ins <- ins[[1]]
    post <- rstan::extract(ins$stanfit)
    b_st <- post$b_st
    norm_st <- ins$norm_st
    if (exists("a", post)) {
      b_st <- b_st*post$a*ins$norm_g
      norm_st <- 1/norm_st
    }
  } else {
    b_st <- ins$b_st
    norm_st <- ins$norm_st
  }
  b_st <- t(t(b_st)*norm_st)
  nz <- length(Z)
  nt <- length(ages)
  T.gyr <- 10^(ages-9)
  tbar <- log10((b_st %*% rep(T.gyr, nz))/rowSums(b_st)) + 9
  tbar_lum <- log10((b_st %*% (rep(T.gyr, nz) * gri.ssp["r",]))/
                   ((b_st %*% gri.ssp["r",]))) + 9
  g_i <- 2.5*log10(b_st %*% gri.ssp["i",]) - 2.5*log10(b_st %*% gri.ssp["g",])
  data.frame(tbar=tbar, tbar_lum=tbar_lum, g_i=g_i)
}

## emission line fluxes, luminosity density, equivalent width

get_em <- function(..., z, fibersinbin=1, ew_width=15) {
  ins <- list(...)
  if (is.list(ins[[1]])) {
    ins <- ins[[1]]
    post <- rstan::extract(ins$stanfit)
    b_st <- post$b_st
    b_em <- post$b_em
    norm_st <- ins$norm_st
    norm_em <- ins$norm_em
    if (exists("norm_g", ins)) {
      b_st <- b_st*ins$norm_g
      b_em <- b_em*ins$norm_g
      norm_st <- 1/norm_st
    }
 } else {
    b_st <- ins$b_st
    b_em <- ins$b_em
    norm_st <- ins$norm_st
    norm_em <- ins$norm_em
  }
  emlines <- emlines[ins$in.em]
  b_st <- t(t(b_st)*norm_st)
  ne <- ncol(b_em)
  nsim <- nrow(b_em)
  binarea <- log10(pi*fibersinbin*cosmo::ascale(z)^2)
  em.mult <- emlines*log(10)/10000

  flux_em <- matrix(NA, nsim, ne)
  sigma_logl_em <- matrix(NA, nsim, ne)
  ew_em <- matrix(NA, nsim, ne)
  
  flux_em <- t(t(b_em)*em.mult)*norm_em
  flux_em[flux_em < 0] <- 0
  sigma_logl_em <- cosmo::loglum.ergs(flux_em, z) - binarea
  
  mu_st <- tcrossprod(b_st, as.matrix(lib.ssp[, -1]))
  il_em <- findInterval(emlines, lib.ssp$lambda)
  for (i in 1:ne) {
    intvl <- (il_em[i]-ew_width):(il_em[i]+ew_width)
    fc <- rowMeans(mu_st[, intvl])
    ew_em[, i] <- flux_em[, i]/fc
  }
  colnames(flux_em) <- names(emlines)
  colnames(sigma_logl_em) <- names(emlines)
  colnames(ew_em) <- names(emlines)
  list(flux_em=flux_em, sigma_logl_em=sigma_logl_em, ew_em=ew_em)
}

batch_em <- function(gdat, sfits, ew_width=15) {
  nsim <- dim(sfits$b_em)[1]
  ne <- length(emlines)
  nf <- length(gdat$xpos)
  norm_st <- sfits$norm_st
  if (exists("norm_g", sfits)) {
    norm_g <- sfits$norm_g
    norm_st <- 1/sfits$norm_st
  } else {
    norm_g <- rep(1, nf)
  }
  
  flux_em <- array(NA, dim=c(nsim, ne, nf))
  sigma_logl_em <- array(NA, dim=c(nsim, ne, nf))
  ew_em <- array(NA, dim=c(nsim, ne, nf))
  
  for (i in 1:nf) {
    if (is.na(sfits$b_st[1, 1, i])) next
    in.em <- sfits$in_em[!is.na(sfits$in_em[,i]), i]
    if(is.null(dim(norm_st))) {
      nst <- norm_st[i]
    } else {
      nst <- norm_st[,i]
    }
    emi <- get_em(b_em=sfits$b_em[,in.em,i]*norm_g[i], b_st=sfits$b_st[,,i]*norm_g[i], 
                  norm_em=sfits$norm_em[i], norm_st=nst, 
                  in.em=in.em, z=gdat$meta$z,
                  fibersinbin=gdat$fibersinbin[i], ew_width=ew_width)
    flux_em[,in.em,i] <- emi$flux_em
    sigma_logl_em[,in.em,i] <- emi$sigma_logl_em
    ew_em[,in.em,i] <- emi$ew_em
  }
  dimnames(flux_em)[[2]] <- names(emlines)
  dimnames(sigma_logl_em)[[2]] <- names(emlines)
  dimnames(ew_em)[[2]] <- names(emlines)
  list(flux_em=flux_em, sigma_logl_em=sigma_logl_em, ew_em=ew_em)
}

## bpt class from [N II]/Halpha

batch_bptclass <- function(flux_em, snthresh=3) {
  nb <- dim(flux_em)[3]
  bpt <- as.factor(rep("NO EM", nb))
  levels(bpt) <- c("NO EM", "EL", "SF", "COMP", "LINER", "AGN")
  f_m <- apply(flux_em, c(2, 3), mean)
  f_sd <- apply(flux_em, c(2, 3), sd)
  for (i in 1:nb) {
    if (all(is.na(f_m[, i]))) {
      bpt[i] <- NA
      next
    }
    if(any(f_m[, i]/f_sd[, i] > snthresh, na.rm=TRUE)) bpt[i] <- "EL"
    if (any(is.na(f_m[c("h_beta", "oiii_5007", "h_alpha", "nii_6584"), i]))) next
    if (f_m["h_beta", i]/f_sd["h_beta", i] > snthresh &&
        f_m["oiii_5007", i]/f_sd["oiii_5007", i] > snthresh &&
        f_m["h_alpha", i]/f_sd["h_alpha", i] > snthresh &&
        f_m["nii_6584", i]/f_sd["nii_6584", i] > snthresh) {
      o3hbeta <- log10(f_m["oiii_5007", i]/f_m["h_beta", i])
      n2halpha <- log10(f_m["nii_6584", i]/f_m["h_alpha", i])
      if ((o3hbeta <= 0.61/(n2halpha-0.05)+1.3) &&
         (n2halpha <= 0.05)) {
        bpt[i] <- "SF"
        next
      }
      if ((o3hbeta > 0.61/(n2halpha-0.05)+1.3 || n2halpha > 0.05) &&
          (o3hbeta <= 0.61/(n2halpha-0.47)+1.19)) {
        bpt[i] <- "COMP"
        next
      }
      if ((o3hbeta > 0.61/(n2halpha-0.47)+1.19 || n2halpha > 0.47) &&
          (o3hbeta > 1.05*n2halpha+0.45)) {
        bpt[i] <- "AGN"
      } else {
        bpt[i] <- "LINER"
      }
    }
  }
  bpt
}


## emission line ratios and various "strong line" metallicity calibrations
  
get_lineratios <- function(flux_em, tauv, tauv_mult=1, alaw=calzetti) {
  o3hbeta <- log10(flux_em[,"oiii_5007"]/flux_em[,"h_beta"])
  o1halpha <- log10(flux_em[,"oi_6300"]/flux_em[,"h_alpha"])
  n2halpha <- log10(flux_em[,"nii_6584"]/flux_em[,"h_alpha"])
  s2halpha <- log10((flux_em[,"sii_6717"]+flux_em[,"sii_6730"])/flux_em[,"h_alpha"])
  
  o2 <- (flux_em[,"oii_3727"]+flux_em[,"oii_3729"])*alaw(3728., -tauv*tauv_mult)
  o3 <- (flux_em[,"oiii_4959"]+flux_em[,"oiii_5007"])*alaw(4980., -tauv*tauv_mult)
  hb <- flux_em[,"h_beta"]*alaw(4863., -tauv*tauv_mult)
  
  r23 <- log10((o2+o3)/hb)
  o3n2 <- o3hbeta-n2halpha

  ## log(O/H) estimates from Pettini & Pagel 2004 or Tremonti et al. 2004
    
  oh_n2 <- 9.37+2.03*n2halpha+1.26*n2halpha^2+0.32*n2halpha^3
  oh_o3n2 <- 8.73-0.32*o3n2
  oh_r23 <- 9.185-0.313*r23-0.264*r23^2-0.321*r23^3
  
  data.frame(o3hbeta=o3hbeta, o1halpha=o1halpha, n2halpha=n2halpha, s2halpha=s2halpha,
             r23=r23, o3n2=o3n2, oh_n2=oh_n2, oh_o3n2=oh_o3n2, oh_r23=oh_r23)
}

  
sum_batchfits <- function(gdat, nnfits, sfits, drpcat, alaw=calzetti, intr_bd=2.86, clim=0.95) {
    
    tauv.bd <- function(flux_em, intr_bd, alaw) {
      bd <- flux_em[,'h_alpha']/flux_em[,'h_beta']
      bd[!is.finite(bd)] <- NA
      tauv <- log(bd/intr_bd)/(log(alaw(6562.8,1))-log(alaw(4861.3,1)))
      tauv[tauv<0] <- 0
      tauv
    }
    
    logl.ha.cor <- function(logl.halpha, tauv, alaw) {
      att <- alaw(lambda=6562.8, tauv)
      logl.halpha - log10(att)
    }
    
    nf <- length(gdat$xpos)
    fibersinbin <- rep(1, nf)
    if (exists("bin.fiber", gdat)) {
      if (!exists("fibersinbin", gdat)) {
        bin.fiber <- gdat$bin.fiber
        bin.no <- unique(bin.fiber[!is.na(bin.fiber)])
        for (i in seq_along(bin.no)) {
          fibersinbin[i] <- length(which(bin.fiber == bin.no[i]))
        }
      } else {
        fibersinbin <- gdat$fibersinbin
      }
    }
    
    nsim <- nrow(sfits$b_st)
    nt <- length(ages)
    fiberarea <- pi*cosmo::ascale(gdat$meta$z)^2
    plateifu <- rep(gdat$meta$plateifu, nf)
    
    ## projected distance in kpc and relative to effective radius
    
    d_kpc <- cosmo::ascale(gdat$meta$z)*sqrt(gdat$xpos^2+gdat$ypos^2)
    d_re <- sqrt(gdat$xpos^2+gdat$ypos^2)/
    drpcat$nsa_petro_th50[match(gdat$meta$plateifu,drpcat$plateifu)]
    
    
    ## stuff taken from nnfits
    
    d4000_n <- nnfits$d4000_n
    d4000_n_err <- nnfits$d4000_n_err
    lick_hd_a <- nnfits$lick[,'HdeltaA']
    lick_hd_a_err <- nnfits$lick.err[,'HdeltaA_err']
    mgfe <- sqrt(nnfits$lick[,'Mg_b']*(0.72*nnfits$lick[,'Fe5270']+0.28*nnfits$lick[,'Fe5335']))
    mgfe[is.nan(mgfe)] <- NA
    
    ## tauv from batch fits
    
    tauv_m <- colMeans(sfits$tauv)
    tauv_std <- apply(sfits$tauv, 2, sd)
    
    mgh_post <- array(NA, dim=c(nsim, nt+1, nf))
    sfh_post <- array(NA, dim=c(nsim, nt, nf))
    totalmg_post <- matrix(0, nsim, nt+1)
    
    varnames <- c("sigma_mstar", "sigma_sfr", "ssfr", "relsfr",
                  "tbar", "tbar_lum", "g_i",
                  "tauv_bd", "sigma_logl_ha", "sigma_logl_ha_ctauv", "sigma_logl_ha_ctauv_bd",
                  "eqw_ha", "o3hbeta", "o1halpha", "n2halpha", "s2halpha",
                  "r23", "o3n2", "oh_n2", "oh_o3n2", "oh_r23")
    suffixes <- c("m", "std", "lo", "hi")
    
    for (i in seq_along(varnames)) {
      for (j in seq_along(suffixes)) {
        assign(paste(varnames[i], suffixes[j], sep="_"), numeric(nf))
      }
    }
    
    sfh_all <- batch_sfh(gdat, sfits)
    
    for (i in 1:4) {
      assign(paste(varnames[i], suffixes[1], sep="_"), colMeans(sfh_all[[varnames[i]]]))
      assign(paste(varnames[i], suffixes[2], sep="_"), apply(sfh_all[[varnames[i]]], 2, sd))
    }
      
        
    em_all <- batch_em(gdat, sfits)
    bpt <- batch_bptclass(em_all$flux_em)
    
    
    for (i in 1:nf) {
      
      ## star formation history, etc.
      
      quants <- hdiofmcmc(sfh_all$sigma_mstar[,i])
      sigma_mstar_lo[i] <- quants[1]
      sigma_mstar_hi[i] <- quants[2]
      
      quants <- hdiofmcmc(sfh_all$sigma_sfr[,i])
      sigma_sfr_lo[i] <- quants[1]
      sigma_sfr_hi[i] <- quants[2]
      
      quants <- hdiofmcmc(sfh_all$ssfr[,i])
      ssfr_lo[i] <- quants[1]
      ssfr_hi[i] <- quants[2]
      
      quants <- hdiofmcmc(sfh_all$relsfr[,i])
      relsfr_lo[i] <- quants[1]
      relsfr_hi[i] <- quants[2]
      
      if (is.null(dim(sfits$norm_st))) {
        norm_st <- sfits$norm_st[i]
      } else {
        norm_st <- 1/sfits$norm_st[,i]
      }
      proxi <- get_proxies(b_st=sfits$b_st[,,i], norm_st=norm_st)
      
      tbar_m[i] <- mean(proxi$tbar)
      tbar_std[i] <- sd(proxi$tbar)
      quants <- hdiofmcmc(proxi$tbar, credmass=clim)
      tbar_lo[i] <- quants[1]
      tbar_hi[i] <- quants[2]
      
      tbar_lum_m[i] <- mean(proxi$tbar_lum)
      tbar_lum_std[i] <- sd(proxi$tbar_lum)
      quants <- hdiofmcmc(proxi$tbar_lum, credmass=clim)
      tbar_lum_lo[i] <- quants[1]
      tbar_lum_hi[i] <- quants[2]
      
      g_i_m[i] <- mean(proxi$g_i)
      g_i_std[i] <- sd(proxi$g_i)
      quants <- hdiofmcmc(proxi$g_i, credmass=clim)
      g_i_lo[i] <- quants[1]
      g_i_hi[i] <- quants[2]
      
      linesi <- get_lineratios(em_all$flux_em[,,i], sfits$tauv[,i], alaw=alaw)
      
      tauv_bd <- tauv.bd(em_all$flux_em[,,i], intr_bd=intr_bd, alaw=alaw)
      
      tauv_bd_m[i] <- mean(tauv_bd)
      tauv_bd_std[i] <- sd(tauv_bd)
      quants <- hdiofmcmc(tauv_bd, credmass=clim)
      tauv_bd_lo[i] <- quants[1]
      tauv_bd_hi[i] <- quants[2]
      
      ## uncorrected Halpha luminosity from stan fits
      
      sigma_logl_ha_m[i] <- mean(em_all$sigma_logl_em[,"h_alpha", i])
      sigma_logl_ha_std[i] <- sd(em_all$sigma_logl_em[,"h_alpha", i])
      quants <- hdiofmcmc(em_all$sigma_logl_em[,"h_alpha", i], credmass=clim)
      sigma_logl_ha_lo[i] <- quants[1]
      sigma_logl_ha_hi[i] <- quants[2]
      
      ## correct from stan fit estimate of tauv
      
      logl_ha_c <- logl.ha.cor(em_all$sigma_logl_em[, "h_alpha", i], sfits$tauv[,i], alaw=alaw)
      
      sigma_logl_ha_ctauv_m[i] <- mean(logl_ha_c)
      sigma_logl_ha_ctauv_std[i] <- sd(logl_ha_c)
      quants <- hdiofmcmc(logl_ha_c, credmass=clim)
      sigma_logl_ha_ctauv_lo[i] <- quants[1]
      sigma_logl_ha_ctauv_hi[i] <- quants[2]
      
      
      ## correct from balmer decrement
      
      logl_ha_c <- logl.ha.cor(em_all$sigma_logl_em[, "h_alpha", i], tauv_bd, alaw=alaw)
      
      sigma_logl_ha_ctauv_bd_m[i] <- mean(logl_ha_c)
      sigma_logl_ha_ctauv_bd_std[i] <- sd(logl_ha_c)
      quants <- hdiofmcmc(logl_ha_c, credmass=clim)
      sigma_logl_ha_ctauv_bd_lo[i] <- quants[1]
      sigma_logl_ha_ctauv_bd_hi[i] <- quants[2]
      
      eqw_ha_m[i] <- mean(em_all$ew_em[,"h_alpha", i])
      eqw_ha_std[i] <- sd(em_all$ew_em[,"h_alpha", i])
      quants <- hdiofmcmc(em_all$ew_em[,"h_alpha", i], credmass=clim)
      eqw_ha_lo[i] <- quants[1]
      eqw_ha_hi[i] <- quants[2]
      
      ## some emission line ratios
      
      o3hbeta_m[i] <- mean(linesi$o3hbeta)
      o3hbeta_std[i] <- sd(linesi$o3hbeta)
      quants <- hdiofmcmc(linesi$o3hbeta, credmass=clim)
      o3hbeta_lo[i] <- quants[1]
      o3hbeta_hi[i] <- quants[2]
      
      o1halpha_m[i] <- mean(linesi$o1halpha)
      o1halpha_std[i] <- sd(linesi$o1halpha)
      quants <- hdiofmcmc(linesi$o1halpha, credmass=clim)
      o1halpha_lo[i] <- quants[1]
      o1halpha_hi[i] <- quants[2]
      
      n2halpha_m[i] <- mean(linesi$n2halpha)
      n2halpha_std[i] <- sd(linesi$n2halpha)
      quants <- hdiofmcmc(linesi$n2halpha, credmass=clim)
      n2halpha_lo[i] <- quants[1]
      n2halpha_hi[i] <- quants[2]
      
      s2halpha_m[i] <- mean(linesi$s2halpha)
      s2halpha_std[i] <- sd(linesi$s2halpha)
      quants <- hdiofmcmc(linesi$s2halpha, credmass=clim)
      s2halpha_lo[i] <- quants[1]
      s2halpha_hi[i] <- quants[2]
      
      r23_m[i] <- mean(linesi$r23)
      r23_std[i] <- sd(linesi$r23)
      quants <- hdiofmcmc(linesi$r23, credmass=clim)
      r23_lo[i] <- quants[1]
      r23_hi[i] <- quants[2]
      
      o3n2_m[i] <- mean(linesi$o3n2)
      o3n2_std[i] <- sd(linesi$o3n2)
      quants <- hdiofmcmc(linesi$o3n2, credmass=clim)
      o3n2_lo[i] <- quants[1]
      o3n2_hi[i] <- quants[2]

      oh_n2_m[i] <- mean(linesi$oh_n2)
      oh_n2_std[i] <- sd(linesi$oh_n2)
      quants <- hdiofmcmc(linesi$oh_n2, credmass=clim)
      oh_n2_lo[i] <- quants[1]
      oh_n2_hi[i] <- quants[2]
      
      oh_o3n2_m[i] <- mean(linesi$oh_o3n2)
      oh_o3n2_std[i] <- sd(linesi$oh_o3n2)
      quants <- hdiofmcmc(linesi$oh_o3n2, credmass=clim)
      oh_o3n2_lo[i] <- quants[1]
      oh_o3n2_hi[i] <- quants[2]
      
      oh_r23_m[i] <- mean(linesi$oh_r23)
      oh_r23_std[i] <- sd(linesi$oh_r23)
      quants <- hdiofmcmc(linesi$oh_r23, credmass=clim)
      oh_r23_lo[i] <- quants[1]
      oh_r23_hi[i] <- quants[2]
      
    }
    
    data.frame(plateifu, d4000_n, d4000_n_err,
               lick_hd_a, lick_hd_a_err, mgfe,
               d_kpc, d_re,
               tauv_m, tauv_std,
               sigma_mstar_m , sigma_mstar_std , sigma_mstar_lo , sigma_mstar_hi , 
               sigma_sfr_m , sigma_sfr_std , sigma_sfr_lo , sigma_sfr_hi , 
               ssfr_m , ssfr_std , ssfr_lo , ssfr_hi ,     
               relsfr_m , relsfr_std , relsfr_lo , relsfr_hi , 
               tbar_m , tbar_std , tbar_lo , tbar_hi , 
               tbar_lum_m , tbar_lum_std , tbar_lum_lo , tbar_lum_hi , 
               g_i_m , g_i_std , g_i_lo , g_i_hi , 
               tauv_bd_m , tauv_bd_std , tauv_bd_lo , tauv_bd_hi , 
               sigma_logl_ha_m , sigma_logl_ha_std , sigma_logl_ha_lo , sigma_logl_ha_hi , 
               sigma_logl_ha_ctauv_m , sigma_logl_ha_ctauv_std , sigma_logl_ha_ctauv_lo , sigma_logl_ha_ctauv_hi , 
               sigma_logl_ha_ctauv_bd_m , sigma_logl_ha_ctauv_bd_std , sigma_logl_ha_ctauv_bd_lo , sigma_logl_ha_ctauv_bd_hi ,
               eqw_ha_m, eqw_ha_std, eqw_ha_lo, eqw_ha_hi,
               o3hbeta_m , o3hbeta_std , o3hbeta_lo , o3hbeta_hi , 
               o1halpha_m , o1halpha_std , o1halpha_lo , o1halpha_hi , 
               n2halpha_m , n2halpha_std , n2halpha_lo , n2halpha_hi , 
               s2halpha_m , s2halpha_std , s2halpha_lo , s2halpha_hi , 
               r23_m , r23_std , r23_lo , r23_hi , 
               o3n2_m , o3n2_std , o3n2_lo , o3n2_hi , 
               oh_n2_m , oh_n2_std , oh_n2_lo , oh_n2_hi , 
               oh_o3n2_m , oh_o3n2_std , oh_o3n2_lo , oh_o3n2_hi , 
               oh_r23_m , oh_r23_std , oh_r23_lo , oh_r23_hi , 
               bpt=bpt)
}

## estimated mean mass fraction in broad ages bins

sum_binnedmass <- function(mgh_post, ages, ages.bins = c(0.1, 2.5, 5)) {
  T.gyr <- 10^(ages-9)
  ind.bins <- findInterval(ages.bins, T.gyr)
  dims <- dim(mgh_post)
  nsim <- dims[1]
  nfib <- dims[3]
  nt <- length(ages.bins)+2
  mgh.binned <- mgh_post[, c(1, ind.bins, dims[2]), ]
  mdiff <- mgh.binned[, 1:(nt-1),] - mgh.binned[, 2:nt, ]
  mb_m <- apply(mdiff, c(2, 3), mean)
  mb_sd <- apply(mdiff, c(2, 3), sd)
  df <- data.frame(cbind(t(mb_m), t(mb_sd)))
  df[df==0] <- NA
  T.ind <- c(0, T.gyr[ind.bins], T.gyr[length(T.gyr)])
  nt <- length(T.ind)
  bnames <- paste(formatC(T.ind[1:(nt-1)], format="f", digits=1), " < T < ",
                     formatC(T.ind[2:nt], format="f", digits=1), sep="")
  names(df) <- c(paste(bnames, "_m", sep=""), paste(bnames, "_sd", sep=""))
  df
}
  

## computes highest density interval from a sample of representative values,
##   estimated as shortest credible interval.
## arguments:
##   samplevec
##     is a vector of representative values from a probability distribution.
##   credmass
##     is a scalar between 0 and 1, indicating the mass within the credible
##     interval that is to be estimated.
## value:
##   hdilim is a vector containing the limits of the hdi


hdiofmcmc <- function(samplevec , credmass=0.95) {
    sortedpts <- sort(samplevec)
    ciidxinc <- floor(credmass * length(sortedpts))
    ncis <- length(sortedpts) - ciidxinc
    ciwidth <- rep(0 , ncis)
    for (i in 1:ncis) {
        ciwidth[i] <- sortedpts[i + ciidxinc] - sortedpts[i]
    }
    hdimin <- sortedpts[which.min(ciwidth)]
    hdimax <- sortedpts[which.min(ciwidth) + ciidxinc]
    hdilim <- c(hdimin , hdimax)
    return(hdilim)
}

