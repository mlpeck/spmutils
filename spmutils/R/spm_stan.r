stanfit_one <- function(gdat, dz, nnfits, which.spax, 
                        lib.ssp, age, Z,
                        emlines=lambda_em,
                        stan_model=NULL,
                        stan_file="spm_dust_norm.stan", stan_filedir="~/spmcode/",
                        iter=500, warmup=250, thin=1, chains=4, 
                        OP=FALSE, ...) {
    
    require(rstan)
    i <- which.spax[1]
    j <- which.spax[2]
    z <- gdat$meta$z+dz[i, j]
    flux <- gdat$flux[i, j, ]
    ivar <- gdat$ivar[i, j, ]
    lambda.rest <- gdat$lambda/(1+z)
    logl <- log10(lambda.rest)
    nsim <- (iter-warmup)*chains
    T.gyr <- 10^(age-9)
        
    lib.ssp$lambda <- airtovac(lib.ssp$lambda)
    lib.st <- regrid(lambda.rest, lib.ssp)
    x.st <- blur.lib(lib.st, nnfits$vdisp.st[i, j])
    n.st <- ncol(x.st)
    nz <- length(Z)
    nt <- n.st/nz
    ind.young <- (0:(nz-1))*nt+1
    allok <- complete.cases(flux, ivar, x.st)
    nl <- length(which(allok))
    
    x.em <- make_emlib(emlines, nnfits$vdisp.em[i, j], logl, allok)
    in.em <- x.em$in_em
    x.em <- x.em$x_em
    n.em <- ncol(x.em)
    
    tauv <- nnfits$tauv[i, j]
    
    b.st <- nnfits$nnfits[i, j, 1:n.st]
    b.em <- nnfits$nnfits[i, j, (n.st+1):(n.st+n.em)]
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

                        
stanfit_batch <- function(gdat, dz, nnfits,
                        lib.ssp, age, Z,
                        emlines=lambda_em,
                        stan_file="spm_dust_norm.stan", stan_filedir="~/spmcode/",
                        iter=500, warmup=250, chains=4,
                        OP=FALSE,
                        start=NULL, end=NULL, fpart="bfits.rda", ...) {
    dims <- dim(gdat$flux)
    if (dims[2] > 1) stop("Use stacked data!")
    dz <- dz$dz
    nsim <- (iter-warmup)*chains
    nt <- length(age)
    nr <- dims[1]
    n_st <- ncol(lib.ssp)-1
    n_em <- length(emlines)
    nl <- length(gdat$lambda)
    smodel <- rstan::stan_model(file.path(stan_filedir, stan_file))
    if (is.null(start) || !file.exists(fpart)) {
      b_st <- array(NA, dim=c(nsim, n_st, nr))
      b_em <- array(NA, dim=c(nsim, n_em, nr))
      tauv <- matrix(NA, nrow=nsim, ncol=nr)
      in_em <- matrix(NA, nrow=n_em, ncol=nr)
      walltime <- rep(NA, nr)
      divergences <- rep(NA, nr)
      max_treedepth <- rep(NA, nr)
      norm_st <- rep(NA, nr)
      norm_em <- rep(NA, nr)
      start <- 1
    } else {
        load(fpart)
    }
    if (is.null(end)) {
      end <- nr
    }
    for (i in start:end) {
        if (is.na(dz[i, 1]) || is.na(nnfits$Mstar[i, 1])) next
        sfit <- stanfit_one(gdat, dz, nnfits, which.spax=c(i, 1),
                               lib.ssp, age=age, Z=Z,
                               emlines=emlines,
                               stan_model=smodel,
                               iter = iter, warmup = warmup, chains = chains, OP=OP, ...)
        plot(plotpp(sfit)+ggtitle(paste("fiber =", i)))
        post <- extract(sfit$stanfit)
        b_st[,,i] <- post$b_st
        b_em[,sfit$in.em,i] <- post$b_em
        tauv[,i] <- post$tauv
        in_em[sfit$in.em,i] <- sfit$in.em
        walltime[i] <- max(rowSums(get_elapsed_time(sfit$stanfit)))
        sp <- get_sampler_params(sfit$stanfit)
        divergences[i] <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
        max_treedepth[i] <- max(sapply(sp, function(x) max(x[, "treedepth__"])))
        norm_st[i] <- sfit$norm.st
        norm_em[i] <- sfit$norm.em
        save(b_st, b_em, tauv, in_em,
                  walltime, divergences, max_treedepth,
                  norm_st, norm_em, file = fpart)
        rm(sfit, post, sp)
    }
    list(b_st=b_st, b_em=b_em, tauv=tauv, in_em=in_em,
                  walltime=walltime, divergences=divergences, max_treedepth=max_treedepth,
                  norm_st=norm_st, norm_em=norm_em)
}

## star formation history, mass growth history, etc.

get_sfh <- function(b_st, norm_st, age, mstar, z, fibersinbin=1, tsf=0.1) {
  nsim <- nrow(b_st)
  nt <- length(age)
  nz <- ncol(b_st)/nt
  T.gyr <- 10^(age-9)
  isf <- which.min(abs(tsf-T.gyr))
  binarea <- log10(pi*fibersinbin*cosmo::ascale(z)^2)
  b_st <- b_st*norm_st*cosmo::lum.sol(1, z)
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

batch_sfh <- function(b_st, norm_st, age, mstar, z, fibersinbin, tsf=0.1) {
  dims <- dim(b_st)
  nsim <- dims[1]
  nt <- length(age)
  nz <- dims[2]/nt
  nf <- length(norm_st)
  
  sfh_post <- array(NA, dim=c(nsim, nt, nf))
  mgh_post <- array(0, dim=c(nsim, nt+1, nf))
  totalmg_post <- matrix(0, nsim, nt+1)
  sigma_mstar <- matrix(NA, nsim, nf)
  sigma_sfr <- matrix(NA, nsim, nf)
  ssfr <- matrix(NA, nsim, nf)
  relsfr <- matrix(NA, nsim, nf)
  
  for (i in 1:nf) {
    if (is.na(b_st[1, 1, i])) next
    sfi <- get_sfh(b_st[,,i], norm_st[i], age=age, mstar=mstar, z=z, fibersinbin=fibersinbin[i])
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

get_proxies <- function(b_st, age, Z, gri.ssp) {
  nz <- length(Z)
  nt <- length(age)
  T.gyr <- 10^(age-9)
  tbar <- log10((b_st %*% rep(T.gyr, nz))/rowSums(b_st)) + 9
  tbar_lum <- log10((b_st %*% (rep(T.gyr, nz) * gri.ssp["r",]))/
                   ((b_st %*% gri.ssp["r",]))) + 9
  g_i <- 2.5*log10(b_st %*% gri.ssp["i",]) - 2.5*log10(b_st %*% gri.ssp["g",])
  data.frame(tbar=tbar, tbar_lum=tbar_lum, g_i=g_i)
}

## emission line fluxes, luminosity density, equivalent width

get_em <- function(b_em, b_st, norm_em, norm_st, emlines, z, lib.ssp, fibersinbin=1, ew_width=15) {
  ne <- ncol(b_em)
  nsim <- nrow(b_em)
  binarea <- log10(pi*fibersinbin*cosmo::ascale(z)^2)
  em.mult <- emlines*log(10)

  flux_em <- matrix(NA, nsim, ne)
  sigma_logl_em <- matrix(NA, nsim, ne)
  ew_em <- matrix(NA, nsim, ne)
  
  flux_em <- t(t(b_em)*em.mult)*norm_em
  flux_em[flux_em < 0] <- 0
  sigma_logl_em <- cosmo::loglum.ergs(flux_em, z) - binarea
  
  mu_st <- tcrossprod(b_st, as.matrix(lib.ssp[, -1])) * norm_st
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

batch_em <- function(b_em, b_st, norm_em, norm_st, in_em, emlines, z, lib.ssp, fibersinbin, ew_width=15) {
  nsim <- dim(b_em)[1]
  ne <- length(emlines)
  nf <- length(fibersinbin)
  
  flux_em <- array(NA, dim=c(nsim, ne, nf))
  sigma_logl_em <- array(NA, dim=c(nsim, ne, nf))
  ew_em <- array(NA, dim=c(nsim, ne, nf))
  
  for (i in 1:nf) {
    if (is.na(b_st[1, 1, i])) next
    in_em_i <- in_em[!is.na(in_em[,i]), i]
    emi <- get_em(b_em[,in_em_i,i], b_st[,,i], norm_em[i], norm_st[i], emlines=emlines[in_em_i], z=z, lib.ssp=lib.ssp,
                  fibersinbin=fibersinbin[i], ew_width=ew_width)
    flux_em[,in_em_i,i] <- emi$flux_em
    sigma_logl_em[,in_em_i,i] <- emi$sigma_logl_em
    ew_em[,in_em_i,i] <- emi$ew_em
  }
  dimnames(flux_em)[[2]] <- names(emlines)
  dimnames(sigma_logl_em)[[2]] <- names(emlines)
  dimnames(ew_em)[[2]] <- names(emlines)
  list(flux_em=flux_em, sigma_logl_em=sigma_logl_em, ew_em=ew_em)
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
  
sum_batchfits <- function(gdat, nnfits, sfits, bptclass, lib.ssp, age, Z, gri.ssp, mstar, 
                            emlines=lambda_em, alaw=calzetti, intr_bd=2.86, clim=0.95) {
    
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
    nt <- length(age)
    fiberarea <- pi*cosmo::ascale(gdat$meta$z)^2
    plateifu <- rep(gdat$meta$plateifu, nf)
    
    ## projected distance in kpc and relative to effective radius
    
    d_kpc <- cosmo::ascale(gdat$meta$z)*sqrt(gdat$xpos^2+gdat$ypos^2)
    d_re <- sqrt(gdat$xpos^2+gdat$ypos^2)/
    drpcat$nsa_petro_th50[match(gdat$meta$plateifu,drpcat$plateifu)]
    
    
    ## stuff taken from nnfits
    
    d4000_n <- nnfits$d4000_n
    d4000_n_err <- nnfits$d4000_n_err
    lick_hd_a <- nnfits$lick[,,'HdeltaA']
    lick_hd_a_err <- nnfits$lick.err[,,'HdeltaA_err']
    mgfe <- sqrt(nnfits$lick[,,'Mg_b']*(0.72*nnfits$lick[,,'Fe5270']+0.28*nnfits$lick[,,'Fe5335']))
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
    
    
    oii_flux <- oii_flux_err <- numeric(nf)
    
    sfh_all <- batch_sfh(sfits$b_st, sfits$norm_st, age=age, mstar=mstar, z=gdat$meta$z, fibersinbin=fibersinbin)
    
    for (i in 1:4) {
      assign(paste(varnames[i], suffixes[1], sep="_"), colMeans(sfh_all[[varnames[i]]]))
      assign(paste(varnames[i], suffixes[2], sep="_"), apply(sfh_all[[varnames[i]]], 2, sd))
    }
      
        
    em_all <- batch_em(sfits$b_em, sfits$b_st, sfits$norm_em, sfits$norm_st, sfits$in_em, emlines=emlines, z=gdat$meta$z, 
                       lib.ssp=lib.ssp, fibersinbin=fibersinbin)
    
    
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
      
      proxi <- get_proxies(sfits$b_st[,,i], age, Z, gri.ssp)
      
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
               bpt=bptclass)
}

## estimated mean mass fraction in broad age bins

sum_binnedmass <- function(mgh.post, ages, ages.bins = c(0.1, 2.5, 5)) {
  T.gyr <- 10^(ages-9)
  ind.bins <- findInterval(ages.bins, T.gyr)
  dims <- dim(mgh.post)
  nsim <- dims[1]
  nfib <- dims[3]
  nt <- length(ages.bins)+2
  mgh.binned <- mgh.post[, c(1, ind.bins, dims[2]), ]
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
  
## useful plots
  
plotsfh <- function(sfh, age, ptype="instsfr", quants=c(0.025,.975), log="", ylim=NULL) {
    require(ggplot2)
    
    nt <- ncol(sfh)
    ns <- nrow(sfh)
    age.years <- 10^age
    yvals <- matrix(0, ns, nt)
    switch(ptype,
        avgsfr = {
            yvals <- t(apply(sfh, 1, cumsum)/age.years)
        },
        instsfr = {
            dt <- diff(c(0,age.years))
            yvals <- t(t(sfh)/dt)
        },
        cumfrac = {
            tm <- rowSums(sfh)
            yvals <- t(apply(sfh, 1, cumsum))/tm
        },
        cumfrac0 = {
            tm <- rowSums(sfh)
            yvals <- 1-t(apply(sfh, 1, cumsum))/tm
        }
    )
    age.gyr <- age.years*10^(-9)
    y <- colMeans(yvals)
    ylims <- apply(yvals, 2, quantile, probs=quants)
    df <- data.frame(age.gyr=age.gyr, y=y, ymin=ylims[1,], ymax=ylims[2,])
    g1 <- ggplot(df, aes(x=age.gyr, y=y)) + geom_line() + 
            geom_ribbon(aes(ymin=ymin, ymax=ymax), color="gray70", alpha=0.5) +
            xlab("T (Gyr)") + ylab(ptype)
    if (!is.null(ylim)) {
        g1 <- g1 + ylim(ylim)
    }
    if (log=="x" || log=="xy") {
        g1 <- g1 + scale_x_log10()
    }
    if (log=="y" || log=="xy") {
        g1 <- g1 + scale_y_log10()
    }
    g1
}

plotmgh <- function(mgh, age, quants=c(.025, .975), log="", ylim=c(0,1)) {
    require(ggplot2)
    age.gyr <- 10^(age-9)
    mgh.mean <- colMeans(mgh)
    if (length(mgh.mean) > length(age.gyr)) {
      age.gyr <- c(1.e-9, age.gyr)
    }
    lims <- apply(mgh, 2, quantile, probs=quants)
    df <- data.frame(T=age.gyr, mgh=mgh.mean, ymin=lims[1,], ymax=lims[2,])
    g1 <- ggplot(df, aes(x=T, y=mgh)) + geom_line() + xlab("T (Gyr)") + 
            ylab("Cumulative mass fraction")
    g1 <- g1 + geom_ribbon(aes(ymin=ymin, ymax=ymax), color="gray70", alpha=0.5) +
            ylim(ylim)
    if (log=="x" || log=="xy") {
        g1 <- g1 + scale_x_log10()
    }
    if (log=="y" || log=="xy") {
        g1 <- g1 + scale_y_log10()
    }
    g1
}

addnnmgh <- function(ggraph, nnfits, which.spax, z, age, mstar, cumfrac=TRUE, color="red", linetype=2) {
  if (length(which.spax)==1) {
    which.spax <- c(which.spax,1)
  }
  ns <- length(mstar)
  nt <- length(age)
  nz <- ns/nt
  rmass <- matrix(nnfits$nnfits[which.spax[1],which.spax[2],1:ns]*mstar, nt, nz)
  rmass <- rowSums(rmass) * cosmo::lum.sol(1, z)
  mgh <- c(sum(rmass), sum(rmass)-cumsum(rmass))
  if (cumfrac) {
    mgh <- mgh/sum(rmass)
  }
  df <- data.frame(t = c(0, 10^(age-9)), mgh=mgh)
  ggraph + geom_line(aes(x=t, y=mgh), data=df, color=color, linetype=linetype)
}
  
  

multimgh <- function(..., age, ids=NULL, quants=c(0.025,0.975),
                     palette="Set1", alpha=0.5, legend=NULL) {
  require(ggplot2)
  ages <- 10^(age-9)
  mghs <- list(...)
  if (dim(mghs[[1]])[2] > length(ages)) {
    ages <- c(1.e-9, ages)
  }
  N <- length(list(...))
  nt <- length(ages)
  x <- numeric(0)
  y <- numeric(0)
  ymin <- numeric(0)
  ymax <- numeric(0)
  if (is.null(ids)) {
    ids <- as.factor(1:N)
  }
  id <- rep(ids, each=nt)
  levels(id) <- ids
  for (i in 1:N) {
    x <- c(x,ages)
    y <- c(y, colMeans(mghs[[i]]))
    ylims <- apply(mghs[[i]], 2, quantile, probs=quants)       
    ymin <- c(ymin, ylims[1,])
    ymax <- c(ymax, ylims[2,])
  }
  df <- data.frame(id=id, x=x, y=y, ymin=ymin, ymax=ymax)
  g1 <- ggplot(df) + geom_line(aes(x=x, y=y, color=id)) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, fill=id), 
                alpha=alpha, show.legend=FALSE)
  g1 <- g1 + scale_color_brewer(type="qual", palette=palette, breaks=ids)
  g1 <- g1 + scale_fill_brewer(type="qual",palette=palette, breaks=ids)
  g1 <- g1 + xlab("Lookback time (Gyr)") + ylab("Cumulative mass fraction")
  if (!is.null(legend)) {
    g1 <- g1 + labs(color=legend)
  }
  plot(g1)
  list(df=df, graph=g1)
}



multimgh2 <- function(mghs, age, ids=NULL, quants=c(0.025,0.975),
                     palette="Set1", alpha=0.5, legend=NULL) {
  require(ggplot2)
  ages <- 10^(age-9)
  if (dim(mghs)[2] > length(ages)) {
    ages <- c(1.e-9, ages)
  }
  ok2plot <- which(!is.na(mghs[1,1,]))
  nt <- length(ages)
  x <- numeric(0)
  y <- numeric(0)
  ymin <- numeric(0)
  ymax <- numeric(0)
  if (is.null(ids)) {
    ids <- as.factor(ok2plot)
  }
  id <- rep(ids, each=nt)
  levels(id) <- ids
  for (i in ok2plot) {
    x <- c(x,ages)
    y <- c(y, colMeans(mghs[,,i]))
    ylims <- apply(mghs[,,i], 2, quantile, probs=quants)       
    ymin <- c(ymin, ylims[1,])
    ymax <- c(ymax, ylims[2,])
  }
  df <- data.frame(id=id, x=x, y=y, ymin=ymin, ymax=ymax)
  g1 <- ggplot(df) + geom_line(aes(x=x, y=y, color=id)) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, fill=id), 
                alpha=alpha, show.legend=FALSE)
  if (length(ok2plot) <= 9) {
    g1 <- g1 + scale_color_brewer(type="qual", palette=palette, breaks=ids)
    g1 <- g1 + scale_fill_brewer(type="qual",palette=palette, breaks=ids)
  } else {
    g1 <- g1 + scale_color_manual(values=viridis::viridis(length(ok2plot)))
    g1 <- g1 + scale_fill_manual(values=viridis::viridis(length(ok2plot)))
  }
  g1 <- g1 + xlab("Lookback time (Gyr)") + ylab("Cumulative mass fraction")
  if (!is.null(legend)) {
    g1 <- g1 + labs(color=legend)
  }
  plot(g1)
  list(df=df, graph=g1)
}

## mghs and ages in lists of same length

multimgh3 <- function(mgh_list, age_list, ids=NULL, quants=c(0.025,0.975),
                     palette="Set1", alpha=0.5, legend=NULL) {
  require(ggplot2)
  N <- length(mgh_list)
  x <- numeric(0)
  y <- numeric(0)
  ymin <- numeric(0)
  ymax <- numeric(0)
  if (is.null(ids)) {
    ids <- as.factor(1:N)
  }
  id <- factor(levels=ids)
  for (i in 1:N) {
    nt <- length(age_list[[i]])
    x <- c(x, 1.e-9, 10^(age_list[[i]]-9))
    y <- c(y, colMeans(mgh_list[[i]]))
    ylims <- apply(mgh_list[[i]], 2, quantile, probs=quants)       
    ymin <- c(ymin, ylims[1,])
    ymax <- c(ymax, ylims[2,])
    id <- c(id, rep(ids[i], nt+1))
  }
  id <- as.factor(id)
  df <- data.frame(id=id, x=x, y=y, ymin=ymin, ymax=ymax)
  g1 <- ggplot(df) + geom_line(aes(x=x, y=y, color=id)) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, fill=id), 
                alpha=alpha, show.legend=FALSE)
  g1 <- g1 + scale_color_brewer(type="qual", palette=palette, breaks=ids)
  g1 <- g1 + scale_fill_brewer(type="qual",palette=palette, breaks=ids)
  g1 <- g1 + xlab("Lookback time (Gyr)") + ylab("Cumulative mass fraction")
  if (!is.null(legend)) {
    g1 <- g1 + labs(color=legend)
  }
  plot(g1)
  list(df=df, graph=g1)
}

plotpp <- function(sfit, quants=c(.025,.975), gcolor="grey70", fcolor="turquoise2") {
    require(ggplot2)
    lambda <- sfit$spm_data$lambda
    gflux <- sfit$spm_data$gflux
    pp <- rstan::extract(sfit$stanfit)$gflux_rep
    ylims <- apply(pp, 2, quantile, probs=quants)
    df <- data.frame(lambda=lambda, gflux=gflux, fitted=colMeans(pp),
                     ymin=ylims[1,], ymax=ylims[2,])
    g1 <- ggplot(df, aes(x=lambda, y=gflux)) + geom_line(color=gcolor)
    g1 <- g1 + xlab(expression(lambda)) + ylab("flux")
    g1 <- g1 + geom_ribbon(aes(ymin=ymin, ymax=ymax), fill=fcolor, alpha=0.33)
    g1 <- g1 + geom_line(aes(y=fitted), color=fcolor)
    g1
}

plotfitted <- function(sfit, quants=c(.025,.975), gcolor="grey70", fcolor="turquoise2") {
    require(ggplot2)
    require(gridExtra)
    lambda <- sfit$spm_data$lambda
    gflux <- sfit$spm_data$gflux
    g_std <- sfit$spm_data$g_std
    fitted <- rstan::extract(sfit$stanfit)$mu_g
    res <- t((gflux-t(fitted))/g_std)
    ylims <- apply(fitted, 2, quantile, probs=quants)
    rlims <- apply(res, 2, quantile, probs=quants)
    df <- data.frame(lambda=lambda, gflux=gflux, fitted=colMeans(fitted),
                     ymin=ylims[1,], ymax=ylims[2,],
                     rmean=colMeans(res),rmin=rlims[1,],rmax=rlims[2,]
                    )
    g1 <- ggplot(df, aes(x=lambda, y=gflux)) + geom_line(color=gcolor)
    g1 <- g1 + xlab(expression(lambda)) + ylab("flux")
    g1 <- g1 + geom_ribbon(aes(ymin=ymin, ymax=ymax), fill=fcolor)
    g1 <- g1 + geom_line(aes(y=fitted), color=fcolor)
    g2 <- ggplot(df, aes(x=lambda, y=rmean)) + geom_line(color=gcolor)
    g2 <- g2 + xlab(expression(lambda)) + ylab("sd")
    g2 <- g2 + geom_line(aes(y=rmean), color=fcolor)
    g2 <- g2 + geom_ribbon(aes(ymin=rmin, ymax=rmax), fill=fcolor)
    plot(grid.arrange(g1, g2, nrow=2, heights=c(3,1)))
    list(g1=g1, g2=g2)
}

## simple 2D spatial interpolation

fillpoly <- function(ra, dec, zvals, dxy=0.5, min_ny=100, usefields=TRUE) {
    ny <- max(round((max(dec)-min(dec))*3600/dxy), min_ny)
    nx <- ny
    if (usefields) {
      xy <- cbind(ra, dec)
      zxy <- fields::Tps(x=xy, Y=zvals, lon.lat=TRUE, miles=FALSE, lambda=0)
      zsurf <- fields::predictSurface(zxy, nx=nx, ny=ny)
    } else {
      allok <- complete.cases(ra, dec, zvals)
      zsurf <- akima::interp(x=ra[allok], y=dec[allok], z=zvals[allok], nx=nx, ny=ny)
    }
    list(ra=zsurf$x, dec=zsurf$y, z=zsurf$z)
}

mapsummaries <- function(gdat, nnfits, stanfits, dz, ages, outnames=
                    c("vrel","d4000_n", "lick_hd_a", "lick_nad", "lick_mgfe",
                      "g_i", "Sigma_mstar_m", "Sigma_sfr_m",
                      "ssfr_m", "tbar_m", "tbar_lum_m", 
                      "zbar_m", "tauv_m", "tauv_bd_m",
                      "logl_halpha", "logl_oii", "logl_oiii", "logl_nii",
                      "oiiihbeta", "oihalpha", "niihalpha", "siihalpha"),
                      which.mgh=c(5,11,14,24), bdint=2.86, alaw=calzetti) {
    t.gyr <- 10^(ages-9)
    cmf <- paste("cmf_", formatC(t.gyr[which.mgh], digits=2, format="f", flag="0"), sep="")
    nsum <- length(outnames) + length(cmf)
    fiberarea <- pi*cosmo::ascale(gdat$meta$z)^2
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, nnfits$d4000_n)
    ra <- zvals$ra
    dec <- zvals$dec
    maps <- array(0, dim=c(length(ra), length(dec), nsum))
    dimnames(maps)[[3]] <- c(outnames, cmf)
    maps[,,"d4000_n"] <- zvals$z
    vf <- cosmo::vrel(gdat$meta$z+dz, gdat$meta$z)
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, vf)
    maps[,,"vrel"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, nnfits$lick[,,'HdeltaA'])
    maps[,,"lick_hd_a"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, nnfits$lick[,,'Na_D'])
    maps[,,"lick_nad"] <- zvals$z
    mgfe <- sqrt(nnfits$lick[,,'Mg_b']*(0.72*nnfits$lick[,,'Fe5270']+0.28*nnfits$lick[,,'Fe5335']))
    mgfe[is.nan(mgfe)] <- NA
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, mgfe)
    maps[,,"lick_mgfe"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, nnfits$gri[,,"g"]-nnfits$gri[,,"i"])
    maps[,,"g_i"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, stanfits$sumpost$Mstar_m)
    maps[,,"Sigma_mstar_m"] <- zvals$z - log10(fiberarea)
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, stanfits$sumpost$sfr_m)
    maps[,,"Sigma_sfr_m"] <- zvals$z - log10(fiberarea)
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, stanfits$sumpost$ssfr_m)
    maps[,,"ssfr_m"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, stanfits$sumpost$tbar_m)
    maps[,,"tbar_m"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, stanfits$sumpost$tbar_lum_m)
    maps[,,"tbar_lum_m"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, stanfits$sumpost$zbar_m)
    maps[,,"zbar_m"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, stanfits$sumpost$tauv_m)
    maps[,,"tauv_m"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(stanfits$logl.em[,'h_alpha',]))
    maps[,,"logl_halpha"] <- zvals$z
    oii <- log10(10^(stanfits$logl.em[,'oii_3727',]) + 10^(stanfits$logl.em[,'oii_3729',]))
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(oii))
    maps[,,"logl_oii"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(stanfits$logl.em[,'oiii_5007',]))
    maps[,,"logl_oiii"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(stanfits$logl.em[,'nii_6584',]))
    maps[,,"logl_nii"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(stanfits$logl.em[,'oiii_5007',] -
                                                    stanfits$logl.em[,'h_beta',]))
    maps[,,"oiiihbeta"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(stanfits$logl.em[,'oi_6300',] -
                                                    stanfits$logl.em[,'h_alpha',]))
    maps[,,"oihalpha"] <- zvals$z
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(stanfits$logl.em[,'nii_6548',] -
                                                    stanfits$logl.em[,'h_alpha',]))
    maps[,,"niihalpha"] <- zvals$z
    sii <- log10(10^(stanfits$logl.em[,'sii_6717',]) + 10^(stanfits$logl.em[,'sii_6730',]))
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(sii -
                                                    stanfits$logl.em[,'h_alpha',]))
    maps[,,"siihalpha"] <- zvals$z
    tauv.bd <- function(logl.em, bdint, alaw) {
        bd <- 10^(logl.em[,'h_alpha',]-logl.em[,'h_beta',])
        bd[!is.finite(bd)] <- NA
        tauv <- log(bd/bdint)/(log(alaw(6562.8,1))-log(alaw(4861.3,1)))
        tauv[tauv<0] <- 0
        tauv
    }
    zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(tauv.bd(stanfits$logl.em, bdint, alaw)))
    maps[,,"tauv_bd_m"] <- zvals$z
    for (i in seq_along(which.mgh)) {
        zvals <- fillpoly(gdat$ra.f, gdat$dec.f, colMeans(stanfits$mgh.post[,which.mgh[i],]))
        maps[,,cmf[i]] <- zvals$z
    }
    returns <- list(ra=ra, dec=dec, maps=maps)
    rra <- rev(range(returns$ra))
    rdec <- range(returns$dec)
    for (i in 1:nsum) {
        g1 <- ggimage(maps[,,i], ra, dec,
                      xlab=expression(alpha), ylab=expression(delta),
                      legend=(dimnames(maps)[[3]])[i], title=(dimnames(maps)[[3]])[i]) +
            xlim(rra) + ylim(rdec)
        plot(g1)
    }
    returns
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

