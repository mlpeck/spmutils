## lines that SDSS tracks

speclines <- c(3727.092,3729.875,3869.857,3890.151,3971.123,4102.892,4341.684,4364.435,
4686.991,4862.683,4960.295,5008.240,5413.024,5578.888,6302.046,6313.806,
6365.535,6549.859,6564.614,6585.268,6718.294,6732.678,7137.757)

names(speclines) <- c("O_II_3725","O_II_3727","Ne_III_3868","H_epsilon",
"Ne_III_3970","H_delta","H_gamma","O_III_4363",
"He_II_4685","H_beta","O_III_4959","O_III_5007",
"He_II_5411","O_I_5577","O_I_6300","S_III_6312",
"O_I_6363","N_II_6548","H_alpha","N_II_6583",
"S_II_6716","S_II_6730","Ar_III_7135")


## major emission lines

lambda_f <- c(3727.092,3729.875,3869.857,3971.123,4960.295,
              5008.240,6302.046,6365.535,6549.859,6585.268,
              6718.294,6732.678)
names(lambda_f) <- c("oii_3727","oii_3729","neiii_3869",
                     "neiii_3970","oiii_4959","oiii_5007",
                     "oi_6300","oi_6363","nii_6548","nii_6584",
                     "sii_6717","sii_6730")

lambda_balmer <- c(3890.158,3971.202,4102.899,4341.691,4862.691,6564.632)
names(lambda_balmer) <- c("h_zeta","h_epsilon","h_delta","h_gamma","h_beta","h_alpha")

lambda_em <- sort(c(lambda_f, lambda_balmer))

## forbidden & recombination lines that are usually weak

lambda_weak <- c(4364.435, 4686.991, 5413.024, 5578.888, 6313.806, 7137.757)
names(lambda_weak) <- c("oiii_4363", "heii_4685", "heii_5411", "oi_5577", "siii_6312", "ariii_7135")

## legacy support for reading an sdss spectrum in a csv file

readsdsspec <- function(fname, z, dname="spectra") {
  gdat <- read.csv(file.path(dname,fname), header=TRUE)
  gdat$lambda <- gdat$lambda/(1+z)
  gdat$flux[gdat$ivar==0] <- NA
  gdat$flux[gdat$flux== -9999.] <- NA
  gdat$ivar[gdat$ivar==0] <- NA
  gdat$ivar[gdat$ivar== -9999.] <- NA
  gdat
}

regrid <- function(lambda.out, lib) {
    lambda.in <- lib$lambda
    nc <- ncol(lib)
    lib.out <- matrix(NA, length(lambda.out), nc)
    lib.out[,1] <- lambda.out
    for (i in 2:nc) {
        lib.out[,i] <- approx(lambda.in, lib[,i], xout=lambda.out)$y
    }
    colnames(lib.out) <- colnames(lib)
    data.frame(lib.out)
}


closest <- function(x, vec) which.min(abs(x-vec))


vactoair <- function(lambda) lambda/(1+2.735182e-4+131.4182/lambda^2+2.76249E8/lambda^4)

airtovac <- function(lambda) {
  s2 <- (1e4/lambda)^2
  fac <- 1 + 5.792105e-2/(238.0185-s2)+1.67917e-3/(57.362-s2)
  lambda*fac
}

## Galactic extinction correction from Fitzpatrick (1998): http://arxiv.org/abs/astro-ph/9809387v1
## This is spline fit portion valid from near-UV to near-IR and R=3.1

elaw <- function(lambda, ebv) {
  il <- c(0,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846)
  al <- c(0,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591)
  fai <- splinefun(il,al)
  10^(0.4*fai(10000/lambda)*ebv)
}

fitz <- function(lambda, tauv) {
  il <- c(0,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846)
  al <- c(0,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591)
  fai <- splinefun(il,al)
  exp(-fai(10000/lambda)*tauv)
}

## Calzetti et al. (2000-2001) extinction curve

calzetti.orig <- function(lambda, tauv) {
  lambda <- lambda/10000
  kl <- numeric(length(lambda))
  lb <- lambda < 0.63
  kl[lb] <- 1+(2.659/4.05)*(-2.151285+1.509/lambda[lb]-0.198/(lambda[lb]^2)+0.011/(lambda[lb]^3))
  kl[!lb] <- 1+(2.659/4.05)*(-1.8561715+1.04/lambda[!lb])
  exp(-kl*tauv)
}

## approximation to above using a single function

calzetti <- function(lambda, tauv) {
  ls <- lambda/10000
  k <- -0.2688+0.7958/ls-4.785e-2/(ls^2)-6.033e-3/(ls^3)+7.163e-04/(ls^4)
  exp(-k*tauv)
}

## modified calzetti relation from Salim et al. 2018

calzetti_mod <- function(lambda, tauv, delta=0) {
  ls <- 5500./lambda
  kl <- -0.101771 + 0.549882 * ls + 1.393039 * ls^2 - 1.098615 * ls^3 + 0.260618 * ls^4
  kl <- kl * ls^delta
  exp(-tauv * kl)
}

## galactic extinction law of Cardelli, Clayton, & Mathis 1989, ApJ 345, 245

ccm <- function(lambda, tauv, rv=3.1) {
  x <- 10000/lambda
  y <- x-1.82
  
  kl <- numeric(length(lambda))
  a <- numeric(length(lambda))
  b <- numeric(length(lambda))
  
  lb <- (x <= 1.1)
  
  a[lb] <- 0.574*x[lb]^1.61
  b[lb] <- -0.527*x[lb]^1.61
  
  a[!lb] <- 1 + 0.17699*y[!lb] - 0.50447*y[!lb]^2 - 0.02427*y[!lb]^3 + 0.72085*y[!lb]^4 +
            0.01979*y[!lb]^5 - 0.7753*y[!lb]^6  + 0.32999*y[!lb]^7
  b[!lb] <- 1.41338*y[!lb] + 2.28305*y[!lb]^2 + 1.07233*y[!lb]^3 - 5.38434*y[!lb]^4 -
            0.62251*y[!lb]^5 + 5.30260*y[!lb]^6 - 2.09002*y[!lb]^7
  kl <- a + b/rv
  exp(-tauv * kl)
}


simplescreen <- function(lambda, tauv, delta = 0.) {
    exp(-tauv*(5500./lambda)^(0.7 + delta))
}


## make a fake emission line library

make_emlib <- function(lambda_em, vdisp_em, loglambda, isok, sthresh=5000) {
    vdisp_p <- vdisp_em/(299792.458*log(10))
    n_em <- length(lambda_em)
    nl <- length(loglambda)
    x_em <- matrix(0, nl, n_em)
    for (i in 1:n_em) {
        x_em[,i] <- dnorm(loglambda, mean=log10(lambda_em[i]), sd=vdisp_p)
    }
    in_em <- which(colSums(x_em[isok, ]) >= sthresh)
    x_em <- x_em[,in_em]/10000
    list(x_em=x_em, in_em=in_em)
}

blur.lib <- function(lib.ssp, vdisp, dlogl=1.e-4) {
  c <- 299792.458
  sigma.p <- vdisp/(dlogl*c*log(10))
  kw <- max(round(6*sigma.p), 3)
  if ((kw%%2)==0) kw <- kw+1
  k0 <- (kw %/% 2) + 1
  kernel <- dnorm((1:kw)-k0, sd=sigma.p)
  kernel <- kernel/sum(kernel)
  sp.st <- filter(lib.ssp[,-1], kernel, method="convolution")
  as.matrix(sp.st)
}


dblineratios <- function(db, dthresh=3) {
    o3hbeta <- log10(db$oiii_flux) - log10(db$h_beta_flux)
    o3hbeta[(db$oiii_flux/db$oiii_flux_err < dthresh) | (db$h_beta_flux/db$h_beta_flux_err < dthresh)] <- NA
    n2halpha <- log10(db$nii_6584_flux) -log10(db$h_alpha_flux)
    n2halpha[(db$nii_6584_flux/db$nii_6584_flux_err < dthresh) | (db$h_alpha_flux/db$h_alpha_flux_err < dthresh)] <- NA
    db <- cbind(db, o3hbeta, n2halpha)
    return(db)
}

## galaxy class by emission line strengths, from Kewley et al. 2006 and Schawinski et al. 2007
## only using n2halpha and o3hbeta for now. 

emclass <- function(db=NULL, snthresh=3, PLOT=TRUE) {
  
  
  if (!is.null(db)) {
    attach(db)
    class <- as.factor(rep("no em.",length(n2halpha)))
    levels(class) <- c("no em.", "SF", "Comp.", "AGN", "LINER", "EL")
    
    #starforming
    
    class[(o3hbeta <= 0.61/(n2halpha-0.05)+1.3) & (n2halpha <= 0.05)] <- "SF"
    
    # composite
    
    class[((o3hbeta > 0.61/(n2halpha-0.05)+1.3) | (n2halpha > 0.05)) & 
    (o3hbeta <= .61/(n2halpha-0.47)+1.19)] <- "Comp."
    
    # Seyfert == AGN
    
    class[(o3hbeta > .61/(n2halpha-0.47)+1.19) & (o3hbeta > 1.05*n2halpha+0.45)] <- "AGN"
    class[(n2halpha > 0.47) & (o3hbeta > 1.05*n2halpha+0.45)] <- "AGN"
    
    # Liners
    
    class[(o3hbeta > .61/(n2halpha-0.47)+1.19) & (o3hbeta <= 1.05*n2halpha+0.45)] <- "LINER"
    class[(n2halpha > 0.47) & (o3hbeta <= 1.05*n2halpha+0.45)] <- "LINER"
    
    # weak emission line 
    
    class[((is.na(n2halpha) | is.na(o3hbeta)) & ((oii_flux/oii_flux_err > snthresh) |
    (oiii_flux/oiii_flux_err > snthresh) | 
    (h_alpha_flux/h_alpha_flux_err > snthresh) |
    (nii_6584_flux/nii_6584_flux_err > snthresh)))] <- "EL"
    
    class <- as.factor(class)
    if (PLOT) {
      plot(n2halpha, o3hbeta, pch=20, col=unclass(class))
    }
    detach(db)
  }
  
  if (PLOT) {
    
    #kauffmann et al. curve
    
    curve(1.3+.61/(x-0.05),from= -1.5,to= -0.2,lty=2,add=TRUE)
    
    #kewley 2001 curve
    
    curve(.61/(x-0.47)+1.19, from= -1.5, to= 0.22, add=TRUE)
    
    #Shawinski et al. curve
    
    curve(0.45+1.05*x, from=-0.2, to= 1, add=TRUE)
  }
  
  
  if (!is.null(db)) return(class)
}



## flux ratio at the 4000 A discontinuity. Note 20130908 - weights have been dropped; improves
## correlation with MPA-JHU pipeline estimates. Also corrected lambda weighting

## d4000_n index

d4000n <- function(lambda, flux, ivar, wl4000=c(3850,3950,4000,4100)) {
    wl4000 <- airtovac(wl4000)
    ni <- which(lambda>=wl4000[3] & lambda<=wl4000[4])
    di <- which(lambda>=wl4000[1] & lambda<=wl4000[2])
    dln <- diff(lambda)[ni]
    dld <- diff(lambda)[di]
    fn <- flux[ni]
    if (length(which(!is.na(fn))) < 2) {
        return(list(d4000_n=NA, d4000_n_err=NA))
    }
    if (any(is.na(fn))) {
        fn <- approx(lambda[ni], fn, xout=lambda[ni])$y
    }
    fd <- flux[di]
    if (length(which(!is.na(fd)))<2) {
        return(list(d4000_n=NA, d4000_n_err=NA))
    }
    if (any(is.na(fd))) {
        fd <- approx(lambda[di], fd, xout=lambda[di])$y
    }
    num <- sum(fn*dln*(lambda[ni]^2))
    vn <- sum((dln^2)*(lambda[ni]^4)/ivar[ni])
    den <- sum(fd*dld*(lambda[di]^2))
    vd <- sum((dld^2)*(lambda[di]^4)/ivar[di])
    d4000_n <- num/den
    d4000_n_err <- sqrt(num^2/den^2*(vn/den^2+vd/num^2))
    list(d4000_n=d4000_n, d4000_n_err=d4000_n_err)
}

## lick indices

lickew <- function(lambda, flux, ivar, which.index=1) {
  
  wl <- t(matrix(airtovac(c(4083.5,4122.25,4041.6,4079.75,4128.5,4161,
                                      4142.125,4177.125,4080.125,4117.625,4244.125,4284.125,
                                      4142.125,4177.125,4083.875,4096.375,4244.125,4284.125,
                                      4222.250,4234.750,4211.000,4219.750,4241.000,4251.000,
                                      4281.375,4316.375,4266.375,4282.625,4318.875,4335.125,
                                      4369.125,4420.375,4359.125,4370.375,4442.875,4455.375,
                                      4452.125,4474.625,4445.875,4454.625,4477.125,4492.125,
                                      4514.250,4559.250,4504.250,4514.250,4560.500,4579.250,
                                      4634.000,4720.250,4611.500,4630.250,4742.750,4756.500,
                                      4758.500,4800.000,4742.750,4756.500,4827.875,4847.875,
                                      5069.125,5134.125,4895.125,4957.625,5301.125,5366.125,
                                      5154.125,5196.625,4895.125,4957.625,5301.125,5366.125,
                                      5160.125,5192.625,5142.625,5161.375,5191.375,5206.375,
                                      5245.650,5285.650,5233.150,5248.150,5285.650,5318.150,
                                      5312.125,5352.125,5304.625,5315.875,5353.375,5363.375,
                                      5387.500,5415.000,5376.250,5387.500,5415.000,5425.000,
                                      5445.000,5600.000,5420.000,5442.000,5630.000,5655.000,
                                      5696.625,5720.375,5672.875,5696.625,5722.875,5736.625,
                                      5776.625,5796.625,5765.375,5775.375,5797.875,5811.625,
                                      5876.875,5909.375,5860.625,5875.625,5922.125,5948.125,
                                      5936.625,5994.125,5816.625,5849.125,6038.625,6103.625,
                                      6189.625,6272.125,6066.625,6141.625,6372.625,6415.125,
                                      6357.500,6401.750,6342.125,6356.500,6408.500,6429.750,
                                      6775.000,6900.000,6510.000,6539.250,7017.000,7064.000)),6,24))
  inames <- c("HdeltaA","CN_1","CN_2","Ca4227","G4300",
              "Fe4383","Ca4455","Fe4531","Fe4668","bTiO",
              "Mg_1","Mg_2","Mg_b","Fe5270","Fe5335",
              "Fe5406","aTiO","Fe5709","Fe5782","Na_D",
              "TiO_1","TiO_2","CaH1", "CaH2")
  imag <- rep(FALSE, length(inames))
  imag[c(2,3,11,12,21,22)] <- TRUE
  imag <- imag[which.index]
  wl <- matrix(wl[which.index,], length(which.index), 6)
  ni <- nrow(wl)
  ew <- rep(NA,ni)
  ew_sd <- rep(NA,ni)
  df <- data.frame(lambda=lambda, flux=flux, ivar=ivar)
  for (j in 1:ni) {
    pb <- which(lambda>=wl[j,1] & lambda<=wl[j,2])
    bsb <- which(lambda>=wl[j,3] & lambda<=wl[j,4])
    rsb <- which(lambda>=wl[j,5] & lambda<=wl[j,6])
    if (all(is.na(flux[pb])) || all(is.na(flux[bsb])) || all(is.na(flux[rsb]))) next
    cfit <- lm(flux~lambda, data=df, subset=union(bsb,rsb))
    cont <- predict(cfit, newdata=df[pb,], interval="prediction", level=0.7)
    sd.cont <- (cont[,3]-cont[,2])/2
    dl <- (diff(lambda))[pb]
    num <- flux[pb]
    den <- cont[,1]
    if (imag[j]) {
      ewd <- sum(num/den*dl)
      ew[j] <- -2.5*log10(ewd/(wl[j,2]-wl[j,1]))
      ew_sd[j] <- sqrt(sum(num^2/den^2*(1/(ivar[pb]*num^2)+sd.cont^2/den^2)*dl^2))/
      ewd*2.5/log(10)
    } else {
      ew[j] <- sum((1-num/den)*dl)
      ew_sd[j] <- sqrt(sum(num^2/den^2*(1/(ivar[pb]*num^2)+sd.cont^2/den^2)*dl^2))
    }
  }
  in_names <- inames[which.index]
  in_err_names <- paste(in_names, "err", sep="_")
  retdat <- c(ew, ew_sd)
  names(retdat) <- c(in_names, in_err_names)
  retdat
}

