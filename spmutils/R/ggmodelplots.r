## useful plots
  
plotsfh <- function(sfh, ages, ptype="instsfr", quants=c(0.025,.975), logx=TRUE, logy=FALSE, ylim=NULL) {
    require(ggplot2)
    
    nt <- length(ages)
    ns <- nrow(sfh)
    ages.years <- 10^ages
    yvals <- matrix(0, ns, nt)
    switch(ptype,
        instsfr = {
          dt <- diff(ages)/2
          dt <- c(dt[1], dt)
          tl <- 10^(ages-dt)
          tu <- 10^(ages+dt)
          yvals <- t(t(sfh)/(tu-tl))
        },
        mgh = {
          yvals <- sfh
        },
        avgsfr = {
          yvals <- t(apply(sfh, 1, cumsum)/ages.years)
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
    ages.gyr <- ages.years*10^(-9)
    y <- colMeans(yvals)
    ylims <- apply(yvals, 2, quantile, probs=quants)
    df <- data.frame(ages.gyr=ages.gyr, y=y, ymin=ylims[1,], ymax=ylims[2,])
    g1 <- ggplot(df, aes(x=ages.gyr, y=y)) + geom_line() + 
            geom_ribbon(aes(ymin=ymin, ymax=ymax), color="gray70", alpha=0.5) +
            xlab("T (Gyr)") + ylab(ptype)
    if (!is.null(ylim)) {
        g1 <- g1 + ylim(ylim)
    }
    if (logx) {
        g1 <- g1 + scale_x_log10()
    }
    if (logy) {
        g1 <- g1 + scale_y_log10()
    }
    plot(g1)
    invisible(list(graph=g1, df=df))
}

plotsfhmgh <- function(sfh_post, ages, which.spax) {
  g1 <- plotsfh(sfh_post$sfh[,,which.spax], ages=ages, ptype="instsfr", log="x")$graph
  g1 <- g1 + ggtitle(paste("Bin", which.spax))
  g2 <- plotsfh(sfh_post$mgh[,,which.spax], ages=ages, ptype="mgh", log="x")$graph
  gridExtra::grid.arrange(g1, g2, ncol=1)
}

plotmgh <- function(mgh, ages, quants=c(.025, .975), logx=TRUE, logy=FALSE) {
    require(ggplot2)
    ages.gyr <- 10^(ages-9)
    mgh.mean <- colMeans(mgh)
    if (length(mgh.mean) > length(ages.gyr)) {
      tl <- 1.5*ages[1] - 0.5*ages[2]
      ages.gyr <- c(10^(tl-9), ages.gyr)
    }
    lims <- apply(mgh, 2, quantile, probs=quants)
    df <- data.frame(T=ages.gyr, mgh=mgh.mean, ymin=lims[1,], ymax=lims[2,])
    g1 <- ggplot(df, aes(x=T, y=mgh)) + geom_line() + xlab("T (Gyr)") + 
            ylab("Cumulative mass fraction")
    g1 <- g1 + geom_ribbon(aes(ymin=ymin, ymax=ymax), color="gray70", alpha=0.5)
    if (logx) {
        g1 <- g1 + scale_x_log10()
    }
    if (logy) {
        g1 <- g1 + scale_y_log10()
    }
    plot(g1)
    invisible(list(graph=g1, df=df))
}

addnnmgh <- function(ggraph, nnfits, which.spax, z, ages, mstar, cumfrac=TRUE, color="red", linetype=2) {
  ns <- length(mstar)
  nt <- length(ages)
  nz <- ns/nt
  rmass <- matrix(nnfits$nnfits[which.spax, 1:ns]*mstar, nt, nz)
  rmass <- rowSums(rmass) * cosmo::lum.sol(1, z)
  mgh <- c(sum(rmass), sum(rmass)-cumsum(rmass))
  if (cumfrac) {
    mgh <- mgh/sum(rmass)
  }
  df <- data.frame(t = c(0, 10^(ages-9)), mgh=mgh)
  ggraph + geom_line(aes(x=t, y=mgh), data=df, color=color, linetype=linetype)
}
  
## two or more mass growth histories with same age breaks  

multimgh <- function(..., ages, ids=NULL, quants=c(0.025,0.975),
                     palette="Set1", alpha=0.5, legend=NULL) {
  require(ggplot2)
  ages <- 10^(ages-9)
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
  invisible(list(df=df, graph=g1))
}

## mass growth histories in an array (as returned by batch_sfh for example

multimgh2 <- function(mghs, ages, ids=NULL, quants=c(0.025,0.975),
                     palette="Set1", alpha=0.5, legend=NULL) {
  require(ggplot2)
  ages <- 10^(ages-9)
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
  invisible(list(df=df, graph=g1))
}

## mghs and ages in lists of same length

multimgh3 <- function(mgh_list, ages_list, ids=NULL, quants=c(0.025,0.975),
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
    nt <- length(ages_list[[i]])
    x <- c(x, 1.e-9, 10^(ages_list[[i]]-9))
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
  invisible(list(df=df, graph=g1))
}

## (re)plot posterior predictive fit of spectrum

plotpp <- function(sfit, title=NULL, 
                   quants=c(.025,.975), gcolor="grey70", fcolor="turquoise2") {
    require(ggplot2)
    lambda <- sfit$spm_data$lambda
    gflux <- sfit$spm_data$gflux
    g_std <- sfit$spm_data$g_std
    pp <- extract_post(sfit$stanfit)$gflux_rep
    if (exists("norm_g", sfit)) {
      gflux <- gflux * sfit$norm_g
      g_std <- g_std * sfit$norm_g
    }
    ylims <- apply(pp, 2, quantile, probs=quants)
    df <- data.frame(lambda=lambda, gflux=gflux, fitted=colMeans(pp),
                     ymin=ylims[1,], ymax=ylims[2,],
                     residual = (gflux-colMeans(pp))/g_std,
                     rmin=(gflux-ylims[1,])/g_std, rmax=(gflux-ylims[2,])/g_std)
    g1 <- ggplot(df, aes(x=lambda, y=gflux)) + geom_line(color=gcolor)
    g1 <- g1 + xlab(expression(lambda)) + ylab("flux")
    g1 <- g1 + geom_ribbon(aes(ymin=ymin, ymax=ymax), fill=fcolor, alpha=0.33)
    g1 <- g1 + geom_line(aes(y=fitted), color=fcolor)
    if (!is.null(title)) {
      g1 <- g1 + ggtitle(title)
    }
    g2 <- ggplot(df, aes(x=lambda, y=residual)) + geom_line(color=gcolor)
    g2 <- g2 + xlab(expression(lambda)) + ylab("residual")
    g2 <- g2 + geom_ribbon(aes(ymin=rmin, ymax=rmax), fill=fcolor, alpha=0.33)
    g2 <- g2 + geom_line(aes(y=residual), color=fcolor)
    g3 <- gridExtra::grid.arrange(g1, g2, nrow=2)
#    plot(g3)
    invisible(list(df=df, g1=g1, g2=g2))
}

plotfitted <- function(sfit, quants=c(.025,.975), gcolor="grey70", fcolor="turquoise2") {
    require(ggplot2)
    require(gridExtra)
    lambda <- sfit$spm_data$lambda
    gflux <- sfit$spm_data$gflux
    g_std <- sfit$spm_data$g_std
    if (exists("norm_g", sfit)) {
      gflux <- gflux * sfit$norm_g
      g_std <- g_std * sfit$norm_g
    }
    fitted <- extract_post(sfit$stanfit)$mu_g
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
    invisible(list(g1=g1, g2=g2))
}

## plot a spectrum and nnls fit

replotnn <- function(gdat, dz, nnfit, lib.ssp, 
                   which.spax, rlaw=calzetti,
                   title=NULL) {
    require(ggplot2)
    require(reshape2)
    options("warn" = -1)
    
    flux <- gdat$flux[which.spax, ]
    ivar <- gdat$ivar[which.spax, ]
    lambda <- gdat$lambda
    lambda.rest <- lambda/(1+gdat$meta$z+dz[which.spax])
    logl <- log10(lambda.rest)
    lib.ssp$lambda <- airtovac(lib.ssp$lambda)
    lib.ssp <- regrid(lambda.rest, lib.ssp)
    x.st <- blur.lib(lib.ssp, nnfit$vdisp.st[which.spax])
    n.st <- ncol(x.st)
    allok <- complete.cases(flux, ivar, x.st)
    x.st <- x.st[allok, ]
    lambda.rest <- lambda.rest[allok]
    lambda.em <- lambda_em
    x_em <- make_emlib(lambda.em, nnfit$vdisp.em[which.spax], logl, allok)
    in.em <- x_em$in_em
    x.em <- x_em$x_em
    b <- nnfit$nnfits[which.spax, ]
    b <- b[!is.na(b)]
    b.st <- b[1:n.st]
    b.em <- b[(n.st+1):length(b)]
    fitted <- cbind(x.st*as.vector(rlaw(lambda.rest,nnfit$tauv[which.spax])), x.em[allok,]) %*% b
    fitted.em <- x.em[allok,] %*% b.em
    gflux.net <- flux[allok]-fitted.em
    residual <- (flux[allok]-fitted)*sqrt(ivar[allok])
    tdat <- data.frame(lambda=lambda.rest, flux=flux[allok], fitted=fitted,
                       fitted.em=fitted.em, residual=residual)
    tlong <- melt(tdat, id.vars="lambda")
    val <- c(rep(" obs",length(lambda.rest)), rep("fitted",length(lambda.rest)),
             rep("em", length(lambda.rest)), rep("residual",length(lambda.rest)))
    tlong <- cbind(tlong, val = val)
    tlong$variable[tlong$variable !="residual"] <- "flux"
    base <- qplot(lambda, value, data=tlong, geom="line", xlab=expression(lambda), 
                  ylab="", col=val)
    if (!is.null(title)) {
        base <- base + ggtitle(title)
    }
    add.resid <- facet_grid(variable ~ ., scale="free_y")
    g1 <- base+add.resid
    plot(g1)
    options("warn" = 0)
    g1
}

## replot spectrum and posterior predictions from a stan model

replotpp <- function(gdat, dz, nnfits, sfits, which.spax,
                     prep_data = prep_data_mod,
                     title=NULL, 
                     quants=c(.025,.975), gcolor="grey70", fcolor="turquoise2") {
    require(ggplot2)
    
    stan_dat <- prep_data(gdat, dz, nnfits, which.spax=which.spax)
    
    ## recreate the input spectrum and posterior predictive fits
    
    lambda <- stan_dat$lambda
    gflux <- stan_dat$gflux
    g_std <- stan_dat$g_std
    gflux <- gflux * stan_dat$norm_g
    g_std <- g_std * stan_dat$norm_g
    
    nsim <- nrow(sfits$tauv)
    nl <- stan_dat$nl
    tauv <- sfits$tauv[,which.spax]
    delta <- sfits$delta[,which.spax]
    att <- matrix(0, nl, nsim)
    for (i in 1:nsim) {
      att[, i] <- calzetti_mod(lambda, tauv[i], delta[i])
    }
    mu_g <- stan_dat$norm_g * tcrossprod(stan_dat$sp_st, sfits$b_st[,,which.spax]) * att + 
                      stan_dat$norm_g * stan_dat$norm_em * tcrossprod(stan_dat$sp_em, sfits$b_em[,stan_dat$in_em,which.spax])
    pp <- mu_g + matrix(rnorm(nl*nsim, sd=rep(g_std, nsim)), nl, nsim)
    
    ylims <- apply(pp, 1, quantile, probs=quants)
    df <- data.frame(lambda=lambda, gflux=gflux, fitted=rowMeans(pp),
                     ymin=ylims[1,], ymax=ylims[2,],
                     residual = (gflux-rowMeans(pp))/g_std,
                     rmin=(gflux-ylims[1,])/g_std, rmax=(gflux-ylims[2,])/g_std)
    g1 <- ggplot(df, aes(x=lambda, y=gflux)) + geom_line(color=gcolor)
    g1 <- g1 + xlab(expression(lambda)) + ylab("flux")
    g1 <- g1 + geom_ribbon(aes(ymin=ymin, ymax=ymax), fill=fcolor, alpha=0.33)
    g1 <- g1 + geom_line(aes(y=fitted), color=fcolor)
    if (!is.null(title)) {
      g1 <- g1 + ggtitle(title)
    }
    g2 <- ggplot(df, aes(x=lambda, y=residual)) + geom_line(color=gcolor)
    g2 <- g2 + xlab(expression(lambda)) + ylab("residual")
    g2 <- g2 + geom_ribbon(aes(ymin=rmin, ymax=rmax), fill=fcolor, alpha=0.33)
    g2 <- g2 + geom_line(aes(y=residual), color=fcolor)
    g3 <- gridExtra::grid.arrange(g1, g2, nrow=2)
    plot(g3)
    invisible(list(df=df, g1=g1, g2=g2))
}
