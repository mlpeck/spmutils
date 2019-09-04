## useful plots
  
plotsfh <- function(sfh, ages, ptype="instsfr", quants=c(0.025,.975), log="", ylim=NULL) {
    require(ggplot2)
    
    nt <- ncol(sfh)
    ns <- nrow(sfh)
    ages.years <- 10^ages
    yvals <- matrix(0, ns, nt)
    switch(ptype,
        avgsfr = {
            yvals <- t(apply(sfh, 1, cumsum)/ages.years)
        },
        instsfr = {
            dt <- diff(c(0,ages.years))
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
    if (log=="x" || log=="xy") {
        g1 <- g1 + scale_x_log10()
    }
    if (log=="y" || log=="xy") {
        g1 <- g1 + scale_y_log10()
    }
    g1
}

plotmgh <- function(mgh, ages, quants=c(.025, .975), log="", ylim=c(0,1)) {
    require(ggplot2)
    ages.gyr <- 10^(ages-9)
    mgh.mean <- colMeans(mgh)
    if (length(mgh.mean) > length(ages.gyr)) {
      ages.gyr <- c(1.e-9, ages.gyr)
    }
    lims <- apply(mgh, 2, quantile, probs=quants)
    df <- data.frame(T=ages.gyr, mgh=mgh.mean, ymin=lims[1,], ymax=lims[2,])
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

addnnmgh <- function(ggraph, nnfits, which.spax, z, ages, mstar, cumfrac=TRUE, color="red", linetype=2) {
  if (length(which.spax)==1) {
    which.spax <- c(which.spax,1)
  }
  ns <- length(mstar)
  nt <- length(ages)
  nz <- ns/nt
  rmass <- matrix(nnfits$nnfits[which.spax[1],which.spax[2],1:ns]*mstar, nt, nz)
  rmass <- rowSums(rmass) * cosmo::lum.sol(1, z)
  mgh <- c(sum(rmass), sum(rmass)-cumsum(rmass))
  if (cumfrac) {
    mgh <- mgh/sum(rmass)
  }
  df <- data.frame(t = c(0, 10^(ages-9)), mgh=mgh)
  ggraph + geom_line(aes(x=t, y=mgh), data=df, color=color, linetype=linetype)
}
  
  

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
  list(df=df, graph=g1)
}



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
  list(df=df, graph=g1)
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
  list(df=df, graph=g1)
}

plotpp <- function(sfit, quants=c(.025,.975), gcolor="grey70", fcolor="turquoise2") {
    require(ggplot2)
    lambda <- sfit$spm_data$lambda
    gflux <- sfit$spm_data$gflux
    pp <- rstan::extract(sfit$stanfit)$gflux_rep
    if (exists("norm_g", sfit)) {
      gflux <- gflux * sfit$norm_g
    }
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
    if (exists("norm_g", sfit)) {
      gflux <- gflux * sfit$norm_g
      g_std <- g_std * sfit$norm_g
    }
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

