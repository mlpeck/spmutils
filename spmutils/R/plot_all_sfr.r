plot_all_sfr <- function(sfh, ages, facetby,
                         ptype="instsfr",
                         quants=c(0.025, 0.5, 0.975),
                         logy=TRUE, alpha=0.5) {
  require(ggplot2)
  T <- 10^(ages-9)
  nt <- length(T)
  if (ptype == "instsfr"){
    dt <- diff(ages)/2
    dt <- c(dt[1], dt)
    tl <- 10^(ages-dt)
    th <- 10^(ages+dt)
    dt <- th-tl
    yvals <- apply(sfh, c(1,3), `/`, dt)
    } else {
      yvals <- aperm(sfh, c(2,1,3))
    }
  hasmodel <- which(!is.na(sfh[1,1,]))
  facets <- as.numeric(format(facetby[hasmodel], digits=10))
  y <- apply(yvals[,,hasmodel], c(1, 3), quantile, probs=quants)
  y_med <- as.vector(y[2,,])
  ymin <- as.vector(y[1,,])
  ymax <- as.vector(y[3,,])
  facet <- rep(facets, each=nt)
  
  df <- data.frame(facet=facet, T=T, y_med=y_med, ymin=ymin, ymax=ymax)
  g1 <- ggplot(df) + geom_line(aes(x=T, y=y_med)) +
                geom_ribbon(aes(x=T, ymin=ymin, ymax=ymax), 
                alpha=alpha, show.legend=FALSE)
  g1 <- g1 + scale_x_log10(breaks=c(0.01, .03, 0.1, 0.3, 1, 3, 10))
  g1 <- g1 + xlab("Lookback time (Gyr)") + ylab("SFR")
  if (logy) {
    g1 <- g1 + scale_y_log10()
    g1 <- g1 + facet_wrap(facets = vars(facet))
  } else {
    g1 <- g1 + facet_wrap(facets = vars(facet), scales="free_y")
  }
  plot(g1)
  invisible(list(df=df, graph=g1))
}

plotsfh_binned <- function(sfh, ages, which.bins, ptype="instsfr", quants=c(0.025,.975), log="", ylim=NULL) {
    require(ggplot2)
    
    sfh <- apply(sfh[,,which.bins], c(1, 2), sum, na.rm=TRUE)
    nt <- length(ages)
    ns <- nrow(sfh)
    ages.years <- 10^ages
    yvals <- matrix(0, ns, nt)
    switch(ptype,
        instsfr = {
          dt <- diff(ages)/2
          dt <- c(dt[1], dt)
          tl <- 10^(ages-dt)
          th <- 10^(ages+dt)
          dt <- th-tl
          yvals <- t(t(sfh)/dt)
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
    if (log=="x" || log=="xy") {
        g1 <- g1 + scale_x_log10(breaks=c(0.01, .03, 0.1, 0.3, 1, 3, 10))
    }
    if (log=="y" || log=="xy") {
        g1 <- g1 + scale_y_log10()
    }
    g1
}

