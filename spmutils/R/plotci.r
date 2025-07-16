plotci <- function(df, x, y, x.sd=NULL, xmin=NULL, xmax=NULL, y.sd=NULL, ymin=NULL, ymax=NULL, color=NULL) {
  require(ggplot2)
  attach(df)
  on.exit(detach(df))

  g1 <- ggplot(df, aes(x={{x}}, y={{y}}, color={{color}})) + geom_point(na.rm=TRUE)
  if (!is.null(x.sd)) {
    g1 <- g1 + geom_errorbarh(aes(y={{y}}, xmin={{x}}-{{x.sd}}, xmax={{x}}+{{x.sd}}, color={{color}}))
  }
  if (!is.null(xmin) && !is.null(xmax)) {
    g1 <- g1 + geom_errorbarh(aes(y={{y}}, xmin={{xmin}}, xmax={{xmax}}, color={{color}}))
  }
  if (!is.null(y.sd)) {
    g1 <- g1 + geom_errorbar(aes(x={{x}}, ymin={{y}}-{{y.sd}}, ymax={{y}}+{{y.sd}}, color={{color}}))
  }
  if (!is.null(ymin) && !is.null(ymax)) {
    g1 <- g1 + geom_errorbar(data=df, mapping=aes(x={{x}}, ymin={{ymin}}, ymax={{ymax}}, color={{color}}))
  }
  g1
}
