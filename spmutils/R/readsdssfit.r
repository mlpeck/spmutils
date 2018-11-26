readsdssfit <- function(fname, dname="spectra", extcor=TRUE, rest=TRUE) {
  require(FITSio)
  
  metastrings <- c("PLUG_RA","PLUG_DEC","SPECOBJID","BESTOBJID","PLATE","MJD","FIBERID",
                   "CLASS", "SUBCLASS","Z","Z_ERR","ZWARNING","VDISP","VDISP_ERR",
                   "WAVEMIN","WAVEMAX")
  meta.out <- tolower(metastrings)
  meta.out[meta.out=="plug_ra"] <- "ra"
  meta.out[meta.out=="plug_dec"] <- "dec"
  
  hd <- readFITS(file.path(dname, fname), hdu=0)
  vals <- hd$col
  names(vals) <- hd$colNames
  vals <- vals[metastrings]
  names(vals) <- meta.out
  meta <- vals
  hd <- readFITS(file.path(dname, fname), hdu=1)
  vals <- hd$col
  names(vals) <- hd$colNames
  lambda <- 10^(vals[["loglam"]])
  flux <- vals[["flux"]]
  ivar <- vals[["ivar"]]
  ivar[ivar <= 0] <- NA
  flux[is.na(ivar)] <- NA
  if (extcor) {
    if (!exists("ebvmap", envir=.GlobalEnv)) {
      data(ebvmap, package="sfdmap", envir=.GlobalEnv)
    }
    lb <- radec2lb(meta$ra, meta$dec)
    ebv <- sfdmap::lookupebv(lb)
    meta$ebv <- ebv
    meta$l <- lb$l
    meta$b <- lb$b
    extc <- elaw(lambda, ebv)
    flux <- flux*extc
    ivar <- ivar/extc^2
  }
  if (rest) {
    lambda <- lambda/(1+meta$z)
  }
  list(meta=meta, lambda=lambda, flux=flux, ivar=ivar)
}
