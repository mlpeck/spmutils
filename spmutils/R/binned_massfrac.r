binned_massfrac <- function(sfh_post, ages=ages, ages.bins=c(0.1, 1.25, 2.25)) {
  dims <- dim(sfh_post)
  ind <- findInterval(ages.bins+0.01, 10^(ages-9))
  ind <- c(0, ind, length(ages))
  mfrac <- array(0, dim=c(dims[1], length(ind)-1, dims[3]))
  tmass <- apply(sfh_post, c(1,3), sum, na.rm=TRUE)
  tmass[tmass==0] <- NA
  for (i in 1:(length(ind)-1)) {
    massi <- apply(sfh_post[ ,(ind[i]+1):ind[i+1], ], c(1,3), sum, na.rm=TRUE)
    mfrac[,i,] <- massi/tmass
  }
  mfrac
}

sum_binned_massfrac <- function(mfrac, ages=ages, ages.bins=c(0.1, 1.25, 2.25)) {
  sum.mfrac <- data.frame(cbind(t(apply(mfrac, c(2,3), mean)), t(apply(mfrac, c(2,3), sd))))
  tnames <- paste("T", formatC(c(ages.bins,10^(ages[length(ages)]-9)), format="f", digits=2), sep="_")
  names(sum.mfrac) <- c(paste(tnames, "m", sep="_"), paste(tnames, "sd", sep="_"))
  sum.mfrac
}
