################################################################################
# Copyright (2016) SoCal Bioinformatics Inc. All rights reserved.
# This script is the confidential and proprietary product of SoCal
# Bioinformatics Inc. Any unauthorized reproduction or transfer of
# the contents herein is strictly prohibited.
#
################################################################################
# AUTH:     Jeff Jones | SoCal Bioinofrmatics Inc (SCBI)
# DATE:     2018.01.01
# OWNER:    SoCal Bioinofrmatics Inc
# PROJECT:  SCBI | hdf5
# DESC:     Estimate the msn precursor ion intenisty of the isolated mz value
################################################################################


#' Estimate the msn precursor ion intenisty of the isolated mz value
#'
#' @param obj mzh5::readMzML object
#' @param verbose Print progress. Default T
#' @keywords precursor intensity
#' @export
#' @examples
#' dat <- readMzML("example.mzML")
#' dat <- estPreInt(dat)

estPreInt <- function(dat,
                      verbose = T ){


  msn <- dat$targets$data
  mz <- dat$spectra$ms1$mz
  lc <- dat$spectra$ms1$lc
  int <- dat$spectra$ms1$int

  if(verbose == T) {
    cat('Estimate precursor intensity\n')
    pb <- txtProgressBar(max=dim(msn)[1], style=3, width=20)
  }

  for ( i in 1:dim(msn)[1] ){
    m <- msn[i,]$pre_mz
    l <- msn[i,]$pre_rt_sec

    t_lc <- which(abs(lc - l) <= 5)
    t_mz <- which(abs(mz - m) <= 0.1)

    t_int <- int[t_lc,t_mz]

    msn$pre_int[i] <- max(t_int)

    setTxtProgressBar(pb, i)
  }

  dat$targets$data$pre_int <- msn$pre_int
  cat("\n")
  return(dat)
}
