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
# DESC:     create an R list object from xcms mzml data
#           rasterize the ms1 lcms to "standardize the matrix"
################################################################################


#' Standardize ms1 data to consistent mz/lc values
#'
#' @param file_name path/fileName.mzml
#' @param mz_fidelity_amu The sampling fidelity for mz axis. Default 0.01.
#' @param lc_fidelity_sec The sampling fidelity for lc axis. Default 1.
#' @param get_ms1features Populate the object with XCMS fetaures. Default F.
#' @param verbose Print progress. Default T
#' @keywords lcms, rasterize
#' @export
#' @examples
#' dat <- readMzML("example.mzML")
#' str(dat)
#' List of 4
#' $ scans      :List of 1
#' ..$ data:'data.frame':	11075 obs. of  4 variables:
#' .. ..$ scan_type  : num [1:11075] 1 1 1 1 1 1 1 1 1 1 ...
#' .. ..$ scan_type_n: int [1:11075] 1 2 3 4 5 6 7 8 9 10 ...
#' .. ..$ scan_id    : int [1:11075] 1 2 3 4 5 6 7 8 9 10 ...
#' .. ..$ scan_rt_sec: num [1:11075] 0.35 0.76 1.17 1.59 2 2.41 2.82 3.24 3.66 4.08 ...
#' $ targets    :List of 1
#' ..$ data:'data.frame':	4773 obs. of  10 variables:
#' .. ..$ msn_n      : int [1:4773] 1 2 3 4 5 6 7 8 9 10 ...
#' .. ..$ msn_id     : int [1:4773] 231 264 265 266 267 289 295 307 310 313 ...
#' .. ..$ pre_mz     : num [1:4773] 470 391 397 436 559 ...
#' .. ..$ pre_rt_sec : num [1:4773] 102 117 117 117 118 ...
#' .. ..$ pre_int    : num [1:4773] 0 0 0 0 0 0 0 0 0 0 ...
#' .. ..$ pre_z      : int [1:4773] 0 0 0 0 0 0 0 0 0 0 ...
#' .. ..$ pre_ce     : num [1:4773] 0 0 0 0 0 0 0 0 0 0 ...
#' .. ..$ msn_int_max: logi [1:4773] NA NA NA NA NA NA ...
#' .. ..$ msn_sn     : logi [1:4773] NA NA NA NA NA NA ...
#' .. ..$ msn_qual   : logi [1:4773] NA NA NA NA NA NA ...
#' $ spectra    :List of 2
#' ..$ ms1:List of 3
#' .. ..$ lc : num [1:7000] 0 0.514 1.029 1.543 2.057 ...
#' .. ..$ mz : num [1:90000] 350 350 350 350 350 ...
#' .. ..$ int: int [1:7000, 1:90000] 0 0 0 0 0 0 0 0 16 16 ...
#' ..$ ms2:List of 2
#' .. ..$ mz : num [1:86000] 70 70 70 70 70 ...
#' .. ..$ int: int [1:4773, 1:86000] NA NA NA NA NA NA NA NA NA NA ...
#' $ ms1features: list()

readMzML <- function(file_name,
                          mz_fidelity_amu = 0.01,
                          lc_fidelity_sec = 1,
                          get_ms1features = F,
                          verbose = T ){

  if(verbose == T) cat('Read MZML file ... ')

  obj <- xcmsReadMzML(file_name)

  if(verbose == T) cat(' done\n')


  dat_scans <- getMSScans(obj)
  dat_targets <- getMS2Targets(obj)
  dat_features <- c()

  if(get_ms1features == T)
    dat_features <- getMS1Peaks(mzml)


  ls <- list()
  ls[['scans']] <- list()
  ls[['targets']] <- list()
  ls[['spectra']] <- list()
  ls[['ms1features']] <- list()

  ls[['scans']][['data']] <- dat_scans
  ls[['targets']][['data']] <- dat_targets
  ls[['ms1features']][['data']] <- dat_features

  for(ms_level in 1:2){
    v_scans <- dat_scans[dat_scans$scan_type == ms_level,]$scan_type_n

    if(verbose == T)
      cat('MS level ', ms_level, " rasterization\n")

    spec_mz <- c()
    spec_int <- list()

    # determine range
    ospec <- getSpectrum(obj, ms_level, 1)
    ms_range <- c(floor(min(ospec$mz)/10)*10,
                  ceiling(max(ospec$mz)/10)*10)

    n_rows <- length(v_scans)

    if(verbose == T){
      cat(' mz\n')
      pb <- txtProgressBar(max=n_rows, style=3, width=20)
    }

    for(i in 1:n_rows){
      scan_n <- v_scans[i]
      ospec <- getSpectrum(obj, ms_level, scan_n)
      spec <- ksmooth(ospec$mz,
                      ospec$intensity,
                      'normal',
                      n.points=diff(ms_range)*100,
                      range.x = ms_range,
                      bandwidth = mz_fidelity_amu)
      spec_mz <- spec$x
      spec_int[[scan_n]] <- as.integer(ceiling(spec$y))

      if(verbose == T)
        setTxtProgressBar(pb, i)
    }
    if(verbose == T) close(pb)

    spec_mtrx <- matrix(unlist(spec_int), byrow=TRUE, nrow=length(spec_int) )

    ls[['spectra']][[paste0('ms',ms_level)]] <- list()

    # only rasterize LC if ms_level == 1
    if(ms_level == 1){
      ms1_rt <- dat_scans[dat_scans$scan_type == 1,]$scan_rt_sec
      n_points <- ceiling(length(ms1_rt)/10^floor(log10(length(ms1_rt))))*10^floor(log10(length(ms1_rt)))
      n_rows <- dim(spec_mtrx)[2]
      spec_int <- list()

      if(verbose == T){
        cat(' lc\n')
        pb <- txtProgressBar(max=n_rows, style=3, width=20)
      }

      for(i in 1:n_rows){
        spec <- ksmooth(ms1_rt,
                        spec_mtrx[,i],
                        'normal',
                        n.points= n_points,
                        range.x = c(0,max(floor(ms1_rt))),
                        bandwidth = lc_fidelity_sec)
        spec_lc <- spec$x
        this_spec_int <- as.numeric(spec$y)
        this_spec_int[which(is.na(this_spec_int))] <- 0
        spec_int[[i]] <- as.integer(floor(this_spec_int))

        if(verbose == T)
          setTxtProgressBar(pb, i)
      }
      if(verbose == T) close(pb)

      spec_mtrx <- t(matrix(unlist(spec_int), byrow=TRUE, nrow=n_rows ))

      ls[['spectra']][[paste0('ms',ms_level)]][['lc']] <- spec_lc
    }

    ls[['spectra']][[paste0('ms',ms_level)]][['mz']] <- spec_mz
    ls[['spectra']][[paste0('ms',ms_level)]][['int']] <- spec_mtrx
  }

  return(ls)
}
