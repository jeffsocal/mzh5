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

library(h5)


#' Save standardize ms1 data to an HDF5 file
#'
#' @param obj mzh5::readMzML object
#' @param file_name The path/filename to save to.
#' @keywords save, HDF5
#' @export
#' @examples
#' saveToHDF5

saveToHDF5 <- function(obj, file_name){

  # ls[['scans']] <- list()
  # ls[['targets']] <- list()
  # ls[['spectra']] <- list()
  # ls[['ms1features']] <- list()
  #
  # ls[['scans']][['data']] <- dat_scans
  # ls[['targets']][['data']] <- dat_targets
  # ls[['ms1features']][['data']] <- dat_features

  new_file <- h5file(file_name)

  new_file['scans/data'] <- as.matrix(obj$dat_scans)
  new_file['scans/header'] <- colnames(obj$dat_scans)

  new_file['targets/data'] <- as.matrix(obj$dat_targets)
  new_file['targets/header'] <- colnames(obj$dat_targets)

  new_file['ms1features/data'] <- as.matrix(obj$dat_features)
  new_file['ms1features/header'] <- colnames(obj$dat_features)

  for(ms_level in 1:2){

    if( !paste0('ms',ms_level) %in% obj[['spectra']] )
      next()

    for(ms_part in c('lc','mz','int')){

      if( !ms_part %in% obj[['spectra']][[paste0('ms',ms_level)]] )
        next()

      new_file[paste0('spectra/raw/ms',ms_level,'/',ms_part)] <-
        obj[['spectra']][[paste0('ms',ms_level)]][[ms_part]]
    }
  }

  h5close(new_file)

}

