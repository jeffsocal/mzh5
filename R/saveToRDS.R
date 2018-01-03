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


#' Save standardize ms1 data to an RDS file
#'
#' @param obj mzh5::readMzML object
#' @param file_name The path/filename to save to.
#' @keywords save, HDF5
#' @export
#' @examples
#' saveToRDS

saveToRDS <- function(obj, file_name){

  saveRDS(obj, file = file_name)

}
