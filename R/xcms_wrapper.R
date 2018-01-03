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


library(xcms)

adf <- function(x) as.data.frame(x)

# Wrapper to XCMS | read MZML
#
# @param file_name path/fileName.mzml
# @keywords mzml
# @export
# @examples
# xcmsReadMzML

xcmsReadMzML <- function(file_name){
  dat <- xcmsRaw(file_name, includeMSn = T)
  return(dat)
}

# Wrapper to XCMS | get an ms1 spectrum
#
# @param obj readMZML object
# @param scan_n scan number 0:N
# @keywords ms1, spectrum
# @export
# @examples
# getMS1Spectrum

getMS1Spectrum <- function(obj, scan_n = 1){
  ms1spec <- adf(getScan(obj,scan_n))
  return(ms1spec)
}

# Wrapper to XCMS | get an ms2 spectrum
#
# @param obj readMZML object
# @param scan_n scan number 0:N
# @keywords ms2, spectrum
# @export
# @examples
# getMS2Spectrum

getMS2Spectrum <- function(obj, scan_n = 1){
  ms2spec <- adf(getMsnScan(obj,scan_n))
  return(ms2spec)
}

# Wrapper to XCMS | get an msN spectrum
#
# @param obj readMZML object
# @param level ms level 1,2
# @param scan_n scan number 0:N
# @keywords ms2, spectrum
# @export
# @examples
# getSpectrum

getSpectrum <- function(obj, level = 1, scan_n = 1){
  if(level == 1)
    spec <- getMS1Spectrum(obj, scan_n)

  if(level == 2)
    spec <- getMS2Spectrum(obj, scan_n)

  return(spec)
}

# Wrapper to XCMS | get ms1 features using XCMS CentWave
#
# @param obj readMZML object
# @keywords ms1, peaks
# @export
# @examples
# getMS1Peaks

getMS1Peaks <- function(obj){
  pks <- adf(findPeaks.centWave(obj))
  # pks <- pks[,c(1,4:10)]
  # pks[,2:4] <- round(pks[,2:4],2)
  # colnames(pks) <- c('pk_mz','pk_rt_sec', 'pk_rt_start', 'pk_rt_end',
  #                 'pk_int', 'pk_int_base', 'pk_int_max', 'pk_sn')
  return(pks)
}

# Wrapper to XCMS | get scan attribute data
#
# @param obj readMZML object
# @keywords scans
# @export
# @examples
# getMSScans

getMSScans <- function(obj){

  scn_ms1 <- data.frame(scan_type=1,
                        scan_type_n = 1:length(obj@scantime),
                        scan_id = obj@acquisitionNum,
                        scan_rt_sec = round(obj@scantime,2))

  scn_ms2 <- data.frame(scan_type=2,
                        scan_type_n = 1:length(obj@msnRt),
                        scan_id = obj@msnAcquisitionNum,
                        scan_rt_sec = round(obj@msnRt,2))

  scn <- rbind(scn_ms1, scn_ms2)
  scn <- scn[order(scn$scan_rt_sec),]
  row.names(scn) <- 1:dim(scn)[1]

  return(scn)
}

# Wrapper to XCMS | get target scan attribute data
#
# @param obj readMZML object
# @keywords scans
# @export
# @examples
# getMS2Targets

getMS2Targets <- function(obj){

  pre <- data.frame(msn_n = 1:length(obj@msnPrecursorMz),
                    msn_id = obj@msnAcquisitionNum,
                    pre_mz = obj@msnPrecursorMz,
                    pre_rt_sec = round(obj@msnRt,2),
                    pre_int = obj@msnPrecursorIntensity,
                    pre_z = obj@msnPrecursorCharge,
                    pre_ce = obj@msnCollisionEnergy,
                    msn_int_max = NA,
                    msn_sn = NA,
                    msn_qual = NA)
  return(pre)
}


