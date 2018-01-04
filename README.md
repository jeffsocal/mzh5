# mzh5
##Wrapper to XCMS to Create a consistent R List Object or HDF5 file

Create a logical object to store and access LCMS data, either as an RDS object or HDF5 file. Conversion creates a rasterized lcms matrix defined by the desired fidelity, allowing (using HDF5 format) access to ms1-'hyperslabs' without loading the entire file. Furthermore, this format preserves the ability to collect profile data without increasing disk-space, and in many cases is more efficient then centroid data.

###Motivation

Current projects aimed at open access data for LCMS are either terribly space inefficient (xml), or a poor attempt at utilizing a innovations in big data (mz5). While there is a clear motivation to adhere to the HUPO Proteomics Standards for data integrity, there are only a few essential data constructs - Measurement Time, Measurement Conditions, and M/Z ~ Intensity. The attempt here is to improve that last construct.

####M/Z ~ Intensity Data Collection

Current data from every instrument manufacture records both m/z and intensity for every scan event - the former being the most redundant. Albeit, some experimenters call for different m/z ranges to be collected, however that shouldn't be the case for data recording. By standardizing the m/z axis of collection, for instance, to 0.01Da intervals for TOF/QQQ/QE proteomics and 0.001Da for FTMS/Orbi proteomics, the m/z axis needs only to be defined once and contains sufficient resolution for complex samples such as plasma. This establishes a _matrix_ of data, which can be further refined by collecting ms1 scan events at regular intervals. Storing this _matrix_ of data in the HDF5 format, is more efficient then the proposed mz5 format, for both read and write operations. In addition, taking advantage of the HDF5 concepts, data can be accessed in part (via an HDF5 hyperslab) without the need to load the entire file into memory. This becomes particularly enticing when dealing with large, multi-experiment datasets.  

## R list object  
```
object/
  |--scans/
    |--data/                  :'data.frame':	n obs. of  4 variables:
      |--scan_type            : num 
      |--scan_type_n          : int
      |--scan_id              : int 
      |--scan_rt_sec          : num 
  |--targets/
    |--data/                  :'data.frame':	n obs. of  10 variables:
      |--msn_n                : int 
      |--msn_id               : int
      |--pre_mz               : num
      |--pre_rt_sec           : num
      |--pre_int              : num
      |--pre_z                : int
      |--pre_ce               : num
      |--msn_int_max          : num
      |--msn_sn               : num
      |--msn_qual             : num
  |--spectra/
    |--ms1/
      |--lc                   : num
      |--mz                   : num
      |--int                  : 'matrix': int 
    |--ms2/
      |--mz                   : num
      |--int                  : 'matrix': int 
  |--ms1features: list()
```

## Disk space esitmations  

|FILE TYPE              | Space |
|-----------------------|------:|
|raw instrument vendor  | 1     |
|.mzML                  | 4x    |
|.mz5                   | 0.5x  |
|.mzh5                  | 0.25x |
|.mzh5.rds              | 0.25x |
