# mzh5
Wrapper to XCMS to Create a consistent R List Object or HDF5 file

Create a logical object to store and access LCMS data, either as an RDS object or HDF5 file. Conversion creates a rasterized lcms matrix defined by the desired fidelity, allowing (using HDF5 format) access to 'hyperslabs' without loading the entire file. Furthermore, this format preserves the ability to collect profile data without increasing disk-space.

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
