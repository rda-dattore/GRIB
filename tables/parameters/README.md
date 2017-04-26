# GRIB Parameter Code Tables
various GRIB1 and GRIB2 parameter code tables

Table names have the format of:
  WMO_GRIB\<E\>.\<C\>-\<S\>.\<P\>.xml

where:
- \<E\> = GRIB edition (1 or 2)
- \<C\> = originating center code from [WMO Common Code Table C-1](https://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/LatestVERSION/WMO306_vI2_CommonTable_en.pdf)
- \<S\> = sub-center identification
- \<P\> = parameter table version:

  for GRIB1, \<P\> is a single value
  
  for GRIB2, \<P\> is \<M\>-\<L\>, where \<M\> = master version and \<L\> = local version
