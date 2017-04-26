# GRIB Parameter Code Tables
various GRIB1 and GRIB2 parameter code tables

* table names have the format of:
     WMO_GRIB<E>.<C>-<S>.<P>.xml

        where E = GRIB edition (1 or 2)
              C = originating center code from WMO Common Code Table C-1
              S = sub-center identification
              P = parameter table version
                  for GRIB1, this is a single code
                  for GRIB2, this is <M>.<L>
                    where M = master version
                          L = local version
