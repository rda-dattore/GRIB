# GRIB Edition 2 Parameter Code Tables

**Master tables** (maintained by WMO) have the format:
- WMO_GRIB2.\<V\>.xml

where \<V\> = \<M\>-0, with \<M\> being the master version number
- e.g. "WMO_GRIB2.18-0.xml" for master version 18

**Local tables** (maintained by individual data centers) have the format:
- WMO_GRIB2.\<C\>.\<V\>.xml

where \<C\> = \<C1\>-\<C2\>, with \<C1\> being the data center code (allocated by WMO) and \<C2\> being the sub-center code (allocated by the data center)

and \<V\> = \<M\>-\<L\>, with \<M\> being the master version number (allocated by WMO) and \<L\> being the local version number (allocated by the data center)
- e.g. "WMO_GRIB2.7-2.5-0.xml" for NCEP (data center 7) Ensemble Products (sub-center 2) using WMO master version 5 with local additions defined by NCEP's version 0
