# GRIB Edition 2 Parameter Code Tables

**Master tables** (_maintained by WMO_) have the format:
- WMO_GRIB2.\<V\>.xml

where _\<V\> = \<M\>-0_, with _\<M\>_ being the master version number
- e.g. "WMO_GRIB2.18-0.xml" for master version 18

**Local tables** (_maintained by individual data centers_) have the format:
- WMO_GRIB2.\<C\>.\<V\>.xml

where _\<C\> = \<C1\>-\<C2\>_, with _\<C1\>_ being the data center code (_allocated by WMO_) and _\<C2\>_ being the sub-center code (_allocated by the data center_)

and _\<V\> = \<M\>-\<L\>_, with _\<M\>_ being the master version number (_allocated by WMO_) and _\<L\>_ being the local version number (_allocated by the data center_)
- e.g. "WMO_GRIB2.7-2.5-0.xml" for NCEP (data center 7) Ensemble Products (sub-center 2) using WMO master version 5 with local additions defined by NCEP's version 0
