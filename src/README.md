# src directory

(documentation for each code can be found at the top of each file)

- unpackgrib1.c
  - C code for decoding GRIB1 messages
  
- unpackgrib2.c
  - C code for decoding GRIB2 messages
  
- grib1to2.c
  - C program for converting from GRIB1 to GRIB2 (also requires unpackgrib1.c)
  
- grib2to1.c
  - C program for converting from GRIB2 to GRIB1 (also requires unpackgrib2.c)

- grib2_read_example.c
  - sample C program to read a GRIB2 file
