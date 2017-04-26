/*
** File: grib2to1.c
**
** Author:  Bob Dattore
**          NCAR/DSS
**          dattore@ucar.edu
**          (303) 497-1825
**
** Purpose: to provide a simple C program for converting from GRIB2 to GRIB1
**
** Revision History:
**          20 Feb 2008 - first version; will convert analyses and forecasts
**            on a regular latitude/longitude grid
**          31 Mar 2008 - added code for Lambert conformal grids (GDS template
**            3.30), statistically processed data (PDS template 4.8), and some
**            NCEP-specific parameter mappings
**          03 Apr 2008 - bug fixes
**          17 Apr 2008 - change reference value packing so that it works on
**            little-endian machines
**          16 Jul 2008 - added more GRIB2 to GRIB1 parameter mappings, added
**            PDS packing for GRIB2 Product Definition Template 4.2
**          22 Jul 2008 - patch for NCEP max/min temperature grids
**          03 May 2010 - patch for NCEP CFSR monthly grids
**          16 Jul 2013 - bug fix in ieee2ibm: sizeof should check the size of
**                        an int, not a size_t
**          29 Nov 2013 - added GRIB2 to GRIB1 parameter mappings for some UK
**            Met Office products, added PDS packing for GRIB2 Product Defintion
**            Template 4.15
**          24 Sep 2014 - when setting resolution and component flags: added
**            parentheses around bitwise '&' operations as they have less
**            precedence than '=='
**          14 Aug 2015 - added a few more NCEP GRIB2->GRIB1 mappings
**          25 Feb 2016 - updated 'mapTimeRange' to set the time range indicator
**            to '10' when the time units are 'minutes' and the time range
**            indicator would otherwise be '0'
**          17 Apr 2017 - added a few more NCEP GRIB2->GRIB1 mappings, changed
**            parameter code mapping to 255 instead of exiting when a code can't
**            be mapped
**          21 Apr 2017 - changed return value of mapParameterData to the
**            ParameterData structure because some parameter codes belong to
**            parameter tables other than version 3
**
** Contact Bob Dattore at dattore@ucar.edu to get conversions for other products
** and grid definitions added.
**
** You will need to download the GRIB2 decoder:
**    https://raw.githubusercontent.com/rda-dattore/GRIB/master/src/unpackgrib2.c
** It must be in the same directory as this program.
**
** If you want to decode jpeg-compressed grids (widely used by NCEP), you will
**   also need to link with the JasPer library, which needs the JPEG-6b library.
**     You can get the JasPer code at
**       http://www.ece.uvic.ca/~mdadams/jasper/
**     You can get the installation code for libjpeg.a at
**       http://www.ijg.org/files/jpegsrc.v6b.tar.gz
**
** Example compile command:
**    % cc -std=c99 -o grib2to1 grib2to1.c
**
** Example compile command with JasPer support:
**    % cc -std=c99 -DJASPER -o grib2to1 -I<jasper include path> -L<jasper library path> grib2to1.c -ljasper
**      where <jasper include path> is the directory path of the jasper header
**              files
**            <jasper library path> is the directory path of the jasper library
**
** If the compiler complains about the "pow" function being an undefined symbol,
**   include the math library in the compile
**      e.g. % cc -std=c99 -o grib2to1 grib2to1.c -lm
**
** To use the program:
**    % grib2to1 <name of GRIB2 file to convert> <name of GRIB1 file to create>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "unpackgrib2.c"

/* setBits sets the contents of the various GRIB octets
**   buf is the GRIB buffer as a stream of bytes
**   src is the value of the octet(s) to set
**   off is the offset in BITS from the beginning of the buffer to the beginning
**       of the octet(s) to be packed
**   bits is the number of BITS to pack - will be a multiple of 8 since GRIB
**       octets are 8 bits long
*/
void setBits(unsigned char *buf,int src,size_t off,size_t bits)
{
// no work to do
  if (bits == 0) {
    return;
  }
  size_t src_size=sizeof(int)*8;
  if (bits > src_size) {
    fprintf(stderr,"Error: packing %d bits from a %d-bit field\n",bits,src_size);
    exit(1);
  }
  else {
// create masks to use when right-shifting (necessary because different
// compilers do different things when right-shifting a signed bit-field)
    unsigned char bmask=1;
    size_t buf_size=sizeof(unsigned char)*8;
    for (size_t n=1; n < buf_size; ++n) {
      bmask<<=1;
      ++bmask;
    }
    int smask=1;
    for (size_t n=1; n < src_size; ++n) {
	smask<<=1;
	++smask;
    }
// get number of words and bits to skip before packing begins
    size_t wskip=off/buf_size;
    size_t bskip=off % buf_size;
    size_t lclear=bskip+bits;
    size_t rclear=buf_size-bskip;
    unsigned char left= (rclear != buf_size) ? (buf[wskip]&(bmask<<rclear)) : 0;
    if (lclear <= buf_size) {
// all bits to be packed are in the current word; clear the field to be
// packed
	unsigned char right= (lclear != buf_size) ? (buf[wskip]&~(bmask<<(buf_size-lclear))) : 0;
// fill the field to be packed
	buf[wskip]= (src_size != bits) ? src&~(smask<<bits) : src;
	buf[wskip]=left|right|(buf[wskip]<<(rclear-bits));
    }
    else {
// bits to be packed cross a byte boundary(ies); clear the bit field to be
// packed
	size_t more=bits-rclear;
	buf[wskip]=left|((src>>more)&~(smask<<(bits-more)));
// clear the next (or part of the next) word and pack those bits
	while (more > buf_size) {
	  more-=buf_size;
	  buf[++wskip]=(src>>more)&~(smask<<(src_size-more));
	}
	++wskip;
	more=buf_size-more;
	unsigned char right= (more != buf_size) ? (buf[wskip]&~(bmask<<more)) : 0;
	buf[wskip]= (buf_size > src_size) ? src&~(bmask<<src_size) : src;
	buf[wskip]=right|(buf[wskip]<<more);
    }
  }
}

typedef struct {
  int table_version;
  int parameter_code;
} ParameterData;

ParameterData mapParameterData(GRIB2Message *msg,int grid_number)
{
  ParameterData pdata;
  switch (msg->disc) {
// meteorological products
    case 0:
    {
	switch (msg->grids[grid_number].md.param_cat) {
// temperature parameters
	  case 0:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=11;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=12;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=13;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=14;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=15;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=16;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=17;
		  return pdata;
		}
		case 7:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=18;
		  return pdata;
		}
		case 8:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=19;
		  return pdata;
		}
		case 9:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=25;
		  return pdata;
		}
		case 10:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=121;
		  return pdata;
		}
		case 11:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=122;
		  return pdata;
		}
		case 21:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=131;
			pdata.parameter_code=193;
			return pdata;
		    }
		  }
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=229;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// moisture parameters
	  case 1:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=51;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=52;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=53;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=54;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=55;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=56;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=57;
		  return pdata;
		}
		case 7:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=59;
		  return pdata;
		}
		case 8:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=61;
		  return pdata;
		}
		case 9:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=62;
		  return pdata;
		}
		case 10:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=63;
		  return pdata;
		}
		case 11:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=66;
		  return pdata;
		}
		case 12:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=64;
		  return pdata;
		}
		case 13:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=65;
		  return pdata;
		}
		case 14:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=78;
		  return pdata;
		}
		case 15:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=79;
		  return pdata;
		}
		case 16:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=99;
		  return pdata;
		}
		case 22:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=153;
			return pdata;
		    }
		  }
		}
		case 39:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=194;
			return pdata;
		    }
		  }
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
			pdata.table_version=3;
			pdata.parameter_code=140;
			return pdata;
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=141;
			return pdata;
		    }
		  }
		}
		case 194:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=142;
			return pdata;
		    }
		  }
		}
		case 195:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=143;
			return pdata;
		    }
		  }
		}
		case 196:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=214;
			return pdata;
		    }
		  }
		}
		case 197:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=135;
			return pdata;
		    }
		  }
		}
		case 199:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=228;
			return pdata;
		    }
		  }
		}
		case 200:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=145;
			return pdata;
		    }
		  }
		}
		case 201:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=238;
			return pdata;
		    }
		  }
		}
		case 206:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=186;
			return pdata;
		    }
		  }
		}
		case 207:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=198;
			return pdata;
		    }
		  }
		}
		case 208:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=239;
			return pdata;
		    }
		  }
		}
		case 213:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=243;
			return pdata;
		    }
		  }
		}
		case 214:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=245;
			return pdata;
		    }
		  }
		}
		case 215:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=249;
			return pdata;
		    }
		  }
		}
		case 216:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=159;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// momentum parameters
	  case 2:
	  {
	    switch(msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=31;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=32;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=33;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=34;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=35;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=36;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=37;
		  return pdata;
		}
		case 7:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=38;
		  return pdata;
		}
		case 8:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=39;
		  return pdata;
		}
		case 9:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=40;
		  return pdata;
		}
		case 10:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=41;
		  return pdata;
		}
		case 11:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=42;
		  return pdata;
		}
		case 12:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=43;
		  return pdata;
		}
		case 13:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=44;
		  return pdata;
		}
		case 14:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=4;
		  return pdata;
		}
		case 15:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=45;
		  return pdata;
		}
		case 16:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=46;
		  return pdata;
		}
		case 17:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=124;
		  return pdata;
		}
		case 18:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=125;
		  return pdata;
		}
		case 19:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=126;
		  return pdata;
		}
		case 20:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=123;
		  return pdata;
		}
		case 22:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=180;
			return pdata;
		    }
		  }
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=136;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=172;
			return pdata;
		    }
		  }
		}
		case 194:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=196;
			return pdata;
		    }
		  }
		}
		case 195:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=197;
			return pdata;
		    }
		  }
		}
		case 196:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=252;
			return pdata;
		    }
		  }
		}
		case 197:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=253;
			return pdata;
		    }
		  }
		}
		case 224:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=129;
			pdata.parameter_code=241;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// mass parameters
	  case 3:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=1;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=2;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=3;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=5;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=6;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=7;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=8;
		  return pdata;
		}
		case 7:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=9;
		  return pdata;
		}
		case 8:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=26;
		  return pdata;
		}
		case 9:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=27;
		  return pdata;
		}
		case 10:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=89;
		  return pdata;
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=130;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=222;
			return pdata;
		    }
		  }
		}
		case 194:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=147;
			return pdata;
		    }
		  }
		}
		case 195:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=148;
			return pdata;
		    }
		  }
		}
		case 196:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=221;
			return pdata;
		    }
		  }
		}
		case 197:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=230;
			return pdata;
		    }
		  }
		}
		case 198:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=129;
			return pdata;
		    }
		  }
		}
		case 199:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=137;
			return pdata;
		    }
		  }
		}
		case 200:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=129;
			pdata.parameter_code=141;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// short-wave radiation parameters
	  case 4:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=111;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=113;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=116;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=117;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=118;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=119;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=120;
		  return pdata;
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=204;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=211;
			return pdata;
		    }
		  }
		}
		case 196:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=161;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// long-wave radiation parameters
	  case 5:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=112;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=114;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=115;
		  return pdata;
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=205;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=212;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// cloud parameters
	  case 6:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=58;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=71;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=72;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=73;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=74;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=75;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=76;
		  return pdata;
		}
		case 25:
		{
		  switch (msg->center_id) {
		    case 74:
		    {
			pdata.table_version=140;
			pdata.parameter_code=174;
			return pdata;
		    }
		  }
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=213;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=146;
			return pdata;
		    }
		  }
		}
		case 201:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=133;
			pdata.parameter_code=191;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// thermodynamic stability index parameters
	  case 7:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=24;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=77;
		  return pdata;
		}
		case 6:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=157;
			return pdata;
		    }
		  }
		}
		case 7:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=156;
			return pdata;
		    }
		  }
		}
		case 8:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=190;
			return pdata;
		    }
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=131;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=132;
			return pdata;
		    }
		  }
		}
		case 194:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=254;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// aerosol parameters
	  case 13:
	  {
	    break;
	  }
// trace gas parameters
	  case 14:
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=10;
		  return pdata;
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=154;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
// radar parameters
	  case 15:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=21;
		  return pdata;
		}
		case 7:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=22;
		  return pdata;
		}
		case 8:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=23;
		  return pdata;
		}
	    }
	    break;
	  }
// forecast radar imagery parameters
	  case 16:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 195:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=129;
			pdata.parameter_code=211;
			return pdata;
		    }
		  }
		}
		case 196:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=129;
			pdata.parameter_code=212;
			return pdata;
		    }
		  }
		}
	    }
	  }
// nuclear/radiology parameters
	  case 18:
	  {
	    break;
	  }
// physical atmospheric property parameters
	  case 19:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=20;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=84;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=60;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=67;
		  return pdata;
		}
		case 20:
		{
		  switch (msg->center_id) {
		    case 74:
		    {
			switch (msg->md.spatial_proc.type) {
			  case 0:
			  {
			    pdata.table_version=3;
			    pdata.parameter_code=168;
			    return pdata;
			  }
			  case 2:
			  {
			    pdata.table_version=3;
			    pdata.parameter_code=169;
			    return pdata;
			  }
			}
		    }
		  }
		}
		case 21:
		{
		  switch (msg->center_id) {
		    case 74:
		    {
			switch (msg->md.spatial_proc.type) {
			  case 0:
			  {
			    pdata.table_version=3;
			    pdata.parameter_code=170;
			    return pdata;
			  }
			  case 2:
			  {
			    pdata.table_version=3;
			    pdata.parameter_code=171;
			    return pdata;
			  }
			}
		    }
		  }
		}
		case 22:
		{
		  switch (msg->center_id) {
		    case 74:
		    {
			switch (msg->md.spatial_proc.type) {
			  case 0:
			  {
			    pdata.table_version=3;
			    pdata.parameter_code=172;
			    return pdata;
			  }
			  case 2:
			  {
			    pdata.table_version=3;
			    pdata.parameter_code=173;
			    return pdata;
			  }
			}
		    }
		  }
		}
		case 204:
		{
		  switch (msg->center_id) {
		    case 7:
			pdata.table_version=3;
			pdata.parameter_code=209;
			return pdata;
		  }
		}
	    }
	    break;
	  }
	}
	break;
    }
// hydrologic products
    case 1:
    {
	switch (msg->grids[grid_number].md.param_cat) {
// hydrology basic products
	  case 0:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=234;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=235;
			return pdata;
		    }
		  }
		}
	    }
	  }
	  case 1:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=195;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=194;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
	}
	break;
    }
// land surface products
    case 2:
    {
	switch (msg->grids[grid_number].md.param_cat) {
// vegetation/biomass
	  case 0:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=81;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=83;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=85;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=86;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=87;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=90;
		  return pdata;
		}
		case 192:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=144;
			return pdata;
		    }
		  }
		}
		case 193:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=155;
			return pdata;
		    }
		  }
		}
		case 194:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=207;
			return pdata;
		    }
		  }
		}
		case 195:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=208;
			return pdata;
		    }
		  }
		}
		case 196:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=223;
			return pdata;
		    }
		  }
		}
		case 197:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=226;
			return pdata;
		    }
		  }
		}
		case 198:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=225;
			return pdata;
		    }
		  }
		}
		case 201:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=130;
			pdata.parameter_code=219;
			return pdata;
		    }
		  }
		}
		case 207:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=3;
			pdata.parameter_code=201;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
	  case 3:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 203:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=130;
			pdata.parameter_code=220;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
	  case 4:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 2:
		{
		  switch (msg->center_id) {
		    case 7:
		    {
			pdata.table_version=129;
			pdata.parameter_code=250;
			return pdata;
		    }
		  }
		}
	    }
	    break;
	  }
	}
	break;
    }
// oceanographic products
    case 10:
    {
	switch (msg->grids[grid_number].md.param_cat) {
// waves parameters
	  case 0:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=28;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=29;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=30;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=100;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=101;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=102;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=103;
		  return pdata;
		}
		case 7:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=104;
		  return pdata;
		}
		case 8:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=105;
		  return pdata;
		}
		case 9:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=106;
		  return pdata;
		}
		case 10:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=107;
		  return pdata;
		}
		case 11:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=108;
		  return pdata;
		}
		case 12:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=109;
		  return pdata;
		}
		case 13:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=110;
		  return pdata;
		}
	    }
	    break;
	  }
// currents parameters
	  case 1:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=47;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=48;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=49;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=50;
		  return pdata;
		}
	    }
	    break;
	  }
// ice parameters
	  case 2:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=91;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=92;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=93;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=94;
		  return pdata;
		}
		case 4:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=95;
		  return pdata;
		}
		case 5:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=96;
		  return pdata;
		}
		case 6:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=97;
		  return pdata;
		}
		case 7:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=98;
		  return pdata;
		}
	    }
	    break;
	  }
// surface properties parameters
	  case 3:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=80;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=82;
		  return pdata;
		}
	    }
	    break;
	  }
// sub-surface properties parameters
	  case 4:
	  {
	    switch (msg->grids[grid_number].md.param_num) {
		case 0:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=69;
		  return pdata;
		}
		case 1:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=70;
		  return pdata;
		}
		case 2:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=68;
		  return pdata;
		}
		case 3:
		{
		  pdata.table_version=3;
		  pdata.parameter_code=88;
		  return pdata;
		}
	    }
	    break;
	  }
	}
	break;
    }
  }
  fprintf(stderr,"Warning: no GRIB1 parameter code for discipline %d, parameter category %d, parameter number %d, center %d; setting to 255\n",msg->disc,msg->grids[grid_number].md.param_cat,msg->grids[grid_number].md.param_num,msg->center_id);
  pdata.table_version=3;
  pdata.parameter_code=255;
  return pdata;
}

int mapLevelData(GRIB2Grid *grid,int *level_type,int *level1,int *level2,int center)
{
  if (grid->md.lvl2_type != 255 && grid->md.lvl1_type != grid->md.lvl2_type) {
    fprintf(stderr,"Unable to indicate a layer bounded by different level types %d and %d in GRIB1\n",grid->md.lvl1_type,grid->md.lvl2_type);
    exit(1);
  }
  *level1=*level2=0;
  switch (grid->md.lvl1_type) {
    case 1:
	*level_type=1;
	break;
    case 2:
	*level_type=2;
	break;
    case 3:
	*level_type=3;
	break;
    case 4:
	*level_type=4;
	break;
    case 5:
	*level_type=5;
	break;
    case 6:
	*level_type=6;
	break;
    case 7:
	*level_type=7;
	break;
    case 8:
	*level_type=8;
	break;
    case 9:
	*level_type=9;
	break;
    case 20:
	*level_type=20;
	break;
    case 100:
	if (grid->md.lvl2_type == 255) {
	  *level_type=100;
	  *level1=grid->md.lvl1/100.;
	}
	else {
	  *level_type=101;
	  *level1=grid->md.lvl1/1000.;
	  *level2=grid->md.lvl2/1000.;
	}
	break;
    case 101:
	*level_type=102;
	break;
    case 102:
	if (grid->md.lvl2_type == 255) {
	  *level_type=103;
	  *level1=grid->md.lvl1;
	}
	else {
	  *level_type=104;
	  *level1=grid->md.lvl1/100.;
	  *level2=grid->md.lvl2/100.;
	}
	break;
    case 103:
	if (grid->md.lvl2_type == 255) {
	  *level_type=105;
	  *level1=grid->md.lvl1;
	}
	else {
	  *level_type=106;
	  *level1=grid->md.lvl1/100.;
	  *level2=grid->md.lvl2/100.;
	}
	break;
    case 104:
	if (grid->md.lvl2_type == 255) {
	  *level_type=107;
	  *level1=grid->md.lvl1*10000.;
	}
	else {
	  *level_type=108;
	  *level1=grid->md.lvl1*100.;
	  *level2=grid->md.lvl2*100.;
	}
	break;
    case 105:
	*level1=grid->md.lvl1;
	if (grid->md.lvl2_type == 255)
	  *level_type=109;
	else {
	  *level_type=110;
	  *level2=grid->md.lvl2;
	}
	break;
    case 106:
	*level1=grid->md.lvl1*100.;
	if (grid->md.lvl2_type == 255)
	  *level_type=111;
	else {
	  *level_type=112;
	  *level2=grid->md.lvl2*100.;
	}
	break;
    case 107:
	if (grid->md.lvl2_type == 255) {
	  *level_type=113;
	  *level1=grid->md.lvl1;
	}
	else {
	  *level_type=114;
	  *level1=475.-grid->md.lvl1;
	  *level2=475.-grid->md.lvl2;
	}
	break;
    case 108:
	*level1=grid->md.lvl1/100.;
	if (grid->md.lvl2_type == 255)
	  *level_type=115;
	else {
	  *level_type=116;
	  *level2=grid->md.lvl2/100.;
	}
	break;
    case 109:
	*level_type=117;
	*level1=grid->md.lvl1*1000000000.;
	break;
    case 111:
	if (grid->md.lvl2_type == 255) {
	  *level_type=119;
	  *level1=grid->md.lvl1*10000.;
	}
	else {
	  *level_type=120;
	  *level1=grid->md.lvl1*100.;
	  *level2=grid->md.lvl2*100.;
	}
	break;
    case 117:
	fprintf(stderr,"There is no GRIB1 level code for 'Mixed layer depth'\n");
	exit(1);
    case 160:
	*level_type=160;
	*level1=grid->md.lvl1;
	break;
    case 200:
	switch (center) {
	  case 7:
	    *level_type=200;
	    break;
	}
	break;
  }
}

int mapStatisticalEndTime(GRIB2Message *msg,GRIB2Grid *grid)
{
  switch (grid->md.time_unit) {
    case 0:
	return (grid->md.stat_proc.etime/100 % 100)-(msg->time/100 % 100);
    case 1:
	 return (grid->md.stat_proc.etime/10000-msg->time/10000);
    case 2:
	return (grid->md.stat_proc.edy-msg->dy);
    case 3:
	return (grid->md.stat_proc.emo-msg->mo);
    case 4:
	return (grid->md.stat_proc.eyr-msg->yr);
    default:
	fprintf(stderr,"Unable to map end time with units %d to GRIB1\n",grid->md.time_unit);
	exit(1);
  }
}

void mapTimeRange(GRIB2Message *msg,GRIB2Grid *grid,int *p1,int *p2,int *t_range,int *n_avg,int *n_missing,int center)
{
  size_t n;

  switch (grid->md.pds_templ_num) {
    case 0:
    case 1:
    case 2:
    case 15:
	if (grid->md.time_unit == 0) {
	  *t_range=10;
	}
	else {
	  *t_range=0;
	}
	*p1=grid->md.fcst_time;
	*p2=0;
	*n_avg=*n_missing=0;
	break;
    case 8:
    case 11:
    case 12:
	if (grid->md.stat_proc.num_ranges > 1) {
	  if (center == 7 && grid->md.stat_proc.num_ranges == 2) {
/* NCEP CFSR monthly grids */
	    *p2=grid->md.stat_proc.incr_length[0];
	    *p1=*p2-grid->md.stat_proc.time_length[1];
	    *n_avg=grid->md.stat_proc.time_length[0];
	    switch (grid->md.stat_proc.proc_code[0]) {
		case 193:
		  *t_range=113;
		  break;
		case 194:
		  *t_range=123;
		  break;
		case 195:
		  *t_range=128;
		  break;
		case 196:
		  *t_range=129;
		  break;
		case 197:
		  *t_range=130;
		  break;
		case 198:
		  *t_range=131;
		  break;
		case 199:
		  *t_range=132;
		  break;
		case 200:
		  *t_range=133;
		  break;
		case 201:
		  *t_range=134;
		  break;
		case 202:
		  *t_range=135;
		  break;
		case 203:
		  *t_range=136;
		  break;
		case 204:
		  *t_range=137;
		  break;
		case 205:
		  *t_range=138;
		  break;
		case 206:
		  *t_range=139;
		  break;
		case 207:
		  *t_range=140;
		  break;
		default:
		  fprintf(stderr,"Unable to map NCEP statistical process code %d to GRIB1\n",grid->md.stat_proc.proc_code[0]);
		  exit(1);
	    }
	  }
	  else {
	    fprintf(stderr,"Unable to map multiple statistical processes to GRIB1\n");
	    exit(1);
	  }
	}
	else {
	  switch (grid->md.stat_proc.proc_code[0]) {
	    case 0:
	    case 1:
	    case 4:
		switch (grid->md.stat_proc.proc_code[0]) {
/* average */
		  case 0:
		    *t_range=3;
		    break;
/* accumulation */
		  case 1:
		    *t_range=4;
		    break;
/* difference */
		  case 4:
		    *t_range=5;
		    break;
		}
		*p1=grid->md.fcst_time;
		*p2=mapStatisticalEndTime(msg,grid);
		if (grid->md.stat_proc.incr_length[0] == 0)
		  *n_avg=0;
		else {
		  fprintf(stderr,"Unable to map discrete processing to GRIB1\n");
		  exit(1);
		}
		break;
// maximum
	    case 2:
// minimum
	    case 3:
		*t_range=2;
		*p1=grid->md.fcst_time;
		*p2=mapStatisticalEndTime(msg,grid);
		if (grid->md.stat_proc.incr_length[0] == 0)
		  *n_avg=0;
		else {
		  fprintf(stderr,"Unable to map discrete processing to GRIB1\n");
		  exit(1);
		}
		break;
	    default:
// patch for NCEP grids
		if (grid->md.stat_proc.proc_code[0] == 255 && center == 7) {
 		  if (msg->disc == 0) {
		    if (grid->md.param_cat == 0) {
			switch (grid->md.param_num) {
			  case 4:
			  case 5:
			    *t_range=2;
			    *p1=grid->md.fcst_time;
			    *p2=mapStatisticalEndTime(msg,grid);
			    if (grid->md.stat_proc.incr_length[0] == 0)
				*n_avg=0;
			    else {
				fprintf(stderr,"Unable to map discrete processing to GRIB1\n");
				exit(1);
			    }
			    break;
			}
		    }
		  }
		}
		else {
		  fprintf(stderr,"Unable to map statistical process %d to GRIB1\n",grid->md.stat_proc.proc_code[0]);
		  exit(1);
		}
	  }
	}
	*n_missing=grid->md.stat_proc.nmiss;
	break;
    default:
	fprintf(stderr,"Unable to map time range for Product Definition Template %d into GRIB1\n",grid->md.pds_templ_num);
	exit(1);
  }
}

void packPDS(GRIB2Message *msg,int grid_number,unsigned char *grib1_buffer,size_t *offset)
{
  int level_type,level1,level2,p1,p2,t_range,n_avg,n_missing,D;
  static short warned=0;

  ParameterData pdata;
  pdata=mapParameterData(msg,grid_number);
// length of the PDS
  setBits(grib1_buffer,28,*offset,24);
// GRIB1 tables version number
  setBits(grib1_buffer,pdata.table_version,*offset+24,8);
// originating center ID
  setBits(grib1_buffer,msg->center_id,*offset+32,8);
// generating process ID
  setBits(grib1_buffer,msg->grids[grid_number].md.gen_proc,*offset+40,8);
// grid definition catalog number - set to 255 because GDS is to be included
  setBits(grib1_buffer,255,*offset+48,8);
// flag
  if (msg->grids[grid_number].md.bitmap == NULL) {
    setBits(grib1_buffer,0x80,*offset+56,8);
  }
  else {
    setBits(grib1_buffer,0xc0,*offset+56,8);
  }
// parameter code
  setBits(grib1_buffer,pdata.parameter_code,*offset+64,8);
  mapLevelData(&msg->grids[grid_number],&level_type,&level1,&level2,msg->center_id);
// level type code
  setBits(grib1_buffer,level_type,*offset+72,8);
  if (msg->grids[grid_number].md.lvl2_type == 255) {
    setBits(grib1_buffer,level1,*offset+80,16);
  }
  else {
    setBits(grib1_buffer,level1,*offset+80,8);
    setBits(grib1_buffer,level2,*offset+88,8);
  }
// year of century
  setBits(grib1_buffer,(msg->yr % 100),*offset+96,8);
// month
  setBits(grib1_buffer,msg->mo,*offset+104,8);
// day
  setBits(grib1_buffer,msg->dy,*offset+112,8);
// hour
  setBits(grib1_buffer,msg->time/10000,*offset+120,8);
// minute
  setBits(grib1_buffer,(msg->time/100 % 100),*offset+128,8);
// second
  if (msg->md.time_unit == 13) {
    fprintf(stderr,"Unable to indicate 'Second' for time unit in GRIB1\n");
  }
  else {
    setBits(grib1_buffer,msg->md.time_unit,*offset+136,8);
  }
  mapTimeRange(msg,&msg->grids[grid_number],&p1,&p2,&t_range,&n_avg,&n_missing,msg->center_id);
  if (t_range == 10) {
    setBits(grib1_buffer,p1,*offset+144,16);
  }
  else {
    setBits(grib1_buffer,p1,*offset+144,8);
    setBits(grib1_buffer,p2,*offset+152,8);
  }
  setBits(grib1_buffer,t_range,*offset+160,8);
// century of year
  setBits(grib1_buffer,(msg->yr/100)+1,*offset+192,8);
// originating sub-center ID
  setBits(grib1_buffer,msg->sub_center_id,*offset+200,8);
// decimal scale factor
  D=msg->md.D;
  if (D < 0) {
    D=-D+0x8000;
  }
  setBits(grib1_buffer,D,*offset+208,16);
  (*offset)+=224;
  if (msg->md.ens_type >= 0) {
// length of the PDS
    setBits(grib1_buffer,43,*offset-224,24);
    setBits(grib1_buffer,msg->md.ens_type,*offset+96,8);
    setBits(grib1_buffer,msg->md.perturb_num,*offset+104,8);
    setBits(grib1_buffer,msg->md.nfcst_in_ensemble,*offset+112,8);
    (*offset)+=120;
    if (warned == 0) {
	fprintf(stderr,"Notice: the 'Ensemble type code', the 'Perturbation Number', and the\n");
	fprintf(stderr,"'Number of forecasts in ensemble' from Product Definition Template 4.1 and/or\n");
	fprintf(stderr,"Product Definition Template 4.12 have been packed in octets 41, 42, and 43 of\n");
	fprintf(stderr,"the GRIB1 Product Definition Section\n");
	warned=1;
    }
  }
  else if (msg->md.derived_fcst_code >= 0) {
// length of the PDS
    setBits(grib1_buffer,42,*offset-224,24);
    setBits(grib1_buffer,msg->md.derived_fcst_code,*offset+96,8);
    setBits(grib1_buffer,msg->md.nfcst_in_ensemble,*offset+104,8);
    (*offset)+=112;
    if (warned == 0) {
	fprintf(stderr,"Notice: the 'Derived forecast code' and the 'Number of forecasts in ensemble'\n");
	fprintf(stderr,"from Product Definition Template 4.2 and/or Product Definition Template 4.12\n");
	fprintf(stderr,"have been packed in octets 41 and 42 of the GRIB1 Product Definition Section\n");
	warned=1;
    }
  }
  else if (msg->md.spatial_proc.type >= 0) {
// length of the PDS
    setBits(grib1_buffer,43,*offset-224,24);
    setBits(grib1_buffer,msg->md.spatial_proc.stat_proc,*offset+96,8);
    setBits(grib1_buffer,msg->md.spatial_proc.type,*offset+104,8);
    setBits(grib1_buffer,msg->md.spatial_proc.num_points,*offset+112,8);
    (*offset)+=120;
    if (warned == 0) {
	fprintf(stderr,"Notice: the Spatial processing codes: 'statistical process', 'type' and\n");
	fprintf(stderr,"'number of data points' from Product Definition Template 4.15 have been\n");
	fprintf(stderr,"packed in octets 41, 42, and 43 of the GRIB1 Product Definition Section\n");
	warned=1;
    }
  }
}

void packGDS(GRIB2Message *msg,int grid_number,unsigned char *grib1_buffer,size_t *offset)
{
  int rescomp=0,sign,value;

// NV
  setBits(grib1_buffer,255,*offset+24,8);
// PV
  setBits(grib1_buffer,255,*offset+32,8);
  switch (msg->md.gds_templ_num) {
    case 0:
// length of the GDS
	setBits(grib1_buffer,32,*offset,24);
// data representation
	setBits(grib1_buffer,0,*offset+40,8);
// Ni
	setBits(grib1_buffer,msg->md.nx,*offset+48,16);
// Nj
	setBits(grib1_buffer,msg->md.ny,*offset+64,16);
// first latitude
	value=msg->md.slat*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+80,1);
	  setBits(grib1_buffer,value,*offset+81,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+80,24);
// first longitude
	value=msg->md.slon*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+104,1);
	  setBits(grib1_buffer,value,*offset+105,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+104,24);
// resolution and component flags
	if ((msg->md.rescomp&0x20) == 0x20)
	  rescomp|=0x80;
	if (msg->md.earth_shape == 2)
	  rescomp|=0x40;
	if ((msg->md.rescomp&0x8) == 0x8)
	  rescomp|=0x8;
	setBits(grib1_buffer,rescomp,*offset+128,8);
// last latitude
	value=msg->md.lats.elat*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+136,1);
	  setBits(grib1_buffer,value,*offset+137,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+136,24);
// last longitude
	value=msg->md.lons.elon*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+160,1);
	  setBits(grib1_buffer,value,*offset+161,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+160,24);
// Di increment
	value=msg->md.xinc.loinc*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+184,1);
	  setBits(grib1_buffer,value,*offset+185,15);
	}
	else
	  setBits(grib1_buffer,value,*offset+184,16);
// Dj increment
	value=msg->md.yinc.lainc*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+200,1);
	  setBits(grib1_buffer,value,*offset+201,15);
	}
	else
	  setBits(grib1_buffer,value,*offset+200,16);
// scanning mode
	setBits(grib1_buffer,msg->md.scan_mode,*offset+216,8);
// reserved
	setBits(grib1_buffer,0,*offset+224,32);
	(*offset)+=256;
	break;
    case 30:
// length of the GDS
	setBits(grib1_buffer,42,*offset,24);
// data representation
	setBits(grib1_buffer,3,*offset+40,8);
// Nx
	setBits(grib1_buffer,msg->md.nx,*offset+48,16);
// Ny
	setBits(grib1_buffer,msg->md.ny,*offset+64,16);
// first latitude
	value=msg->md.slat*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+80,1);
	  setBits(grib1_buffer,value,*offset+81,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+80,24);
// first longitude
	value=msg->md.slon*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+104,1);
	  setBits(grib1_buffer,value,*offset+105,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+104,24);
// resolution and component flags
	if ((msg->md.rescomp&0x20) == 0x20)
	  rescomp|=0x80;
	if (msg->md.earth_shape == 2)
	  rescomp|=0x40;
	if ((msg->md.rescomp&0x8) == 0x8)
	  rescomp|=0x8;
	setBits(grib1_buffer,rescomp,*offset+128,8);
// LoV
	value=msg->md.lons.lov*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+136,1);
	  setBits(grib1_buffer,value,*offset+137,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+136,24);
// Dx
	value=msg->md.xinc.dxinc+0.5;
	setBits(grib1_buffer,value,*offset+160,24);
// Dy
	value=msg->md.yinc.dyinc+0.5;
	setBits(grib1_buffer,value,*offset+184,24);
// projection center flag
	setBits(grib1_buffer,msg->md.proj_flag,*offset+208,8);
// scanning mode
	setBits(grib1_buffer,msg->md.scan_mode,*offset+216,8);
// latin1
	value=msg->md.latin1*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+224,1);
	  setBits(grib1_buffer,value,*offset+225,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+224,24);
// latin2
	value=msg->md.latin2*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+248,1);
	  setBits(grib1_buffer,value,*offset+249,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+248,24);
// latitude of southern pole of projection
	value=msg->md.splat*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+272,1);
	  setBits(grib1_buffer,value,*offset+273,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+272,24);
// longitude of southern pole of projection
	value=msg->md.splon*1000.;
	if (value < 0.) {
	  value=-value;
	  setBits(grib1_buffer,1,*offset+296,1);
	  setBits(grib1_buffer,value,*offset+297,23);
	}
	else
	  setBits(grib1_buffer,value,*offset+296,24);
// reserved
	setBits(grib1_buffer,0,*offset+320,16);
	(*offset)+=336;
	break;
    default:
	fprintf(stderr,"Unable to map Grid Definition Template %d into GRIB1\n",msg->md.gds_templ_num);
	exit(1);
  }
}

void packBMS(GRIB2Message *msg,int grid_number,unsigned char *grib1_buffer,size_t *offset,size_t num_points)
{
  int length=6+(num_points+7)/8;
  int ub=8-(num_points % 8);
  size_t n,off;

// length of the BMS
  setBits(grib1_buffer,length,*offset,24);
// unused bits at end of section
  setBits(grib1_buffer,ub,*offset+24,8);
// table reference
  setBits(grib1_buffer,0,*offset+32,16);
// the bitmap
  off=*offset+48;
  for (n=0; n < num_points; ++n) {
    setBits(grib1_buffer,msg->grids[grid_number].md.bitmap[n],off++,1);
  }
  (*offset)+=length*8;
}

int ieee2ibm(double ieee)
{
  int ibm_real=0;
  unsigned char *ir=(unsigned char *)&ibm_real;
  int sign=0,fr=0;
  int exp=64;
  const double full=0xffffff;
  size_t size=sizeof(int)*8,off=0;

  if (ieee != 0.) {
    if (ieee < 0.) {
      sign=1;
      ieee=-ieee;
    }
    ieee/=pow(2.,-24.);
    while (exp > 0 && ieee < full) {
      ieee*=16.;
      exp--;
    }
    while (ieee > full) {
      ieee/=16.;
      exp++;
    }
    fr=ieee+0.5;
    if (size > 32) {
      off=size-32;
      setBits(ir,0,0,off);
    }
    setBits(ir,sign,off,1);
    setBits(ir,exp,off+1,7);
    setBits(ir,fr,off+8,24);
  }
  return ibm_real;
}

void packBDS(GRIB2Message *msg,int grid_number,unsigned char *grib1_buffer,size_t *offset,int *pvals,size_t num_to_pack,size_t pack_width)
{
  int length=11+(num_to_pack*pack_width+7)/8;
  size_t m,off;
  int E,ibm_rep;

// length of the BDS
  setBits(grib1_buffer,length,*offset,24);
// flag
  setBits(grib1_buffer,0,*offset+24,4);
// unused bits
  setBits(grib1_buffer,(length-11)*8-(num_to_pack*pack_width),*offset+28,4);
// scale factor E
  E=msg->grids[grid_number].md.E;
  if (E < 0)
    E=-E+0x8000;
  setBits(grib1_buffer,E,*offset+32,16);
// Reference value
  ibm_rep=ieee2ibm(msg->grids[grid_number].md.R*pow(10.,msg->grids[grid_number].md.D));
  memcpy(&grib1_buffer[(*offset+48)/8],&ibm_rep,4);
// width in bits of each packed value
  setBits(grib1_buffer,pack_width,*offset+80,8);
// packed data values
  off=*offset+88;
  for (m=0; m < num_to_pack; m++) {
    setBits(grib1_buffer,pvals[m],off,pack_width);
    off+=pack_width;
  }
}

int main(int argc,char **argv)
{
  if (argc != 3) {
    fprintf(stderr,"usage: %s GRIB2_file_name GRIB1_file_name\n",argv[0]);
    exit(1);
  }
  FILE *fp;
  if ( (fp=fopen(argv[1],"rb")) == NULL) {
    fprintf(stderr,"Error opening input file %s\n",argv[1]);
    exit(1);
  }
  FILE *ofp;
  if ( (ofp=fopen(argv[2],"wb")) == NULL) {
    fprintf(stderr,"Error opening output file %s\n",argv[2]);
    exit(1);
  }
  GRIB2Message grib2_msg;
  initialize(&grib2_msg);
  int status;
  size_t nmsg=0;
  size_t ngrid=0;
  int max_length=0;
  unsigned char *grib1_buffer=NULL;
  char *head="GRIB",*tail="7777";
  while ( (status=unpackgrib2(fp,&grib2_msg)) == 0) {
    ++nmsg;
    for (size_t n=0; n < grib2_msg.num_grids; ++n) {
// calculate the octet length of the GRIB1 grid (minus the Indicator and End
// Sections, which are both fixed in length
	int length;
	switch (grib2_msg.md.pds_templ_num) {
	  case 0:
	  case 8:
	  {
	    length=28;
	    break;
	  }
	  case 1:
	  case 11:
	  {
	    length=43;
	    break;
	  }
	  case 2:
	  case 12:
	  {
	    length=42;
	    break;
	  }
	  case 15:
	  {
	    length=43;
	    break;
	  }
	  default:
	  {
	    fprintf(stderr,"Unable to map Product Definition Template %d into GRIB1\n",grib2_msg.md.pds_templ_num);
	    exit(1);
	  }
	}
	size_t num_points;
	switch (grib2_msg.md.gds_templ_num) {
	  case 0:
	  {
	    length+=32;
	    num_points=grib2_msg.md.nx*grib2_msg.md.ny;
	    break;
	  }
	  case 30:
	  {
	    length+=42;
	    num_points=grib2_msg.md.nx*grib2_msg.md.ny;
	    break;
	  }
	  default:
	  {
	    fprintf(stderr,"Unable to map Grid Definition Template %d into GRIB1\n",grib2_msg.md.gds_templ_num);
	    exit(1);
	  }
	}
	size_t num_to_pack;
	if (grib2_msg.grids[n].md.bitmap != NULL) {
	  length+=6+(num_points+7)/8;
	  num_to_pack=0;
	  for (size_t m=0; m < num_points; ++m) {
	    if (grib2_msg.grids[n].md.bitmap[m] == 1) {
		++num_to_pack;
	    }
	  }
	}
	else {
	  num_to_pack=num_points;
	}
	int *pvals=(int *)malloc(sizeof(int)*num_to_pack);
	size_t max_pack=0;
	size_t cnt=0;
	for (size_t m=0; m < num_points; ++m) {
	  if (grib2_msg.grids[n].gridpoints[m] != GRIB_MISSING_VALUE) {
	    if (cnt == num_to_pack) {
		fprintf(stderr,"Error: conflicting number of missing gridpoints\n");
		exit(1);
	    }
	    pvals[cnt]=lround((grib2_msg.grids[n].gridpoints[m]-grib2_msg.grids[n].md.R)*pow(10.,grib2_msg.grids[n].md.D)/pow(2.,grib2_msg.grids[n].md.E));
	    if (pvals[cnt] > max_pack) {
		max_pack=pvals[cnt];
	    }
	    ++cnt;
	  }
	}
	size_t pack_width=1;
	while (pow(2.,pack_width)-1 < max_pack) {
	  ++pack_width;
	}
	length+=11+(num_to_pack*pack_width+7)/8;
// allocate enough memory for the GRIB1 buffer
	if (length > max_length) {
	  if (grib1_buffer != NULL) {
	    free(grib1_buffer);
	  }
	  grib1_buffer=(unsigned char *)malloc(length*sizeof(unsigned char));
	  max_length=length;
	}
	size_t offset=0;
// pack the Product Definition Section
	packPDS(&grib2_msg,n,grib1_buffer,&offset);
// pack the Grid Definition Section
	packGDS(&grib2_msg,n,grib1_buffer,&offset);
// pack the Bitmap Section, if it exists
	if (grib2_msg.grids[n].md.bitmap != NULL) {
	  packBMS(&grib2_msg,n,grib1_buffer,&offset,num_points);
	}
// pack the Binary Data Section
	packBDS(&grib2_msg,n,grib1_buffer,&offset,pvals,num_to_pack,pack_width);
	free(pvals);
// output the GRIB1 grid
	fwrite(head,1,4,ofp);
	unsigned char dum[3];
	setBits(dum,length+12,0,24);
	fwrite(dum,1,3,ofp);
	dum[0]=1;
	fwrite(dum,1,1,ofp);
	fwrite(grib1_buffer,1,length,ofp);
	fwrite(tail,1,4,ofp);
	++ngrid;
    }
  }
  if (status != -1) {
    printf("Read error after %d messages\n",nmsg);
  }
  printf("Number of GRIB1 grids written to output: %d\n",ngrid);
  fclose(fp);
  fclose(ofp);
}
