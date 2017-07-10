/*
** File: grib1to2.c
**
** Author:  Bob Dattore
**          NCAR/DSS
**          dattore@ucar.edu
**          (303) 497-1825
**
** Purpose: to provide a simple C program for converting from GRIB1 to GRIB2
**
** Revision History:
**   20 May 2017 - first version
**   10 Jul 2017 - convert Mercator grids; always include bitmap section (6)
**
** Example compile command:
**    % cc -std=c99 -o grib1to2 grib1to2.c
**
** If the compiler complains about the "pow" function being an undefined symbol,
**   include the math library in the compile
**      e.g. % cc -std=c99 -o grib1to2 grib1to2.c -lm
**
** To use the program:
**    % grib1to2 <name of GRIB1 file to convert> <name of GRIB2 file to create>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "unpackgrib1.c"

/* set_bits sets the contents of the various GRIB octets
**   buf is the GRIB buffer as a stream of bytes
**   src is the value of the octet(s) to set
**   off is the offset in BITS from the beginning of the buffer to the beginning
**       of the octet(s) to be packed
**   bits is the number of BITS to pack - will be a multiple of 8 since GRIB
**       octets are 8 bits long
*/
void set_bits(unsigned char *buf,int src,size_t off,size_t bits)
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

void map_GRIB1_parameter_into_GRIB2(GRIBMessage *grib_msg,int *discipline,int *parameter_category,int *parameter_number)
{
  *discipline=255;*parameter_category=255,*parameter_number=255;
  switch (grib_msg->param) {
    case 1:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 228) {
	  *discipline=0;*parameter_category=7;*parameter_number=7;
	}
	else {
	  *discipline=0;*parameter_category=3;*parameter_number=0;
	}
	return;
    }
    case 2:
    {
	*discipline=0;*parameter_category=3;*parameter_number=1;
	return;
    }
    case 3:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 228) {
	  *discipline=10;*parameter_category=0;*parameter_number=17;
	}
	else {
	  *discipline=0;*parameter_category=3;*parameter_number=2;
	}
	return;
    }
    case 4:
    {
	*discipline=0;*parameter_category=2;*parameter_number=14;
	return;
    }
    case 5:
    {
	*discipline=0;*parameter_category=3;*parameter_number=3;
	return;
    }
    case 6:
    {
	*discipline=0;*parameter_category=3;*parameter_number=4;
	return;
    }
    case 7:
    {
	*discipline=0;*parameter_category=3;*parameter_number=5;
	return;
    }
    case 8:
    {
	if (grib_msg->center_id == 78 && grib_msg->table_ver == 174) {
	  *discipline=2;*parameter_category=0;*parameter_number=34;
	}
	else {
	  *discipline=0;*parameter_category=3;*parameter_number=6;
	}
	return;
    }
    case 9:
    {
	*discipline=0;*parameter_category=3;*parameter_number=7;
	return;
    }
    case 10:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 200) {
	  *discipline=0;*parameter_category=14;*parameter_number=2;
	}
	else {
	  *discipline=0;*parameter_category=14;*parameter_number=0;
	}
	return;
    }
    case 11:
    {
	*discipline=0;*parameter_category=0;*parameter_number=0;
	return;
    }
    case 12:
    {
	*discipline=0;*parameter_category=0;*parameter_number=1;
	return;
    }
    case 13:
    {
	*discipline=0;*parameter_category=0;*parameter_number=2;
	return;
    }
    case 14:
    {
	*discipline=0;*parameter_category=0;*parameter_number=3;
	return;
    }
    case 15:
    {
	*discipline=0;*parameter_category=0;*parameter_number=4;
	return;
    }
    case 16:
    {
	*discipline=0;*parameter_category=0;*parameter_number=5;
	return;
    }
    case 17:
    {
	*discipline=0;*parameter_category=0;*parameter_number=6;
	return;
    }
    case 18:
    {
	*discipline=0;*parameter_category=0;*parameter_number=7;
	return;
    }
    case 19:
    {
	*discipline=0;*parameter_category=0;*parameter_number=8;
	return;
    }
    case 20:
    {
	*discipline=0;*parameter_category=19;*parameter_number=0;
	return;
    }
    case 21:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=0;*parameter_number=28;
	}
	else {
	  *discipline=0;*parameter_category=15;*parameter_number=6;
	}
	return;
    }
    case 22:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=3;*parameter_number=31;
	}
	else {
	  *discipline=0;*parameter_category=15;*parameter_number=7;
	}
	return;
    }
    case 23:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=2;*parameter_number=45;
	}
	else {
	  *discipline=0;*parameter_category=15;*parameter_number=8;
	}
	return;
    }
    case 24:
    {
	*discipline=0;*parameter_category=7;*parameter_number=0;
	return;
    }
    case 25:
    {
	*discipline=0;*parameter_category=0;*parameter_number=9;
	return;
    }
    case 26:
    {
	*discipline=0;*parameter_category=3;*parameter_number=8;
	return;
    }
    case 27:
    {
	*discipline=0;*parameter_category=3;*parameter_number=9;
	return;
    }
    case 28:
    {
	*discipline=10;*parameter_category=0;*parameter_number=0;
	return;
    }
    case 29:
    {
	*discipline=10;*parameter_category=0;*parameter_number=1;
	return;
    }
    case 30:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 203) {
	  *discipline=0;*parameter_category=7;*parameter_number=8;
	}
	else {
	  *discipline=10;*parameter_category=0;*parameter_number=2;
	}
	return;
    }
    case 31:
    {
	*discipline=0;*parameter_category=2;*parameter_number=0;
	return;
    }
    case 32:
    {
	*discipline=0;*parameter_category=2;*parameter_number=1;
	return;
    }
    case 33:
    {
	if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 201) {
	    *discipline=0;*parameter_category=1;*parameter_number=82;
	  }
	  else if (grib_msg->table_ver == 203) {
	    *discipline=0;*parameter_category=2;*parameter_number=46;
	  }
	}
	else {
	  *discipline=0;*parameter_category=2;*parameter_number=2;
	}
	return;
    }
    case 34:
    {
	*discipline=0;*parameter_category=2;*parameter_number=3;
	return;
    }
    case 35:
    {
	*discipline=0;*parameter_category=2;*parameter_number=4;
	return;
    }
    case 36:
    {
	*discipline=0;*parameter_category=2;*parameter_number=5;
	return;
    }
    case 37:
    {
	*discipline=0;*parameter_category=2;*parameter_number=6;
	return;
    }
    case 38:
    {
	*discipline=0;*parameter_category=2;*parameter_number=7;
	return;
    }
    case 39:
    {
	*discipline=0;*parameter_category=2;*parameter_number=8;
	return;
    }
    case 40:
    {
	*discipline=0;*parameter_category=2;*parameter_number=9;
	return;
    }
    case 41:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=1;*parameter_number=78;
	}
	else {
	  *discipline=0;*parameter_category=2;*parameter_number=10;
	}
	return;
    }
    case 42:
    {
	*discipline=0;*parameter_category=2;*parameter_number=11;
	return;
    }
    case 43:
    {
	*discipline=0;*parameter_category=2;*parameter_number=12;
	return;
    }
    case 44:
    {
	*discipline=0;*parameter_category=2;*parameter_number=13;
	return;
    }
    case 45:
    {
	*discipline=0;*parameter_category=2;*parameter_number=15;
	return;
    }
    case 46:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 202) {
	  *discipline=0;*parameter_category=3;*parameter_number=20;
	}
	else {
	  *discipline=0;*parameter_category=2;*parameter_number=16;
	}
	return;
    }
    case 47:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 202) {
	  *discipline=0;*parameter_category=3;*parameter_number=24;
	}
	else {
	  *discipline=10;*parameter_category=1;*parameter_number=0;
	}
	return;
    }
    case 48:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 202) {
	  *discipline=0;*parameter_category=3;*parameter_number=21;
	}
	else {
	  *discipline=10;*parameter_category=1;*parameter_number=1;
	}
	return;
    }
    case 49:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 202) {
	  *discipline=0;*parameter_category=3;*parameter_number=22;
	}
	else {
	  *discipline=10;*parameter_category=1;*parameter_number=2;
	}
	return;
    }
    case 50:
    {
	*discipline=10;*parameter_category=1;*parameter_number=3;
	return;
    }
    case 51:
    {
	*discipline=0;*parameter_category=1;*parameter_number=0;
	return;
    }
    case 52:
    {
	*discipline=0;*parameter_category=1;*parameter_number=1;
	return;
    }
    case 53:
    {
	*discipline=0;*parameter_category=1;*parameter_number=2;
	return;
    }
    case 54:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=64;
	}
	else {
	  *discipline=0;*parameter_category=1;*parameter_number=3;
	}
	return;
    }
    case 55:
    {
	*discipline=0;*parameter_category=1;*parameter_number=4;
	return;
    }
    case 56:
    {
	*discipline=0;*parameter_category=1;*parameter_number=5;
	return;
    }
    case 57:
    {
	*discipline=0;*parameter_category=1;*parameter_number=6;
	return;
    }
    case 58:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=70;
	}
	else {
	  *discipline=0;*parameter_category=6;*parameter_number=0;
	}
	return;
    }
    case 59:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=7;*parameter_number=6;
	}
	else {
	  *discipline=0;*parameter_category=1;*parameter_number=7;
	}
	return;
    }
    case 60:
    {
	*discipline=0;*parameter_category=19;*parameter_number=2;
	return;
    }
    case 61:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 202) {
	  *discipline=2;*parameter_category=0;*parameter_number=28;
	}
	else {
	  *discipline=0;*parameter_category=1;*parameter_number=8;
	}
	return;
    }
    case 62:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 202) {
	  *discipline=2;*parameter_category=0;*parameter_number=32;
	}
	else {
	  *discipline=0;*parameter_category=1;*parameter_number=9;
	}
	return;
    }
    case 63:
    {
	*discipline=0;*parameter_category=1;*parameter_number=10;
	return;
    }
    case 64:
    {
	*discipline=0;*parameter_category=1;*parameter_number=12;
	return;
    }
    case 65:
    {
	*discipline=0;*parameter_category=1;*parameter_number=13;
	return;
    }
    case 66:
    {
	*discipline=0;*parameter_category=1;*parameter_number=11;
	return;
    }
    case 67:
    {
	*discipline=0;*parameter_category=19;*parameter_number=3;
	return;
    }
    case 68:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=6;*parameter_number=26;
	}
	else {
	  *discipline=10;*parameter_category=4;*parameter_number=2;
	}
	return;
    }
    case 69:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=6;*parameter_number=27;
	}
	else {
	  *discipline=10;*parameter_category=4;*parameter_number=0;
	}
	return;
    }
    case 70:
    {
	*discipline=10;*parameter_category=4;*parameter_number=1;
	return;
    }
    case 71:
    {
	*discipline=0;*parameter_category=6;*parameter_number=1;
	return;
    }
    case 72:
    {
	*discipline=0;*parameter_category=6;*parameter_number=2;
	return;
    }
    case 73:
    {
	*discipline=0;*parameter_category=6;*parameter_number=3;
	return;
    }
    case 74:
    {
	*discipline=0;*parameter_category=6;*parameter_number=4;
	return;
    }
    case 75:
    {
	if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 128) {
	    *discipline=0;*parameter_category=1;*parameter_number=85;
	  }
	  else if (grib_msg->table_ver == 202) {
	    *discipline=2;*parameter_category=0;*parameter_number=29;
	  }
	}
	else {
	  *discipline=0;*parameter_category=6;*parameter_number=5;
	}
	return;
    }
    case 76:
    {
	if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=1;*parameter_number=69;
	  }
	  else if (grib_msg->table_ver == 128) {
	    *discipline=0;*parameter_category=1;*parameter_number=86;
	  }
	  else if (grib_msg->table_ver == 202) {
	    *discipline=2;*parameter_category=0;*parameter_number=30;
	  }
	}
	else {
	  *discipline=0;*parameter_category=6;*parameter_number=6;
	}
	return;
    }
    case 77:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=2;*parameter_number=32;
	}
	else {
	  *discipline=0;*parameter_category=7;*parameter_number=1;
	}
	return;
    }
    case 78:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 202) {
	  *discipline=2;*parameter_category=0;*parameter_number=31;
	}
	else {
	  *discipline=0;*parameter_category=1;*parameter_number=14;
	}
	return;
    }
    case 79:
    {
	*discipline=0;*parameter_category=1;*parameter_number=15;
	return;
    }
    case 80:
    {
	*discipline=10;*parameter_category=3;*parameter_number=0;
	return;
    }
    case 81:
    {
	*discipline=1;*parameter_category=2;*parameter_number=8;
	return;
    }
    case 82:
    {
	*discipline=10;*parameter_category=3;*parameter_number=1;
	return;
    }
    case 83:
    {
	*discipline=2;*parameter_category=0;*parameter_number=1;
	return;
    }
    case 84:
    {
	*discipline=0;*parameter_category=19;*parameter_number=1;
	return;
    }
    case 85:
    {
	*discipline=2;*parameter_category=0;*parameter_number=2;
	return;
    }
    case 86:
    {
	*discipline=2;*parameter_category=0;*parameter_number=3;
	return;
    }
    case 87:
    {
	*discipline=2;*parameter_category=0;*parameter_number=4;
	return;
    }
    case 88:
    {
	*discipline=10;*parameter_category=4;*parameter_number=3;
	return;
    }
    case 89:
    {
	*discipline=0;*parameter_category=3;*parameter_number=10;
	return;
    }
    case 90:
    {
	*discipline=2;*parameter_category=0;*parameter_number=5;
	return;
    }
    case 91:
    {
	*discipline=1;*parameter_category=2;*parameter_number=7;
	return;
    }
    case 92:
    {
	*discipline=10;*parameter_category=2;*parameter_number=1;
	return;
    }
    case 93:
    {
	*discipline=10;*parameter_category=2;*parameter_number=2;
	return;
    }
    case 94:
    {
	*discipline=10;*parameter_category=2;*parameter_number=3;
	return;
    }
    case 95:
    {
	*discipline=10;*parameter_category=2;*parameter_number=4;
	return;
    }
    case 96:
    {
	*discipline=10;*parameter_category=2;*parameter_number=5;
	return;
    }
    case 97:
    {
	*discipline=10;*parameter_category=2;*parameter_number=6;
	return;
    }
    case 98:
    {
	*discipline=10;*parameter_category=2;*parameter_number=7;
	return;
    }
    case 99:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 203) {
	  *discipline=0;*parameter_category=19;*parameter_number=25;
	}
	else {
	  *discipline=0;*parameter_category=1;*parameter_number=16;
	}
	return;
    }
    case 100:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=1;*parameter_number=77;
	}
	else {
	  *discipline=10;*parameter_category=0;*parameter_number=3;
	}
	return;
    }
    case 101:
    {
	*discipline=10;*parameter_category=0;*parameter_number=4;
	return;
    }
    case 102:
    {
	*discipline=10;*parameter_category=0;*parameter_number=5;
	return;
    }
    case 103:
    {
	*discipline=10;*parameter_category=0;*parameter_number=6;
	return;
    }
    case 104:
    {
	*discipline=10;*parameter_category=0;*parameter_number=7;
	return;
    }
    case 105:
    {
	*discipline=10;*parameter_category=0;*parameter_number=8;
	return;
    }
    case 106:
    {
	*discipline=10;*parameter_category=0;*parameter_number=9;
	return;
    }
    case 107:
    {
	*discipline=10;*parameter_category=0;*parameter_number=10;
	return;
    }
    case 108:
    {
	*discipline=10;*parameter_category=0;*parameter_number=11;
	return;
    }
    case 109:
    {
	if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 162) {
	    *discipline=0;*parameter_category=0;*parameter_number=20;
	  }
	  else if (grib_msg->table_ver == 228) {
	    *discipline=0;*parameter_category=6;*parameter_number=13;
	  }
	}
	else {
	  *discipline=10;*parameter_category=0;*parameter_number=12;
	}
	return;
    }
    case 110:
    {
	*discipline=10;*parameter_category=0;*parameter_number=13;
	return;
    }
    case 111:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=1;*parameter_number=76;
	}
	else {
	  *discipline=0;*parameter_category=4;*parameter_number=0;
	}
	return;
    }
    case 112:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=1;*parameter_number=55;
	}
	else {
	  *discipline=0;*parameter_category=5;*parameter_number=0;
	}
	return;
    }
    case 113:
    {
	*discipline=0;*parameter_category=4;*parameter_number=1;
	return;
    }
    case 114:
    {
	*discipline=0;*parameter_category=5;*parameter_number=1;
	return;
    }
    case 115:
    {
	*discipline=0;*parameter_category=5;*parameter_number=2;
	return;
    }
    case 116:
    {
	*discipline=0;*parameter_category=4;*parameter_number=2;
	return;
    }
    case 117:
    {
	*discipline=0;*parameter_category=4;*parameter_number=3;
	return;
    }
    case 118:
    {
	*discipline=0;*parameter_category=4;*parameter_number=4;
	return;
    }
    case 119:
    {
	*discipline=0;*parameter_category=4;*parameter_number=5;
	return;
    }
    case 120:
    {
	*discipline=0;*parameter_category=4;*parameter_number=6;
	return;
    }
    case 121:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 228) {
	  *discipline=0;*parameter_category=7;*parameter_number=2;
	}
	else {
	  *discipline=0;*parameter_category=0;*parameter_number=10;
	}
	return;
    }
    case 122:
    {
	*discipline=0;*parameter_category=0;*parameter_number=11;
	return;
    }
    case 123:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 228) {
	  *discipline=0;*parameter_category=7;*parameter_number=4;
	}
	else {
	  *discipline=0;*parameter_category=2;*parameter_number=20;
	}
	return;
    }
    case 124:
    {
	*discipline=0;*parameter_category=2;*parameter_number=17;
	return;
    }
    case 125:
    {
	*discipline=0;*parameter_category=2;*parameter_number=18;
	return;
    }
    case 126:
    {
	*discipline=0;*parameter_category=2;*parameter_number=19;
	return;
    }
    case 131:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=1;*parameter_number=70;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=43;
	  }
	}
	return;
    }
    case 132:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=7;*parameter_number=11;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=6;*parameter_number=21;
	  }
	}
	return;
    }
    case 133:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=7;*parameter_number=2;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=44;
	  }
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=1;*parameter_number=61;
	}
	return;
    }
    case 134:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=7;*parameter_number=5;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=6;*parameter_number=16;
	  }
	}
	return;
    }
    case 135:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=1;*parameter_number=38;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=21;
	  }
	}
	return;
    }
    case 136:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=2;*parameter_number=25;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=69;
	  }
	}
	return;
    }
    case 137:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=70;
	  }
	  else if (grib_msg->table_ver == 131) {
	    *discipline=0;*parameter_category=17;*parameter_number=0;
	  }
	}
	return;
    }
    case 138:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 129) {
	  *discipline=0;*parameter_category=1;*parameter_number=45;
	}
	return;
    }
    case 139:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 129) {
	  *discipline=0;*parameter_category=1;*parameter_number=46;
	}
	return;
    }
    case 140:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=1;*parameter_number=33;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=6;*parameter_number=20;
	  }
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 203) {
	  *discipline=0;*parameter_category=7;*parameter_number=3;
	}
	return;
    }
    case 141:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=34;
	}
	return;
    }
    case 142:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=35;
	}
	return;
    }
    case 143:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=36;
	}
	return;
    }
    case 144:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=2;*parameter_category=0;*parameter_number=9;
	  }
	  else if (grib_msg->table_ver == 128) {
	    *discipline=10;*parameter_category=3;*parameter_number=2;
	  }
	}
	return;
    }
    case 145:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 131) {
	  *discipline=0;*parameter_category=1;*parameter_number=41;
	}
	return;
    }
    case 146:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=6;*parameter_number=15;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 200) {
	  *discipline=0;*parameter_category=6;*parameter_number=15;
	}
	return;
    }
    case 147:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=3;*parameter_number=16;
	}
	else if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 201) {
	    *discipline=0;*parameter_category=19;*parameter_number=24;
	  }
	  else if (grib_msg->table_ver == 254) {
	    *discipline=0;*parameter_category=2;*parameter_number=27;
	  }
	}
	return;
    }
    case 148:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=3;*parameter_number=17;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 254) {
	  *discipline=0;*parameter_category=2;*parameter_number=28;
	}
	return;
    }
    case 152:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=19;*parameter_number=11;
	}
	return;
    }
    case 153:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=22;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=2;*parameter_number=31;
	}
	return;
    }
    case 154:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=14;*parameter_number=1;
	}
	return;
    }
    case 155:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=2;*parameter_category=0;*parameter_number=10;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 200) {
	  *discipline=2;*parameter_category=0;*parameter_number=10;
	}
	return;
    }
    case 156:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=7;*parameter_number=7;
	}
	return;
    }
    case 157:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=7;*parameter_number=6;
	}
	return;
    }
    case 158:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=19;*parameter_number=11;
	}
	return;
    }
    case 159:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=0;*parameter_category=19;*parameter_number=17;
	}
	return;
    }
    case 160:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=4;*parameter_number=53;
	  }
	  else if (grib_msg->table_ver == 130) {
	    *discipline=2;*parameter_category=3;*parameter_number=5;
	  }
	}
	return;
    }
    case 163:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=5;*parameter_number=8;
	}
	return;
    }
    case 170:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=1;*parameter_number=24;
	  }
	  else if (grib_msg->table_ver == 130) {
	    *discipline=0;*parameter_category=19;*parameter_number=18;
	  }
	}
	return;
    }
    case 171:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=1;*parameter_number=25;
	  }
	  else if (grib_msg->table_ver == 130) {
	    *discipline=2;*parameter_category=3;*parameter_number=6;
	  }
	}
	else if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 201) {
	    *discipline=0;*parameter_category=0;*parameter_number=19;
	  }
	  else if (grib_msg->table_ver == 228) {
	    *discipline=2;*parameter_category=0;*parameter_number=26;
	  }
	}
	return;
    }
    case 172:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=2;*parameter_number=26;
	}
	return;
    }
    case 174:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 140) {
	  *discipline=0;*parameter_category=6;*parameter_number=25;
	}
	return;
    }
    case 178:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=23;
	}
	return;
    }
    case 180:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=0;*parameter_category=1;*parameter_number=17;
	}
	else if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 128) {
	    *discipline=0;*parameter_category=2;*parameter_number=38;
	  }
	  else if (grib_msg->table_ver == 202) {
	    *discipline=0;*parameter_category=14;*parameter_number=1;
	  }
	}
	return;
    }
    case 181:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=0;*parameter_number=15;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=2;*parameter_number=37;
	}
	return;
    }
    case 182:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=0;*parameter_number=28;
	}
	return;
    }
    case 184:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=0;*parameter_category=19;*parameter_number=19;
	}
	return;
    }
    case 189:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=0;*parameter_number=15;
	}
	return;
    }
    case 190:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=7;*parameter_number=8;
	}
	return;
    }
    case 191:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 133) {
	  *discipline=0;*parameter_category=6;*parameter_number=33;
	}
	return;
    }
    case 192:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 133) {
	  *discipline=10;*parameter_category=191;*parameter_number=1;
	}
	return;
    }
    case 193:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 131) {
	  *discipline=0;*parameter_category=0;*parameter_number=21;
	}
	return;
    }
    case 194:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1,*parameter_number=39;
	}
	return;
    }
    case 195:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 128) {
	  *discipline=10;*parameter_category=4,*parameter_number=4;
	}
	return;
    }
    case 196:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=2;*parameter_number=27;
	  }
	  else if (grib_msg->table_ver == 128) {
	    *discipline=10;*parameter_category=4,*parameter_number=5;
	  }
	  else if (grib_msg->table_ver == 130) {
	    *discipline=2;*parameter_category=0;*parameter_number=7;
	  }
	}
	return;
    }
    case 197:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=2;*parameter_number=28;
	  }
	  else if (grib_msg->table_ver == 128) {
	    *discipline=10;*parameter_category=4,*parameter_number=6;
	  }
	}
	return;
    }
    case 200:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=2;*parameter_category=0;*parameter_number=13;
	}
	return;
    }
    case 202:
    {
	if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 133) {
	    *discipline=0;*parameter_category=3;*parameter_number=27;
	  }
	  else if (grib_msg->table_ver == 200) {
	    *discipline=2;*parameter_category=0;*parameter_number=6;
	  }
	}
	return;
    }
    case 203:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=0;*parameter_number=16;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=0;*parameter_number=18;
	}
	return;
    }
    case 204:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=4,*parameter_number=7;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 200) {
	  *discipline=0;*parameter_category=4,*parameter_number=7;
	}
	return;
    }
    case 205:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=5;*parameter_number=3;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 200) {
	  *discipline=0;*parameter_category=5;*parameter_number=3;
	}
	return;
    }
    case 206:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=0;*parameter_category=15;*parameter_number=3;
	}
	return;
    }
    case 207:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=2;*parameter_category=0;*parameter_number=11;
	}
	return;
    }
    case 208:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=2;*parameter_category=0;*parameter_number=12;
	}
	return;
    }
    case 209:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 133) {
	  *discipline=0;*parameter_category=3;*parameter_number=28;
	}
	return;
    }
    case 211:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=4,*parameter_number=8;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 200) {
	  *discipline=0;*parameter_category=4,*parameter_number=8;
	}
	return;
    }
    case 212:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=5;*parameter_number=4;
	}
	else if (grib_msg->center_id == 98) {
	  if (grib_msg->table_ver == 200) {
	    *discipline=0;*parameter_category=5;*parameter_number=4;
	  }
	  else if (grib_msg->table_ver == 201) {
	    *discipline=2;*parameter_category=0;*parameter_number=16;
	  }
	}
	return;
    }
    case 214:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=37;
	}
	return;
    }
    case 218:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 129) {
	  *discipline=0;*parameter_category=1;*parameter_number=27;
	}
	return;
    }
    case 219:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=6;*parameter_number=13;
	  }
	  else if (grib_msg->table_ver == 130) {
	    *discipline=2;*parameter_category=0;*parameter_number=17;
	  }
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 200) {
	  *discipline=0;*parameter_category=2;*parameter_number=21;
	}
	return;
    }
    case 221:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=3;*parameter_number=18;
	}
	return;
    }
    case 222:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=3;*parameter_number=15;
	}
	return;
    }
    case 223:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=2;*parameter_category=0;*parameter_number=13;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=65;
	  }
	}
	return;
    }
    case 224:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=2;*parameter_category=3;*parameter_number=0;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=66;
	  }
	}
	return;
    }
    case 225:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 129) {
	  *discipline=0;*parameter_category=1;*parameter_number=67;
	}
	return;
    }
    case 226:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=2;*parameter_category=0;*parameter_number=14;
	  }
	  else if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=1;*parameter_number=68;
	  }
	}
	return;
    }
    case 227:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 129) {
	  *discipline=0;*parameter_category=7;*parameter_number=15;
	}
	return;
    }
    case 228:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=40;
	}
	return;
    }
    case 229:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=0;*parameter_number=16;
	}
	return;
    }
    case 230:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 2) {
	    *discipline=0;*parameter_category=3;*parameter_number=19;
	  }
	  else if (grib_msg->table_ver == 130) {
	    *discipline=2;*parameter_category=3;*parameter_number=7;
	  }
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 201) {
	  *discipline=0;*parameter_category=15;*parameter_number=1;
	}
	return;
    }
    case 231:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=3;*parameter_number=8;
	}
	return;
    }
    case 233:
    {
	if (grib_msg->center_id == 98 && grib_msg->table_ver == 140) {
	  *discipline=10;*parameter_category=0;*parameter_number=16;
	}
	return;
    }
    case 234:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=1;*parameter_category=0;*parameter_number=5;
	}
	return;
    }
    case 235:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=1;*parameter_category=0;*parameter_number=6;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=0;*parameter_number=17;
	}
	return;
    }
    case 238:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=1;*parameter_number=42;
	}
	return;
    }
    case 239:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=0;*parameter_number=18;
	}
	return;
    }
    case 240:
    {
	if (grib_msg->center_id == 7) {
	  if (grib_msg->table_ver == 129) {
	    *discipline=0;*parameter_category=16;*parameter_number=3;
	  }
	  else if (grib_msg->table_ver == 130) {
	    *discipline=2;*parameter_category=3;*parameter_number=9;
	  }
	}
	return;
    }
    case 246:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=0;*parameter_number=18;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=1;*parameter_number=83;
	}
	return;
    }
    case 247:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=0;*parameter_number=19;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=1;*parameter_number=84;
	}
	return;
    }
    case 248:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=0;*parameter_number=20;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 128) {
	  *discipline=0;*parameter_category=6;*parameter_number=32;
	}
	return;
    }
    case 249:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 130) {
	  *discipline=2;*parameter_category=0;*parameter_number=21;
	}
	return;
    }
    case 250:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 129) {
	  *discipline=2;*parameter_category=4,*parameter_number=2;
	}
	return;
    }
    case 252:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=2;*parameter_number=29;
	}
	return;
    }
    case 253:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=2;*parameter_number=30;
	}
	else if (grib_msg->center_id == 98 && grib_msg->table_ver == 140) {
	  *discipline=10;*parameter_category=0;*parameter_number=44;
	}
	return;
    }
    case 254:
    {
	if (grib_msg->center_id == 7 && grib_msg->table_ver == 2) {
	  *discipline=0;*parameter_category=7;*parameter_number=12;
	}
	return;
    }
  }
}

void pack_IDS(GRIBMessage *msg,unsigned char *grib2_buffer,size_t *offset)
{
// length of the IDS
  size_t length=21;
  set_bits(grib2_buffer,length,*offset,32);
// section number
  set_bits(grib2_buffer,1,*offset+32,8);
// originating center code
  set_bits(grib2_buffer,msg->center_id,*offset+40,16);
// sub-center code
  set_bits(grib2_buffer,msg->sub_center_id,*offset+56,16);
// master table version number
  set_bits(grib2_buffer,18,*offset+72,8);
// local table version number
  set_bits(grib2_buffer,0,*offset+80,8);
// significance of reference time
  set_bits(grib2_buffer,1,*offset+88,8);
// year, month, day, hour, minute, second of reference time
  set_bits(grib2_buffer,msg->yr,*offset+96,16);
  set_bits(grib2_buffer,msg->mo,*offset+112,8);
  set_bits(grib2_buffer,msg->dy,*offset+120,8);
  set_bits(grib2_buffer,msg->time/100,*offset+128,8);
  set_bits(grib2_buffer,(msg->time % 100),*offset+136,8);
  set_bits(grib2_buffer,0,*offset+144,8);
// production status
  set_bits(grib2_buffer,255,*offset+152,8);
// type of data
  set_bits(grib2_buffer,255,*offset+160,8);
  (*offset)+=length*8;
}

void pack_GDS(GRIBMessage *msg,unsigned char *grib2_buffer,size_t *offset)
{
  size_t length,template_num;
  switch (msg->data_rep) {
    case 0:
    {
// Latitude-Longitude
	length=72;
	template_num=0;
	break;
    }
    case 1:
    {
// Mercator
	length=72;
	template_num=10;
	break;
    }
    case 4:
    {
// Gaussian Latitude-Longitude
	length=72;
	template_num=40;
	break;
    }
    case 5:
    {
// Polar Stereographic
	length=65;
	template_num=20;
	break;
    }
    default:
    {
	fprintf(stderr,"Unable to convert grid %d\n",msg->data_rep);
	exit(1);
    }
  }
// length of the GDS
  set_bits(grib2_buffer,length,*offset,32);
// section number
  set_bits(grib2_buffer,3,*offset+32,8);
// source of grid definition
  set_bits(grib2_buffer,0,*offset+40,8);
// number of data points
  set_bits(grib2_buffer,(msg->nx*msg->ny),*offset+48,32);
// octets 11 and 12
  set_bits(grib2_buffer,0,*offset+80,16);
// template number
  set_bits(grib2_buffer,template_num,*offset+96,16);
  switch (template_num) {
    case 0:
    {
// Latitude-Longitude
// shape of the Earth
	set_bits(grib2_buffer,6,*offset+112,8);
// next six fields relate to shape of earth, which are not used because we used
//  a standard radius in octet 15
	set_bits(grib2_buffer,0,*offset+120,8);
	set_bits(grib2_buffer,0,*offset+128,32);
	set_bits(grib2_buffer,0,*offset+160,8);
	set_bits(grib2_buffer,0,*offset+168,32);
	set_bits(grib2_buffer,0,*offset+200,8);
	set_bits(grib2_buffer,0,*offset+208,32);
// number of points in i-direction
	set_bits(grib2_buffer,msg->nx,*offset+240,32);
// number of points in j-direction
	set_bits(grib2_buffer,msg->ny,*offset+272,32);
// octets 39-42 and 43-46 can be ignored because we don't have high-precision
//  latitudes and longitudes
	set_bits(grib2_buffer,0,*offset+304,32);
	set_bits(grib2_buffer,0,*offset+336,32);
// latitude of first gridpoint
	if (msg->slat < 0.) {
	  set_bits(grib2_buffer,1,*offset+368,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+368,1);
	}
	set_bits(grib2_buffer,abs(msg->slat*1000000.),*offset+369,31);
// longitude of first gridpoint
	if (msg->slon < 0.) {
	  set_bits(grib2_buffer,1,*offset+400,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+400,1);
	}
	set_bits(grib2_buffer,abs(msg->slon*1000000.),*offset+401,31);
// resolution and component flags
	unsigned char rcflg=((msg->rescomp & 0x80) >> 2) | ((msg->rescomp & 0x80) >> 3) | (msg->rescomp & 0xf);
	set_bits(grib2_buffer,rcflg,*offset+432,8);
// latitude of last gridpoint
	if (msg->elat < 0.) {
	  set_bits(grib2_buffer,1,*offset+440,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+440,1);
	}
	set_bits(grib2_buffer,abs(msg->elat*1000000.),*offset+441,31);
// longitude of last gridpoint
	if (msg->elon < 0.) {
	  set_bits(grib2_buffer,1,*offset+472,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+472,1);
	}
	set_bits(grib2_buffer,abs(msg->elon*1000000.),*offset+473,31);
// i-direction increment
	set_bits(grib2_buffer,lround(msg->loinc*1000000.),*offset+504,32);
// j-direction increment
	set_bits(grib2_buffer,lround(msg->lainc*1000000.),*offset+536,32);
// scanning mode
	set_bits(grib2_buffer,msg->scan_mode,*offset+568,8);
	break;
    }
    case 10:
    {
// Mercator
// shape of the Earth
	set_bits(grib2_buffer,6,*offset+112,8);
// next six fields relate to shape of earth, which are not used because we used
//  a standard radius in octet 15
	set_bits(grib2_buffer,0,*offset+120,8);
	set_bits(grib2_buffer,0,*offset+128,32);
	set_bits(grib2_buffer,0,*offset+160,8);
	set_bits(grib2_buffer,0,*offset+168,32);
	set_bits(grib2_buffer,0,*offset+200,8);
	set_bits(grib2_buffer,0,*offset+208,32);
// number of points in i-direction
	set_bits(grib2_buffer,msg->nx,*offset+240,32);
// number of points in j-direction
	set_bits(grib2_buffer,msg->ny,*offset+272,32);
// latitude of first gridpoint
	if (msg->slat < 0.) {
	  set_bits(grib2_buffer,1,*offset+304,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+304,1);
	}
	set_bits(grib2_buffer,abs(msg->slat*1000000.),*offset+305,31);
// longitude of first gridpoint
	if (msg->slon < 0.) {
	  set_bits(grib2_buffer,1,*offset+336,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+336,1);
	}
	set_bits(grib2_buffer,abs(msg->slon*1000000.),*offset+337,31);
// resolution and component flags
	unsigned char rcflg=((msg->rescomp & 0x80) >> 2) | ((msg->rescomp & 0x80) >> 3) | (msg->rescomp & 0xf);
	set_bits(grib2_buffer,rcflg,*offset+368,8);
// standard latitude (latitude where projection intersects the Earth)
	if (msg->std_lat1 < 0.) {
	  set_bits(grib2_buffer,1,*offset+376,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+376,1);
	}
	set_bits(grib2_buffer,abs(msg->std_lat1*1000000.),*offset+377,31);
// latitude of last gridpoint
	if (msg->elat < 0.) {
	  set_bits(grib2_buffer,1,*offset+408,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+408,1);
	}
	set_bits(grib2_buffer,abs(msg->elat*1000000.),*offset+409,31);
// longitude of last gridpoint
	if (msg->elon < 0.) {
	  set_bits(grib2_buffer,1,*offset+440,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+440,1);
	}
	set_bits(grib2_buffer,abs(msg->elon*1000000.),*offset+441,31);
// scanning mode
	set_bits(grib2_buffer,msg->scan_mode,*offset+472,8);
// orientation of the grid
	set_bits(grib2_buffer,0,*offset+480,32);
// i-direction increment
	set_bits(grib2_buffer,lround(msg->xlen*1000.),*offset+512,32);
// j-direction increment
	set_bits(grib2_buffer,lround(msg->ylen*1000.),*offset+544,32);
	break;
    }
    case 20:
    {
// Polar Stereographic
// shape of the Earth
	set_bits(grib2_buffer,6,*offset+112,8);
// next six fields relate to shape of earth, which are not used because we used
//  a standard radius in octet 15
	set_bits(grib2_buffer,0,*offset+120,8);
	set_bits(grib2_buffer,0,*offset+128,32);
	set_bits(grib2_buffer,0,*offset+160,8);
	set_bits(grib2_buffer,0,*offset+168,32);
	set_bits(grib2_buffer,0,*offset+200,8);
	set_bits(grib2_buffer,0,*offset+208,32);
// number of points in i-direction
	set_bits(grib2_buffer,msg->nx,*offset+240,32);
// number of points in j-direction
	set_bits(grib2_buffer,msg->ny,*offset+272,32);
// latitude of first gridpoint
	if (msg->slat < 0.) {
	  set_bits(grib2_buffer,1,*offset+304,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+304,1);
	}
	set_bits(grib2_buffer,abs(msg->slat*1000000.),*offset+305,31);
// longitude of first gridpoint
	if (msg->slon < 0.) {
	  set_bits(grib2_buffer,1,*offset+336,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+336,1);
	}
	set_bits(grib2_buffer,abs(msg->slon*1000000.),*offset+337,31);
// resolution and component flags
	unsigned char rcflg=((msg->rescomp & 0x80) >> 2) | ((msg->rescomp & 0x80) >> 3) | (msg->rescomp & 0xf);
	set_bits(grib2_buffer,rcflg,*offset+368,8);
// latitude at which dx and dy are valid
	if (msg->proj == 1) {
	  set_bits(grib2_buffer,1,*offset+376,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+376,1);
	}
	set_bits(grib2_buffer,40000000,*offset+377,31);
// longitude of orientation
	if (msg->olon < 0.) {
	  set_bits(grib2_buffer,1,*offset+408,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+408,1);
	}
	set_bits(grib2_buffer,lround(fabs(msg->olon)*1000000.),*offset+409,31);
// dx
	set_bits(grib2_buffer,msg->xlen*1000,*offset+440,32);
// dy
	set_bits(grib2_buffer,msg->ylen*1000,*offset+472,32);
// projection center flag
	set_bits(grib2_buffer,msg->proj,*offset+504,8);
// scanning mode
	unsigned char scan_mode=(msg->scan_mode | 0x10);
	set_bits(grib2_buffer,scan_mode,*offset+512,8);
	break;
    }
    case 40:
    {
// Gaussian Latitude-Longitude
// shape of the Earth
	set_bits(grib2_buffer,6,*offset+112,8);
// next six fields relate to shape of earth, which are not used because we used
//  a standard radius in octet 15
	set_bits(grib2_buffer,0,*offset+120,8);
	set_bits(grib2_buffer,0,*offset+128,32);
	set_bits(grib2_buffer,0,*offset+160,8);
	set_bits(grib2_buffer,0,*offset+168,32);
	set_bits(grib2_buffer,0,*offset+200,8);
	set_bits(grib2_buffer,0,*offset+208,32);
// number of points in i-direction
	set_bits(grib2_buffer,msg->nx,*offset+240,32);
// number of points in j-direction
	set_bits(grib2_buffer,msg->ny,*offset+272,32);
// octets 39-42 and 43-46 can be ignored because we don't have high-precision
//  latitudes and longitudes
	set_bits(grib2_buffer,0,*offset+304,32);
	set_bits(grib2_buffer,0,*offset+336,32);
// latitude of first gridpoint
	if (msg->slat < 0.) {
	  set_bits(grib2_buffer,1,*offset+368,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+368,1);
	}
	set_bits(grib2_buffer,abs(msg->slat*1000000.),*offset+369,31);
// longitude of first gridpoint
	if (msg->slon < 0.) {
	  set_bits(grib2_buffer,1,*offset+400,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+400,1);
	}
	set_bits(grib2_buffer,abs(msg->slon*1000000.),*offset+401,31);
// resolution and component flags
	unsigned char rcflg=((msg->rescomp & 0x80) >> 2) | ((msg->rescomp & 0x80) >> 3) | (msg->rescomp & 0xf);
	set_bits(grib2_buffer,rcflg,*offset+432,8);
// latitude of last gridpoint
	if (msg->elat < 0.) {
	  set_bits(grib2_buffer,1,*offset+440,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+440,1);
	}
	set_bits(grib2_buffer,abs(msg->elat*1000000.),*offset+441,31);
// longitude of last gridpoint
	if (msg->elon < 0.) {
	  set_bits(grib2_buffer,1,*offset+472,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+472,1);
	}
	set_bits(grib2_buffer,abs(msg->elon*1000000.),*offset+473,31);
// i-direction increment
	set_bits(grib2_buffer,lround(msg->lainc*1000000.),*offset+504,32);
// number of parallels between equator and pole
	set_bits(grib2_buffer,lround(msg->lainc),*offset+536,32);
// scanning mode
	unsigned char scan_mode=(msg->scan_mode | 0x10);
	set_bits(grib2_buffer,scan_mode,*offset+568,8);
	break;
    }
  }
  (*offset)+=length*8;
}

int mdays[]={0,31,28,31,30,31,30,31,31,30,31,30,31};
void add_time(int time_to_add,int time_units,int *yr,int *mo,int *dy,int *time)
{
  int hr=*time/100;
  int min=(*time % 100);
  switch (time_units) {
    case 0:
    {
	min+=time_to_add;
	break;
    }
    case 1:
    {
	min+=time_to_add*60;
	break;
    }
    case 2:
    {
	min+=time_to_add*1440;
	break;
    }
    default:
    {
	fprintf(stderr,"Unable to add time for unit %d\n",time_units);
	exit(1);
    }
  }
  if (min >= 60) {
    hr+=min/60;
    min=(min % 60);
    if (hr >= 24) {
	*dy+=hr/24;
	hr=(hr % 24);
	if ( (*yr % 4) == 0 && ( (*yr % 100 != 0) || (*yr % 400) == 0)) {
	  mdays[2]=29;
	}
	else {
	  mdays[2]=28;
	}
	while (*dy > mdays[*mo]) {
	  *dy-=mdays[*mo];
	  ++(*mo);
	  if (*mo > 12) {
	    ++(*yr);
	    *mo=1;
	    if ( (*yr % 4) == 0 && ( (*yr % 100 != 0) || (*yr % 400) == 0)) {
		mdays[2]=29;
	    }
	    else {
		mdays[2]=28;
	    }
	  }
	}
    } 
  }
  *time=hr*100+min;
}

void pack_PDS(GRIBMessage *msg,int parameter_category,int parameter_number,unsigned char *grib2_buffer,size_t *offset)
{
  size_t length,template_num;
  switch (msg->t_range) {
    case 0:
    case 1:
    case 10:
    {
	length=34;
	template_num=0;
	break;
    }
    case 2:
    case 3:
    case 4:
    {
	length=58;
	template_num=8;
	break;
    }
    default:
    {
	fprintf(stderr,"Unable to convert time range %d\n",msg->t_range);
	exit(1);
    }
  }
// length of the PDS
  set_bits(grib2_buffer,length,*offset,32);
// section number
  set_bits(grib2_buffer,4,*offset+32,8);
// number of coordinates after template
  set_bits(grib2_buffer,0,*offset+40,16);
// template number
  set_bits(grib2_buffer,template_num,*offset+56,16);
  switch (template_num) {
    case 0:
    case 8:
    {
// parameter category
	set_bits(grib2_buffer,parameter_category,*offset+72,8);
// parameter number
	set_bits(grib2_buffer,parameter_number,*offset+80,8);
// generating process
	set_bits(grib2_buffer,255,*offset+88,8);
// background generating process
	set_bits(grib2_buffer,msg->gen_proc,*offset+96,8);
// analysis or forecast generating process
	set_bits(grib2_buffer,255,*offset+104,8);
// hours of observational data cutoff
	set_bits(grib2_buffer,65535,*offset+112,16);
// minutes of observational data cutoff
	set_bits(grib2_buffer,255,*offset+128,8);
// time range units
	set_bits(grib2_buffer,msg->fcst_units,*offset+136,8);
// forecast time
	switch (msg->t_range) {
	  case 0:
	  case 10:
	  {
	    set_bits(grib2_buffer,msg->p1,*offset+144,32);
	    break;
	  }
	  case 1:
	  {
	    set_bits(grib2_buffer,0,*offset+144,32);
	    break;
	  }
	}
// type of first surface
// scale factor of first surface
// scaled value of first surface
// type of second surface
// scale factor of second surface
// scaled value of second surface
	int lvl1_type=msg->level_type,lvl2_type=255;
	int lvl1_scale=0,lvl2_scale=255;
	int lvl1_value=msg->lvl1,lvl2_value=msg->lvl2;
	switch (msg->level_type) {
	  case 20:
	  case 100:
	  {
	    lvl1_scale=-2;
	    break;
	  }
	  case 101:
	  {
	    lvl1_type=lvl2_type=100;
	    lvl1_scale=lvl2_scale=-3;
	    break;
	  }
	  case 102:
	  {
	    lvl1_type=101;
	    break;
	  }
	  case 103:
	  {
	    lvl1_type=102;
	    break;
	  }
	  case 104:
	  {
	    lvl1_type=lvl2_type=102;
	    lvl1_scale=lvl2_scale=-2;
	  }
	  case 105:
	  {
	    lvl1_type=103;
	    break;
	  }
	  case 106:
	  {
	    lvl1_type=lvl2_type=103;
	    lvl1_scale=lvl2_scale=-2;
	    break;
	  }
	  case 107:
	  {
	    lvl1_type=104;
	    lvl1_scale=4;
	    break;
	  }
	  case 108:
	  {
	    lvl1_type=lvl2_type=104;
	    lvl1_scale=lvl2_scale=2;
	    break;
	  }
	  case 109:
	  {
	    lvl1_type=105;
	    break;
	  }
	  case 110:
	  {
	    lvl1_type=lvl2_type=105;
	    break;
	  }
	  case 111:
	  {
	    lvl1_type=106;
	    lvl1_scale=2;
	    break;
	  }
	  case 112:
	  {
	    lvl1_type=lvl2_type=106;
	    lvl1_scale=lvl2_scale=2;
	    break;
	  }
	  case 113:
	  {
	    lvl1_type=107;
	    break;
	  }
	  case 114:
	  {
	    lvl1_type=lvl2_type=107;
	    lvl1_value=475-lvl1_value;
	    lvl2_value=475-lvl2_value;
	    break;
	  }
	  case 115:
	  {
	    lvl1_type=108;
	    lvl1_scale=-2;
	    break;
	  }
	  case 116:
	  {
	    lvl1_type=lvl2_type=108;
	    lvl1_scale=lvl2_scale=-2;
	    break;
	  }
	  case 117:
	  {
	    lvl1_type=109;
	    lvl1_scale=9;
	    break;
	  }
	  case 119:
	  {
	    lvl1_type=111;
	    lvl1_scale=4;
	    break;
	  }
	  case 120:
	  {
	    lvl1_type=lvl2_type=111;
	    lvl1_scale=lvl2_scale=2;
	    break;
	  }
	  case 121:
	  {
	    lvl1_type=lvl2_type=100;
	    lvl1_scale=lvl2_scale=-2;
	    lvl1_value=1100-lvl1_value;
	    lvl2_value=1100-lvl2_value;
	    break;
	  }
	  case 125:
	  {
	    lvl1_type=103;
	    lvl1_scale=2;
	    break;
	  }
	  case 128:
	  {
	    lvl1_type=lvl2_type=104;
	    lvl1_scale=lvl2_scale=3;
	    lvl1_value=1100-lvl1_value;
	    lvl2_value=1100-lvl2_value;
	    break;
	  }
	  case 141:
	  {
	    lvl1_type=lvl2_type=100;
	    lvl1_scale=-3;
	    lvl2_scale=-2;
	    lvl2_value=1100-lvl2_value;
	    break;
	  }
	}
	set_bits(grib2_buffer,lvl1_type,*offset+176,8);
	if (lvl1_scale < 0) {
	  set_bits(grib2_buffer,1,*offset+184,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+184,1);
	}
	set_bits(grib2_buffer,abs(lvl1_scale),*offset+185,7);
	set_bits(grib2_buffer,lvl1_value,*offset+192,32);
	set_bits(grib2_buffer,lvl2_type,*offset+224,8);
	if (lvl2_scale < 0) {
	  set_bits(grib2_buffer,1,*offset+232,1);
	}
	else {
	  set_bits(grib2_buffer,0,*offset+232,1);
	}
	set_bits(grib2_buffer,abs(lvl2_scale),*offset+233,7);
	set_bits(grib2_buffer,lvl2_value,*offset+240,32);
	if (template_num == 8) {
	  int yr=msg->yr,mo=msg->mo,dy=msg->dy,time=msg->time;
	  add_time(msg->p2,msg->fcst_units,&yr,&mo,&dy,&time);
// year of overall end time
	  set_bits(grib2_buffer,yr,*offset+272,16);
// month of overall end time
	  set_bits(grib2_buffer,mo,*offset+288,8);
// day of overall end time
	  set_bits(grib2_buffer,dy,*offset+296,8);
// hour of overall end time
	  set_bits(grib2_buffer,time/100,*offset+304,8);
// minute of overall end time
	  set_bits(grib2_buffer,(time % 100),*offset+312,8);
// second of overall end time
	  set_bits(grib2_buffer,0,*offset+320,8);
// number of time range specifications
	  set_bits(grib2_buffer,1,*offset+328,8);
// number missing from statistical process
	  set_bits(grib2_buffer,0,*offset+336,32);
// statistical process
	  int process=255,time_incr;
	  switch (msg->t_range) {
	    case 4:
	    {
		process=1;
		time_incr=2;
		break;
	    }
	  }
	  if (process == 255) {
	    switch (msg->param) {
		case 15:
		{
// Maximum temperature
		  process=2;
		  time_incr=2;
		  break;
		}
		case 16:
		{
// Minimum temperature
		  process=3;
		  time_incr=2;
		  break;
		}
	    }
	  }
	  if (process == 255) {
	    fprintf(stderr,"Unable to determine statistical process type for parameter code %d\n",msg->param);
	    exit(1);
	  }
	  else {
	    set_bits(grib2_buffer,process,*offset+368,8);
	  }
// type of time increment between successive fields used in statistic
	  set_bits(grib2_buffer,time_incr,*offset+376,8);
// indicator of time unit for statistical process
	  set_bits(grib2_buffer,msg->fcst_units,*offset+384,8);
// length of time range of statistical process
	  set_bits(grib2_buffer,(msg->p2-msg->p1),*offset+392,32);
// indicator of time unit of increment between successive fields
	  set_bits(grib2_buffer,msg->fcst_units,*offset+424,8);
// time increment between successive fields
	  set_bits(grib2_buffer,0,*offset+432,32);
	}
	break;
    }
  }
  (*offset)+=length*8;
}

void pack_DRS(GRIBMessage *msg,unsigned char *grib2_buffer,size_t *offset)
{
// length of the DRS
  size_t length=21;
  set_bits(grib2_buffer,length,*offset,32);
// section number
  set_bits(grib2_buffer,5,*offset+32,8);
// number of data points
  set_bits(grib2_buffer,msg->nx*msg->ny,*offset+40,32);
// template number
  set_bits(grib2_buffer,0,*offset+72,16);
// reference value
  union {
    float f;
    int i;
  } u;
  u.f=msg->ref_val*pow(10.,msg->D);
  set_bits(grib2_buffer,u.i,*offset+88,32);
// binary scale factor
  float E=msg->E;
  if (E < 0) {
    E=0x8000-E;
  }
  set_bits(grib2_buffer,E,*offset+120,16);
// decimal scale factor
  float D=msg->D;
  if (D < 0) {
    D=0x8000-D;
  }
  set_bits(grib2_buffer,D,*offset+136,16);
// packed value width
  set_bits(grib2_buffer,msg->pack_width,*offset+152,8);
// type of original values
  set_bits(grib2_buffer,0,*offset+160,8);
  (*offset)+=length*8;
}

void pack_BMS(GRIBMessage *msg,unsigned char *grib2_buffer,size_t *offset)
{
// length of the BMS
  size_t length=6;
  if (msg->bms_included) {
    length+=(msg->bitmap_len+7)/8;
  }
  set_bits(grib2_buffer,length,*offset,32);
// section number
  set_bits(grib2_buffer,6,*offset+32,8);
  if (length > 6) {
// bitmap indicator
    set_bits(grib2_buffer,0,*offset+40,8);
// bitmap
    size_t off=*offset+48;
    for (size_t n=0; n < msg->bitmap_len; ++n) {
	int bval=msg->bitmap[n];
	set_bits(grib2_buffer,bval,off,1);
	++off;
    }
  }
  else {
// bitmap indicator
    set_bits(grib2_buffer,255,*offset+40,8);
  }
  (*offset)+=length*8;
}

void pack_DS(GRIBMessage *msg,unsigned char *grib2_buffer,size_t *offset)
{
// length of the DS
  size_t length=5+(msg->nx*msg->ny*msg->pack_width+7)/8;
  set_bits(grib2_buffer,length,*offset,32);
// section number
  set_bits(grib2_buffer,7,*offset+32,8);
  size_t off=*offset+40;
  float d=pow(10.,msg->D);
  float e=pow(2.,msg->E);
  for (size_t n=0; n < msg->nx*msg->ny; ++n) {
    if (msg->gridpoints[n] != GRIB_MISSING_VALUE) {
	int pval=lround((msg->gridpoints[n]-msg->ref_val)*d/e);
	set_bits(grib2_buffer,pval,off,msg->pack_width);
	off+=msg->pack_width;
    }
  }
  (*offset)+=length*8;
}

int main(int argc,char **argv)
{
  if (argc != 3) {
    fprintf(stderr,"usage: %s GRIB1_file_name GRIB2_file_name\n",argv[0]);
    exit(1);
  }
  FILE *ifile;
  if ( (ifile=fopen(argv[1],"rb")) == NULL) {
    fprintf(stderr,"Error opening input file %s\n",argv[1]);
    exit(1);
  }
  FILE *ofile;
  if ( (ofile=fopen(argv[2],"wb")) == NULL) {
    fprintf(stderr,"Error opening output file %s\n",argv[2]);
    exit(1);
  }
  GRIBMessage grib_msg;
  initialize(&grib_msg);
  int status;
  size_t nmsg=0;
  size_t max_buffer_length=0;
  unsigned char *grib2_buffer=NULL;
  char head[]={'G','R','I','B',0,0,0,2,0,0,0,0,0,0,0,0};
  char tail[]={'7','7','7','7'};
  long long length;
  while ( (status=unpackgrib1(ifile,&grib_msg)) == 0) {
    ++nmsg;
// Identification Section
    length=21;
// Grid Definition Section
    if (grib_msg.data_rep == 3) {
	length+=81;
    }
    else {
	length+=72;
    }
// Product Definition Section
    if (grib_msg.t_range <= 1 || grib_msg.t_range == 10) {
	length+=34;
    }
    else if (grib_msg.t_range >= 2 && grib_msg.t_range <= 4) {
	length+=58;
    }
    else {
	fprintf(stderr,"Unable to convert time range indicator %d\n",grib_msg.t_range);
	exit(1);
    }
// Data Representation Section
    length+=21;
// Bit-map Section
    length+=6;
    if (grib_msg.bms_included) {
	length+=(grib_msg.bitmap_len+7)/8;
    }
// Data Section
    length+=5+(grib_msg.nx*grib_msg.ny*grib_msg.pack_width+7)/8;
// allocate enough memory for the GRIB2 buffer
    if (length > max_buffer_length) {
	if (grib2_buffer != NULL) {
	  free(grib2_buffer);
	}
	max_buffer_length=length;
	grib2_buffer=(unsigned char *)malloc(max_buffer_length*sizeof(unsigned char));
    }
    int discipline,parameter_category,parameter_number;
    map_GRIB1_parameter_into_GRIB2(&grib_msg,&discipline,&parameter_category,&parameter_number);
    size_t offset=0;
// pack the Identification Section
    pack_IDS(&grib_msg,grib2_buffer,&offset);
// pack the Grid Definition Section
    pack_GDS(&grib_msg,grib2_buffer,&offset);
// pack the Product Definition Section
    pack_PDS(&grib_msg,parameter_category,parameter_number,grib2_buffer,&offset);
// pack the Data Representation Section
    pack_DRS(&grib_msg,grib2_buffer,&offset);
// pack the Bit Map Section
    pack_BMS(&grib_msg,grib2_buffer,&offset);
// pack the Data Section
    pack_DS(&grib_msg,grib2_buffer,&offset);
// output the GRIB2 message
    size_t l=length+20;
    set_bits(&head[8],(l >> 32),0,32);
    set_bits(&head[12],(l & 0xffffffff),0,32);
    fwrite(head,1,16,ofile);
    fwrite(grib2_buffer,1,length,ofile);
    fwrite(tail,1,4,ofile);
  }
  if (status != -1) {
    fprintf(stderr,"Read error after %d messages\n",nmsg);
  }
  printf("Number of GRIB2 messages written to output: %d\n",nmsg);
  fclose(ifile);
  fclose(ofile);
}
