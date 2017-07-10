/*
** File: unpackgrib1.c
**
** Author:  Bob Dattore
**          NCAR/DSS
**          dattore@ucar.edu
**          (303) 497-1825
**
** Revision History:
**   21 Feb 2001 - modified unpack_IS to search for next GRIB message in files
**                 where the GRIB messages are not contiguous
**   22 Feb 2001 - modified routine to handle grids with data representation of
**                 10
**   26 Feb 2001 - modified to run under Linux - added byteSwap routine to
**                 handle byte-swapping
**   15 Feb 2008 - changed GRIB buffer stream handling to byte rather than word
**                 to remove the need for byte-swapping on little-endian
**                 machines, removed byteSwap routine
**   21 Nov 2008 - bug fix
**   18 Jun 2010 - code fix to handle constant fields
**   23 Oct 2012 - fixed a memory leak in unpack_BDS
**   20 May 2017 - code refactoring, added bitmap to the GRIBMessage structure
**                 to reduce memory reallocations
**   10 Jul 2017 - decode Mercator grid definition
**
** Purpose: to provide a single C-routine for unpacking GRIB grids
**
** Notes:   1) There are several routines defined in this file, any of which can
**             be called independently. However, to unpack a GRIB record with
**             one function call, use only "unpackgrib1".
**
**          2) The user is expected to have some understanding of the GRIB
**             format, as some of the information consists of codes, etc. that
**             are defined by the GRIB format. You can find a description of
**             the GRIB format at
**             https://rda.ucar.edu/docs/formats/grib/gribdoc/
**
**          3) As of this version, these routines unpack the following types of
**             grids:
**               Latitude/Longitude (data representation = 0)
**               Gaussian Latitude/Longitude (data representation = 4)
**               Polar Stereographic (data representation = 5)
**               Rotated Latitude/Longitude (data representation = 10)
**             please contact dattore@ucar.edu to get a grid type added
**
**
** example C syntax for using unpackgrib1:
**    FILE *fp;
**    GRIBMessage grib_msg;
**    int status;
**
**    initialize(&grib_msg);
**    fp=fopen("my_GRIB_file","rb");
**    while ( (status=unpackgrib1(fp,&grib_msg)) == 0) {
**      if (status == -1) {
**        printf("Found EOF\n");
**      }
**      else if (status == 1) {
**        printf("Error reading GRIB record\n");
**      }
**      ...
**    }
** 
** where:
**   fp            is a FILE pointer to an open GRIB data file
**   GRIBMessage   is a structure used to hold the data in a GRIB message
**
** On return:
**   unpackgrib1 returns 0 for a successful read, -1 for an EOF, and 1 for a
**   read error
**
** 
** Overview of GRIBMessage:
**   total_len:     Total length of the GRIB record, in octets (8-bit bytes)
**   pds_len:       Length in octets of the Product Definition Section (PDS)
**   pds_ext_len:   Length in octets of any extension to the PDS
**   gds_len:       Length in octets of the Grid Definition Section (GDS)
**   bds_len:       Length in octets of the Binary Data Section (BDS)
**   ed_num:        GRIB Edition number
**   table_ver:     GRIB Parameter Table version
**   center_id:     Center ID
**   gen_proc:      Generating Process ID number
**   grid_type:     Grid Identification
**   param:         GRIB Parameter code
**   level_type:    Level Type Code
**   lvl1:          Value for First Level
**   lvl2:          Value for Second Level (or 0)
**   fcst_units:    Forecast Time Unit
**   p1:            P1 (or 0)
**   p2:            P2 (or 0)
**   t_range:       Time Range Indicator
**   navg:          Number included in Average (or 0)
**   nmiss:         Number of Missing Grids in Average
**   sub_center_id: Sub-center ID
**   bds_flag:      Binary Data Section flag
**   pack_width:    Width in bits of a packed data point
**   gds_included:  Indication of GDS (0 = no GDS; 1 = GDS in GRIB record)
**   bms_included:  Indication of Bitmap Section (BMS)
**                    (0 = no BMS; 1 = BMS in GRIB record)
**   yr:            Year (4-digits - YYYY)
**   mo:            Month
**   dy:            Day
**   time:          Time (HHMM - HH=hour, MM=minutes)
**   offset:        For Internal Use Only (offset in bytes to next GRIB section
**                    from the beginning of current section)
**   E:             Binary scale factor
**   D:             Decimal scale factor
**   data_rep:      Data representation type
** 
**   For Latitude/Longitude and Gaussian Lat/Lon Grids:
**     nx:         Number of points along a latitude circle
**     ny:         Number of points along a longitude meridian
**     slat:       Latitude of the first gridpoint (*1000)
**     slon:       Longitude of the first gridpoint (*1000)
**     rescomp:    Resolution and component flags
**     elat:       Latitude of the last gridpoint (*1000)
**     elon:       Longitude of the last gridpoint (*1000)
**     lainc:      Latitude increment (*1000) for Lat/Lon grid
**                 -OR-
**                 Number of latitude circles between equator and pole for
**                   Gaussian Lat/Lon grid
**     loinc:      Longitude increment (*1000)
**     scan_mode:  Scanning mode flags
**
**   For Polar Stereographic Grids:
**     nx:         Number of points in the X-direction
**     ny:         Number of points in the Y-direction
**     slat:       Latitude of the first gridpoint (*1000)
**     slon:       Longitude of the first gridpoint (*1000)
**     rescomp:    Resolution and component flags
**     olon:       Longitude of grid orientation (*1000)
**     xlen:       X-direction grid length in meters
**     ylen:       Y-direction grid length in meters
**     proj:       Projection center flag
**     scan_mode:  Scanning mode flags
**
**   For Lambert Conformal Grids:
**     nx:         Number of points in the X-direction
**     ny:         Number of points in the Y-direction
**     slat:       Latitude of the first gridpoint (*1000)
**     slon:       Longitude of the first gridpoint (*1000)
**     rescomp:    Resolution and component flags
**     olon:       Longitude of grid orientation (*1000)
**     xlen:       X-direction grid length in meters
**     ylen:       Y-direction grid length in meters
**     proj:       Projection center flag
**     scan_mode:  Scanning mode flags
**     std_lat1:   First standard parallel (*1000)
**     std_lat2:   Second standard parallel (*1000)
**
**   For Mercator Grids:
**     nx:         Number of points along a latitude circle
**     ny:         Number of points along a longitude meridian
**     slat:       Latitude of the first gridpoint (*1000)
**     slon:       Longitude of the first gridpoint (*1000)
**     rescomp:    Resolution and component flags
**     elat:       Latitude of the last gridpoint (*1000)
**     elon:       Longitude of the last gridpoint (*1000)
**     std_lat1:   Standard parallel (*1000)
**     xlen:       X-direction grid length in meters
**     ylen:       Y-direction grid length in meters
**     scan_mode:  Scanning mode flags
** 
** 
**   buffer:          For internal use only (used to hold the GRIB record that
**                      was read from the GRIB data file)
**   buffer_capacity: For internal use only (the capacity of 'buffer', used to
**                      minimize memory allocations)
**   pds_ext:         This array is free-form and contains any 8-bit values that
**                      were found after the end of the standard PDS, but before
**                      the beginning of the next GRIB seciton.
**   bitmap:          For internal use only (used to hold the data point bitmap,
**                      if it exists
**   bcapacity:       For internal use only (the capacity of 'bitmap', used to
**                      minimize memory allocations)
**   bitmap_len:      For internal use only (the length of the bitmap)
**   ref_val:         GRIB Reference Value
**   gridpoints:      The array of gridpoints as a single stream - you will need
**                      to use the grid definition parameters (dimensions,
**                      scanning mode, etc.) to interpret the gridpoints
**                      properly
**   gcapacity:       For internal use only (the capacity of 'gridpoints', used
**                      to minimize memory allocations)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const double GRIB_MISSING_VALUE=1.e30;
typedef struct {
  int total_len,pds_len,pds_ext_len,gds_len,bds_len;
  int ed_num,table_ver,center_id,gen_proc,grid_type,param,level_type,lvl1,lvl2,fcst_units,p1,p2,t_range,navg,nmiss,sub_center_id,bds_flag,pack_width;
  int gds_included,bms_included;
  int yr,mo,dy,time;
  int offset;  /* offset in bytes to next GRIB section */
  int E,D;
  int data_rep,nx,ny,rescomp,scan_mode,proj;
  double slat,slon,elat,elon,lainc,loinc,olon,std_lat1,std_lat2;
  int xlen,ylen;
  unsigned char *buffer,*pds_ext,*bitmap;
  size_t buffer_capacity,bcapacity,bitmap_len;
  double ref_val,*gridpoints;
  int gcapacity;
} GRIBMessage;

/* get_bits gets the contents of the various GRIB octets
**   buf is the GRIB buffer as a stream of bytes
**   loc is the variable to hold the octet contents
**   off is the offset in BITS from the beginning of the buffer to the beginning
**       of the octet(s) to be unpacked
**   bits is the number of BITS to unpack - will be a multiple of 8 since GRIB
**       octets are 8 bits long
*/
void get_bits(unsigned char *buf,int *loc,size_t off,size_t bits)
{
/* no work to do */
  if (bits == 0) {
    return;
  }
  size_t loc_size=sizeof(int)*8;
  if (bits > loc_size) {
    fprintf(stderr,"Error: unpacking %d bits into a %d-bit field\n",bits,loc_size);
    exit(1);
  }
  else {
/* create masks to use when right-shifting (necessary because different
   compilers do different things when right-shifting a signed bit-field) */
    unsigned char bmask=1;
    size_t buf_size=sizeof(unsigned char)*8;
    for (size_t n=1; n < buf_size; ++n) {
	bmask<<=1;
	bmask++;
    }
    int lmask=1;
    for (size_t n=1; n < loc_size; ++n) {
	lmask<<=1;
	lmask++;
    }
/* get number of words to skip before unpacking begins */
    size_t wskip=off/buf_size;
/* right shift the bits in the packed buffer "word" to eliminate unneeded
   bits */
    int rshift=buf_size-(off % buf_size)-bits;
/* check for a packed field spanning multiple "words" */
    if (rshift < 0) {
	*loc=0;
	while (rshift < 0) {
	  int temp=buf[wskip++];
	  *loc+=(temp<<-rshift);
	  rshift+=buf_size;
	}
	if (rshift != 0) {
	  *loc+=(buf[wskip]>>rshift)&~(bmask<<(buf_size-rshift));
	}
	else {
	  *loc+=buf[wskip];
	}
    }
    else {
	*loc=(buf[wskip]>>rshift);
    }
/* remove any unneeded leading bits */
    if (bits != loc_size) {
	*loc&=~(lmask<<bits);
    }
  }
}

double ibm2real(unsigned char *buf,size_t off)
{
  int sign;
  get_bits(buf,&sign,off,1);
  int exp;
  get_bits(buf,&exp,off+1,7);
  exp-=64;
  int fr;
  get_bits(buf,&fr,off+8,24);
  double native_real=pow(2.,-24.)*(double)fr*pow(16.,(double)exp);
  return (sign == 1) ? -native_real : native_real;
}

void initialize(GRIBMessage *grib_msg)
{
  grib_msg->buffer=NULL;
  grib_msg->buffer_capacity=0;
  grib_msg->bitmap=NULL;
  grib_msg->bcapacity=0;
  grib_msg->bitmap_len=0;
  grib_msg->gridpoints=NULL;
  grib_msg->gcapacity=0;
}

int unpack_IS(FILE *fp,GRIBMessage *grib_msg)
{
  unsigned char temp[8];
  int status;
  if ( (status=fread(temp,1,4,fp)) != 4) {
    if (status == 0) {
	return -1;
    }
    else {
	return 1;
    }
  }
/* search for the beginning of the next GRIB message */
  if (strncmp((char *)temp,"GRIB",4) != 0) {
    while (temp[0] != 0x47 || temp[1] != 0x52 || temp[2] != 0x49 || temp[3] != 0x42) {
	switch (temp[1]) {
	  case 0x47:
	  {
	    for (size_t n=0; n < 3; ++n) {
		temp[n]=temp[n+1];
	    }
	    if ( (status=fread(&temp[3],1,1,fp)) == 0) {
		return -1;
	    }
	    break;
	  }
	  default:
	  {
	    switch(temp[2]) {
		case 0x47:
		{
		  for (size_t n=0; n < 2; ++n) {
		    temp[n]=temp[n+2];
		  }
		  if ( (status=fread(&temp[2],1,2,fp)) == 0) {
		    return -1;
		  }
		  break;
		}
		default:
		{
		  switch(temp[3]) {
		    case 0x47:
		    {
			temp[0]=temp[3];
			if ( (status=fread(&temp[1],1,3,fp)) == 0) {
			  return -1;
			}
			break;
		    }
		    default:
		    {
			if ( (status=fread(temp,1,4,fp)) == 0) {
			  return -1;
			}
		    }
		  }
		}
	    }
	  }
	}
    }
  }
  if ( (status=fread(&temp[4],1,4,fp)) == 0) {
    return 1;
  }
  get_bits(temp,&grib_msg->total_len,32,24);
  if (grib_msg->total_len == 24) {
    grib_msg->ed_num=0;
    grib_msg->pds_len=grib_msg->total_len;

/* add the four bytes for 'GRIB' + 3 bytes for the length of the section
** following the PDS */
    grib_msg->total_len+=7;
  }
  else {
    grib_msg->ed_num=1;
  }
  grib_msg->nx=grib_msg->ny=0;
  size_t required_size=grib_msg->total_len+4;
  if (required_size > grib_msg->buffer_capacity) {
    if (grib_msg->buffer != NULL) {
	free(grib_msg->buffer);
    }
    grib_msg->buffer_capacity=required_size;
    grib_msg->buffer=(unsigned char *)malloc(grib_msg->buffer_capacity*sizeof(unsigned char));
  }
  memcpy(grib_msg->buffer,temp,8);
  size_t num=grib_msg->total_len-8;
  status=fread(&grib_msg->buffer[8],1,num,fp);
  if (status != num) {
    return 1;
  }
  else {
    if (strncmp(&((char *)grib_msg->buffer)[grib_msg->total_len-4],"7777",4) != 0) {
	fprintf(stderr,"Warning: no end section found\n");
    }
    return 0;
  }
}

void unpack_PDS(GRIBMessage *grib_msg)
{
  if (grib_msg->ed_num == 0) {
    grib_msg->offset=32;
  }
  else {
    grib_msg->offset=64;
/* length of PDS */
    get_bits(grib_msg->buffer,&grib_msg->pds_len,grib_msg->offset,24);
/* table version */
    get_bits(grib_msg->buffer,&grib_msg->table_ver,grib_msg->offset+24,8);
  }
/* center ID */
  get_bits(grib_msg->buffer,&grib_msg->center_id,grib_msg->offset+32,8);
/* generating process */
  get_bits(grib_msg->buffer,&grib_msg->gen_proc,grib_msg->offset+40,8);
/* grid type */
  get_bits(grib_msg->buffer,&grib_msg->grid_type,grib_msg->offset+48,8);
  int flag;
  get_bits(grib_msg->buffer,&flag,grib_msg->offset+56,8);
/* indication of GDS */
  grib_msg->gds_included= ( (flag & 0x80) == 0x80) ? 1 : 0;
/* indication of BMS */
  grib_msg->bms_included= ( (flag & 0x40) == 0x40) ? 1 : 0;
/* parameter code */
  get_bits(grib_msg->buffer,&grib_msg->param,grib_msg->offset+64,8);
/* level type */
  get_bits(grib_msg->buffer,&grib_msg->level_type,grib_msg->offset+72,8);
  switch (grib_msg->level_type) {
    case 100:
    case 103:
    case 105:
    case 107:
    case 109:
    case 111:
    case 113:
    case 115:
    case 125:
    case 160:
    case 200:
    case 201:
    {
/* first level */
	get_bits(grib_msg->buffer,&grib_msg->lvl1,grib_msg->offset+80,16);
	grib_msg->lvl2=0;
	break;
    }
    default:
    {
/* first level */
	get_bits(grib_msg->buffer,&grib_msg->lvl1,grib_msg->offset+80,8);
/* second level */
	get_bits(grib_msg->buffer,&grib_msg->lvl2,grib_msg->offset+88,8);
    }
  }
/* year of the century */
  get_bits(grib_msg->buffer,&grib_msg->yr,grib_msg->offset+96,8);
/* month */
  get_bits(grib_msg->buffer,&grib_msg->mo,grib_msg->offset+104,8);
/* day */
  get_bits(grib_msg->buffer,&grib_msg->dy,grib_msg->offset+112,8);
/* hour */
  get_bits(grib_msg->buffer,&grib_msg->time,grib_msg->offset+120,8);
/* minutes */
  int min;
  get_bits(grib_msg->buffer,&min,grib_msg->offset+128,8);
/* complete time */
  grib_msg->time=grib_msg->time*100+min;
/* forecast time units */
  get_bits(grib_msg->buffer,&grib_msg->fcst_units,grib_msg->offset+136,8);
/* P1 */
  get_bits(grib_msg->buffer,&grib_msg->p1,grib_msg->offset+144,8);
/* P2 */
  get_bits(grib_msg->buffer,&grib_msg->p2,grib_msg->offset+152,8);
/* time range indicator*/
  get_bits(grib_msg->buffer,&grib_msg->t_range,grib_msg->offset+160,8);
  switch (grib_msg->p2) {
    case 3:
    case 4:
    case 51:
    case 113:
    case 114:
    case 115:
    case 116:
    case 117:
    case 123:
    case 124:
    {
/* number in average */
	get_bits(grib_msg->buffer,&grib_msg->navg,grib_msg->offset+168,16);
	break;
    }
    default:
    {
/* number in average */
	grib_msg->navg=0;
    }
  }
/* missing grids in average */
  get_bits(grib_msg->buffer,&grib_msg->nmiss,grib_msg->offset+184,8);
/* if GRIB0, done with PDS */
  if (grib_msg->ed_num == 0) {
    grib_msg->pds_ext_len=0;
    grib_msg->offset+=192;
    return;
  }
  int cent;
  get_bits(grib_msg->buffer,&cent,grib_msg->offset+192,8);  /* century */
  grib_msg->yr+=(cent-1)*100;
/* sub-center ID */
  get_bits(grib_msg->buffer,&grib_msg->sub_center_id,grib_msg->offset+200,8);
  int sign;
  get_bits(grib_msg->buffer,&sign,grib_msg->offset+208,1);
/* decimal scale factor */
  get_bits(grib_msg->buffer,&grib_msg->D,grib_msg->offset+209,15);
  if (sign == 1) {
    grib_msg->D=-grib_msg->D;
  }
  grib_msg->offset+=224;
  if (grib_msg->pds_len > 28) {
    if (grib_msg->pds_ext != NULL) {
	free(grib_msg->pds_ext);
	grib_msg->pds_ext=NULL;
    }
    char *c_buf=(char *)grib_msg->buffer;
    if (grib_msg->pds_len < 40) {
	fprintf(stderr,"Warning: PDS extension is in wrong location\n");
	grib_msg->pds_ext_len=grib_msg->pds_len-28;
	grib_msg->pds_ext=(unsigned char *)malloc(grib_msg->pds_ext_len);
	for (size_t n=0; n < grib_msg->pds_ext_len; ++n) {
	  grib_msg->pds_ext[n]=c_buf[36+n];
	}
	grib_msg->offset+=grib_msg->pds_ext_len*8;
    }
    else {
	grib_msg->pds_ext_len=grib_msg->pds_len-40;
	grib_msg->pds_ext=(unsigned char *)malloc(grib_msg->pds_ext_len);
	for (size_t n=0; n < grib_msg->pds_ext_len; ++n) {
	  grib_msg->pds_ext[n]=c_buf[48+n];
	}
	grib_msg->offset+=(grib_msg->pds_ext_len+12)*8;
    }
  }
  else {
    grib_msg->pds_ext_len=0;
  }
}

void unpack_GDS(GRIBMessage *grib_msg)
{
/* length of the GDS */
  get_bits(grib_msg->buffer,&grib_msg->gds_len,grib_msg->offset,24);
  if (grib_msg->ed_num == 0) {
    grib_msg->total_len+=grib_msg->gds_len;
  }
/* data representation type */
  get_bits(grib_msg->buffer,&grib_msg->data_rep,grib_msg->offset+40,8);
  switch (grib_msg->data_rep) {
/* Latitude/Longitude grid */
    case 0:
/* Gaussian Lat/Lon grid */
    case 4:
/* Rotated Lat/Lon grid */
    case 10:
    {
/* number of latitudes */
	get_bits(grib_msg->buffer,&grib_msg->nx,grib_msg->offset+48,16);
/* number of longitudes */
	get_bits(grib_msg->buffer,&grib_msg->ny,grib_msg->offset+64,16);
	int sign;
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+80,1);
	int ival;
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+81,23);
/* latitude of first gridpoint */
	grib_msg->slat=ival*0.001;
	if (sign == 1) {
	  grib_msg->slat=-grib_msg->slat;
	}
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+104,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+105,23);
/* longitude of first gridpoint */
	grib_msg->slon=ival*0.001;
	if (sign == 1) {
	  grib_msg->slon=-grib_msg->slon;
	}
/* resolution and component flags */
	get_bits(grib_msg->buffer,&grib_msg->rescomp,grib_msg->offset+128,8);
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+136,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+137,23);
/* latitude of last gridpoint */
	grib_msg->elat=ival*0.001;
	if (sign == 1) {
	  grib_msg->elat=-grib_msg->elat;
	}
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+160,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+161,23);
/* longitude of last gridpoint */
	grib_msg->elon=ival*0.001;
	if (sign == 1) {
	  grib_msg->elon=-grib_msg->elon;
	}
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+184,16);
/* longitude increment */
	grib_msg->loinc=ival*0.001;
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+200,16);
/* latitude increment */
	if (grib_msg->data_rep == 0) {
	  grib_msg->lainc=ival*0.001;
	}
	else {
	  grib_msg->lainc=ival;
	}
/* scanning mode flag */
	get_bits(grib_msg->buffer,&grib_msg->scan_mode,grib_msg->offset+216,8);
	break;
    }
/* Mercator grid */
    case 1:
    {
/* number of latitudes */
	get_bits(grib_msg->buffer,&grib_msg->nx,grib_msg->offset+48,16);
/* number of longitudes */
	get_bits(grib_msg->buffer,&grib_msg->ny,grib_msg->offset+64,16);
	int sign;
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+80,1);
	int ival;
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+81,23);
/* latitude of first gridpoint */
	grib_msg->slat=ival*0.001;
	if (sign == 1) {
	  grib_msg->slat=-grib_msg->slat;
	}
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+104,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+105,23);
/* longitude of first gridpoint */
	grib_msg->slon=ival*0.001;
	if (sign == 1) {
	  grib_msg->slon=-grib_msg->slon;
	}
/* resolution and component flags */
	get_bits(grib_msg->buffer,&grib_msg->rescomp,grib_msg->offset+128,8);
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+136,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+137,23);
/* latitude of last gridpoint */
	grib_msg->elat=ival*0.001;
	if (sign == 1) {
	  grib_msg->elat=-grib_msg->elat;
	}
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+160,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+161,23);
/* longitude of last gridpoint */
	grib_msg->elon=ival*0.001;
	if (sign == 1) {
	  grib_msg->elon=-grib_msg->elon;
	}
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+184,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+185,23);
/* standard parallel */
	grib_msg->std_lat1=ival*0.001;
	if (sign == 1) {
	  grib_msg->std_lat1=-grib_msg->std_lat1;
	}
/* scanning mode flag */
	get_bits(grib_msg->buffer,&grib_msg->scan_mode,grib_msg->offset+216,8);
/* x-direction grid length */
	get_bits(grib_msg->buffer,&grib_msg->xlen,grib_msg->offset+224,24);
/* y-direction grid length */
	get_bits(grib_msg->buffer,&grib_msg->ylen,grib_msg->offset+248,24);
	break;
    }
/* Lambert Conformal grid */
    case 3:
/* Polar Stereographic grid */
    case 5:
    {
/* number of x-points */
	get_bits(grib_msg->buffer,&grib_msg->nx,grib_msg->offset+48,16);
/* number of y-points */
	get_bits(grib_msg->buffer,&grib_msg->ny,grib_msg->offset+64,16);
	int sign;
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+80,1);
	int ival;
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+81,23);
/* latitude of first gridpoint */
	grib_msg->slat=ival*0.001;
	if (sign == 1) {
	  grib_msg->slat=-grib_msg->slat;
	}
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+104,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+105,23);
/* longitude of first gridpoint */
	grib_msg->slon=ival*0.001;
	if (sign == 1) {
	  grib_msg->slon=-grib_msg->slon;
	}
/* resolution and component flags */
	get_bits(grib_msg->buffer,&grib_msg->rescomp,grib_msg->offset+128,8);
	get_bits(grib_msg->buffer,&sign,grib_msg->offset+136,1);
	get_bits(grib_msg->buffer,&ival,grib_msg->offset+137,23);
/* longitude of grid orientation */
	grib_msg->olon=ival*0.001;
	if (sign == 1) {
	  grib_msg->olon=-grib_msg->olon;
	}
/* x-direction grid length */
	get_bits(grib_msg->buffer,&grib_msg->xlen,grib_msg->offset+160,24);
/* y-direction grid length */
	get_bits(grib_msg->buffer,&grib_msg->ylen,grib_msg->offset+184,24);
/* projection center flag */
	get_bits(grib_msg->buffer,&grib_msg->proj,grib_msg->offset+208,8);
/* scanning mode flag */
	get_bits(grib_msg->buffer,&grib_msg->scan_mode,grib_msg->offset+216,8);
	if (grib_msg->data_rep == 3) {
	  get_bits(grib_msg->buffer,&sign,grib_msg->offset+224,1);
	  get_bits(grib_msg->buffer,&ival,grib_msg->offset+225,23);
/* first standard parallel */
	  grib_msg->std_lat1=ival*0.001;
	  if (sign == 1) {
	    grib_msg->std_lat1=-grib_msg->std_lat1;
	  }
	  get_bits(grib_msg->buffer,&sign,grib_msg->offset+248,1);
	  get_bits(grib_msg->buffer,&ival,grib_msg->offset+249,23);
/* second standard parallel */
	  grib_msg->std_lat2=ival*0.001;
	  if (sign == 1) {
	    grib_msg->std_lat2=-grib_msg->std_lat2;
	  }
	}
	break;
    }
    default:
    {
	fprintf(stderr,"Grid type %d is not understood\n",grib_msg->data_rep);
	exit(1);
    }
  }
  grib_msg->offset+=grib_msg->gds_len*8;
}

void unpack_BDS(GRIBMessage *grib_msg)
{
  if (grib_msg->bms_included == 1) {
    int bms_length;
    get_bits(grib_msg->buffer,&bms_length,grib_msg->offset,24);
    if (grib_msg->ed_num == 0) {
	grib_msg->total_len+=bms_length;
    }
    int ub;
    get_bits(grib_msg->buffer,&ub,grib_msg->offset+24,8);
    int tref;
    get_bits(grib_msg->buffer,&tref,grib_msg->offset+32,16);
    if (tref != 0) {
	fprintf(stderr,"Aborting: unknown pre-defined bit-map %d\n",tref);
	exit(1);
    }
    grib_msg->bitmap_len=(bms_length-6)*8-ub;
    if (grib_msg->bitmap_len > grib_msg->bcapacity) {
	if (grib_msg->bitmap != NULL) {
	  free(grib_msg->bitmap);
	}
	grib_msg->bcapacity=grib_msg->bitmap_len;
	grib_msg->bitmap=(unsigned char *)malloc(grib_msg->bcapacity*sizeof(unsigned char));
    }
    size_t boff=grib_msg->offset+48;
    for (size_t n=0; n < grib_msg->bcapacity; ++n) {
	int bval;
	get_bits(grib_msg->buffer,&bval,boff,1);
	grib_msg->bitmap[n]=bval;
	++boff;
    }
    grib_msg->offset+=bms_length*8;
  }
  else {
    grib_msg->bitmap_len=0;
  }
/* length of the BDS */
  get_bits(grib_msg->buffer,&grib_msg->bds_len,grib_msg->offset,24);
  if (grib_msg->ed_num == 0) {
    grib_msg->total_len+=(grib_msg->bds_len+1);
  }
/* flag */
  get_bits(grib_msg->buffer,&grib_msg->bds_flag,grib_msg->offset+24,4);
  int ub;
  get_bits(grib_msg->buffer,&ub,grib_msg->offset+28,4);
/* bit width of the packed data points */
  get_bits(grib_msg->buffer,&grib_msg->pack_width,grib_msg->offset+80,8);
  int sign;
  get_bits(grib_msg->buffer,&sign,grib_msg->offset+32,1);
/* binary scale factor */
  get_bits(grib_msg->buffer,&grib_msg->E,grib_msg->offset+33,15);
  if (sign == 1) {
    grib_msg->E=-grib_msg->E;
  }
  double e=pow(2.,grib_msg->E);
/* reference value */
  double d=pow(10.,grib_msg->D);
  grib_msg->ref_val=ibm2real(grib_msg->buffer,grib_msg->offset+48)/d;
  if ((grib_msg->bds_flag & 0x40) == 0) {
/* simple packing */
    int *packed=NULL;
    grib_msg->offset+=88;
    size_t num_packed=0;
    if (grib_msg->pack_width > 0) {
	num_packed=(grib_msg->bds_len*8-88-ub)/grib_msg->pack_width;
	packed=(int *)malloc(sizeof(int)*num_packed);
    }
    switch (grib_msg->data_rep) {
	case 0:
/* Latitude/Longitude grid */
	case 4:
/* Gaussian Lat/Lon grid */
	case 10:
/* Rotated Lat/Lon grid */
	{
	  switch (grib_msg->grid_type) {
	    case 23:
	    case 24:
	    case 26:
	    case 63:
	    case 64:
	    {
		grib_msg->offset+=grib_msg->pack_width;
		break;
	    }
	  }
	}
	case 1:
/* Mercator grid */
	case 3:
/* Lambert Conformal grid */
	case 5:
/* Polar Stereographic grid */
	{
	  for (size_t n=0; n < num_packed; ++n) {
	    get_bits(grib_msg->buffer,&packed[n],grib_msg->offset,grib_msg->pack_width);
	    grib_msg->offset+=grib_msg->pack_width;
	  }
	  size_t num_points=grib_msg->ny*grib_msg->nx;
	  if (num_points > grib_msg->gcapacity) {
	    if (grib_msg->gridpoints != NULL) {
		free(grib_msg->gridpoints);
	    }
	    grib_msg->gcapacity=num_points;
	    grib_msg->gridpoints=(double *)malloc(grib_msg->gcapacity*sizeof(double));
	  }
	  size_t bcnt=0;
	  size_t pcnt=0;
	  for (size_t n=0; n < num_points; ++n) {
	    if (grib_msg->bitmap_len == 0 || (grib_msg->bitmap_len > 0 && grib_msg->bitmap[bcnt++] == 1)) {
		if (packed == NULL) {
/* constant field */
		  grib_msg->gridpoints[n]=grib_msg->ref_val;
		}
		else {
		  grib_msg->gridpoints[n]=grib_msg->ref_val+packed[pcnt++]*e/d;
		}
	    }
	    else {
		grib_msg->gridpoints[n]=GRIB_MISSING_VALUE;
	    }
	  }
	  break;
	}
	default:
/* no recognized GDS, so just unpack the stream of gridpoints */
	{
	  size_t num_points= (num_packed > grib_msg->bitmap_len) ? num_packed : grib_msg->bitmap_len;
	  if (num_points > grib_msg->gcapacity) {
	    if (grib_msg->gridpoints != NULL) {
		free(grib_msg->gridpoints);
	    }
	    grib_msg->gcapacity=num_points;
	    grib_msg->gridpoints=(double *)malloc(grib_msg->gcapacity*sizeof(double));
	  }
	  size_t bcnt=0;
	  size_t pcnt=0;
	  for (size_t n=0; n < num_points; ++n) {
	    if (grib_msg->bitmap_len == 0 || (grib_msg->bitmap_len > 0 && grib_msg->bitmap[bcnt++] == 1)) {
		grib_msg->gridpoints[n]=grib_msg->ref_val+packed[pcnt++]*e/d;
	    }
	    else {
		grib_msg->gridpoints[n]=GRIB_MISSING_VALUE;
	    }
	  }
	}
    }
    if (num_packed > 0) {
	free(packed);
    }
  }
  else {
/* second-order packing */
    fprintf(stderr,"Aborting: complex packing not currently supported\n");
    exit(1);
  }
}

int unpackgrib1(FILE *fp,GRIBMessage *grib_msg)
{
  int status;
  if ( (status=unpack_IS(fp,grib_msg)) != 0) {
    return status;
  }
  unpack_PDS(grib_msg);
  if (grib_msg->gds_included == 1) {
    unpack_GDS(grib_msg);
  }
  unpack_BDS(grib_msg);
  return 0;
}
