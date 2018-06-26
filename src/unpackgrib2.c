/*
** File: unpackgrib2.c
**
** Author:  Bob Dattore
**          NCAR/DSS
**          dattore@ucar.edu
**          (303) 497-1825
**
** Revision History:
**          12 Feb 2008 - first version, unpacks:
**             GDS Template 3.0 (Latitude/Longitude grid)
**             GDS Template 3.40 (Gaussian Latitude/Longitude grid)
**             PDS Template 4.0 (Analysis/Forecast on level/layer at a point in
**                               time)
**             DRS Template 5.0 (Simple grid point packing)
**             DRS Template 5.40 (JPEG 2000 packing)
**              **NOTE: JPEG 2000 unpacking requires libjasper.a and libjpeg.a
**                      You can get the JasPer code at
**                        http://www.ece.uvic.ca/~mdadams/jasper/
**                      You can get the installation code for libjpeg.a at
**                        http://www.ijg.org/files/jpegsrc.v6b.tar.gz
**          31 Mar 2008:
**             GDS Template 3.30 (Lambert conformal grid)
**             PDS Template 4.8 (Average, accumulation, etc.)
**          03 Apr 2008:
**             Some bug fixes
**          16 Jul 2008:
**             PDS Template 4.2 (Derived forecast based on all ensemble members)
**          24 Oct 2008:
**             Some bug fixes
**             PDS Template 4.1 (Individual ensemble forecast at a point in
**                               time)
**             PDS Template 4.11 (Individual ensemble forecast of average,
**                                accumulation, etc.)
**          12 Nov 2008:
**             DRS Template 5.40000 (same as 5.40, the local value for JPEG 2000
**                                   packing before it was approved by WMO; some
**                                   centers still use this value)
**          29 Nov 2013:
**             PDS Template 4.15 (U.K. Met Office)
**          14 Aug 2015:
**             DRS Template 5.3 (complex packing and spatial differencing)
**          17 Apr 2017:
**             removed restriction that grid must be oriented north-to-south to
**               unpack DRS template 5.3
**          13 Jul 2017:
**             GDS Template 3.10 (Mercator grid)
**
** Purpose: to provide a single C-routine for unpacking GRIB2 messages
**
** Notes:   1) There are several routines defined in this file, any of which can
**             be called independently.  However, to unpack a GRIB2 message with
**             one function call, use only "unpackgrib2".
**
**          2) The user is expected to have some understanding of the GRIB2
**             format, as some of the information consists of codes, etc. that
**             are defined by the GRIB2 format.  You can find a description of
**             the GRIB2 format at
**             https://rda.ucar.edu/docs/formats/grib2/grib2doc/.
**
**          3) please report any problems to dattore@ucar.edu.
**
** example C syntax for using unpackgrib2:
**    FILE *fp;
**    GRIB2Message grib2_msg;
**    int status;
**
**    initialize(&grib2_msg);
**    fp=fopen("my_GRIB2_file","rb");
**    while ( (status=unpackgrib2(fp,&grib2_msg)) == 0) {
**      if (status == -1) {
**        printf("Found EOF\n");
**      }
**      else if (status == 1) {
**        printf("Error reading GRIB2 message\n");
**      }
**      ...
**    }
** 
** where:
**   fp            is a FILE pointer to an open GRIB data file
**   GRIB2Message  is a structure used to hold the data in a GRIB2 message
**
** On return:
**   unpackgrib2 returns 0 for a successful read, -1 for an EOF, and 1 for a
**   read error
**
** Overview of the GRIB2Message structure:
**   buffer:          For internal use only (used to hold the GRIB2 message that
**                      was read from the GRIB2 data file)
**   buffer_capacity: For internal use only (the capacity of 'buffer', used to
**                      minimize memory allocations
**   offset:          For internal use only (offset in bytes to next GRIB2
**                      section from the beginning of the message)
**   total_len:       Total length of the GRIB2 message, in octets (8-bit bytes)
**   disc:            Discipline number
**   ed_num:          Edition number (2)
**   center_id:       Center ID
**   sub_center_id:   Sub-center ID
**   table_ver:       GRIB2 parameter table version
**   local_table_ver: Local parameter table version
**   ref_time_type:   Signficance of reference time
**   yr:              Year of the reference time (4-digits - YYYY)
**   mo:              Month of the reference time
**   dy:              Day of the reference time
**   time:            Time of the reference time (HHMMSS; hours/minutes/seconds)
**   prod_status:     Production status of data
**   data_type:       Type of data
**   md:        Metadata that is common to all grids in the message
**   num_grids:       Number of individual grids in the GRIB2 message
**   grids:           Array of individual grids
**   grid_capacity:   For internal use only (the capacity of 'grids', used to
**                      minimize memory allocations
**
** Overview of the GRIB2Metadata structure:
**   gds_templ_num:   Grid definition template number
**   earth_shape:     Code for shape of the earth
**   For Latitude/longitude grids:
**     nx:              Number of points along a latitude circle
**     ny:              Number of points along a longitude meridian
**     slat:            Latitude of the first gridpoint
**     slon:            Longitude of the first gridpoint
**     rescomp:         Resolution and component flags
**     elat:            Latitude of the last gridpoint
**     elon:            Longitude of the last gridpoint
**     loinc:           Longitude increment
**     lainc:           Latitude increment
**     scan_mode:       Scanning mode flags
**
**   For Mercator grids:
**     nx:              Number of points along a latitude circle
**     ny:              Number of points along a longitude meridian
**     slat:            Latitude of the first gridpoint
**     slon:            Longitude of the first gridpoint
**     rescomp:         Resolution and component flags
**     elat:            Latitude of the last gridpoint
**     elon:            Longitude of the last gridpoint
**     latin1:          Standard parallel
**     dxinc:           X-direction grid length in meters
**     dyinc:           Y-direction grid length in meters
**     scan_mode:       Scanning mode flags
**
**   For Lambert conformal grids:
**     nx:              Number of points along a latitude circle
**     ny:              Number of points along a longitude meridian
**     slat:            Latitude of the first gridpoint
**     slon:            Longitude of the first gridpoint
**     rescomp:         Resolution and component flags
**     lad:             Latitude at which x- and y-direction grid lenths are
**                      valid
**     lov:             Longitude of meridian parallel to y-axis along which
**                      latitude increases as y increases
**     dxinc:           X-direction grid length
**     dyinc:           Y-direction grid length
**     proj_flag:       Projection center flag
**     scan_mode:       Scanning mode flags
**     latin1:          First latitude from pole at which secant cone cuts the
**                      sphere
**     latin2:          Second latitude from pole at which secant cone cuts the
**                      sphere
**     splat:           Latitude of southern pole of projection
**     splon:           Longitude of southern pole of projection
**
**   pds_templ_num:   Production definition template number
**   param_cat:       GRIB2 parameter category
**   param_num:       GRIB2 parameter number
**   gen_proc:        Type of generating process
**   time_unit:       Indicator of unit of time range
**   fcst_time:       Forecast time
**   lvl1_type:       Type of first level
**   lvl1:            Value of first level
**   lvl2_type:       Type of second level
**   lvl2:            Value of second level
**   drs_templ_num:   Data representatio template number
**   R:               Reference value
**   E:               Binary scale factor
**   D:               Decimal scale factor
**   num_packed:      Number of packed values in the Data Section
**   pack_width:      Number of bits used for each packed data value
**   bms_ind:         Bit map indicator
**   bitmap:          Buffer to hold the bitmap
**
** Overview of the GRIB2Grid structure:
**   md:          Metadata that is common to all grids in the message
**   gridpoints:  The array of gridpoints as a single stream - you will need to
**                  use the grid definition parameters (dimensions, scanning
**                  mode, etc.) to interpret the gridpoints properly
**   gcapacity:   For internal use only (the capacity of 'gridpoints', used to
**                  minimize memory allocations)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef JASPER
#include <jasper/jasper.h>
#endif

#ifdef JASPER
int dec_jpeg2000(char *injpc,int bufsize,int *outfld)
/*$$$  SUBPROGRAM DOCUMENTATION BLOCK
*                .      .    .                                       .
* SUBPROGRAM:    dec_jpeg2000      Decodes JPEG2000 code stream
*   PRGMMR: Gilbert          ORG: W/NP11     DATE: 2002-12-02
*
* ABSTRACT: This Function decodes a JPEG2000 code stream specified in the
*   JPEG2000 Part-1 standard (i.e., ISO/IEC 15444-1) using JasPer 
*   Software version 1.500.4 (or 1.700.2) written by the University of British
*   Columbia and Image Power Inc, and others.
*   JasPer is available at http://www.ece.uvic.ca/~mdadams/jasper/.
*
* PROGRAM HISTORY LOG:
* 2002-12-02  Gilbert
*
* USAGE:     int dec_jpeg2000(char *injpc,int bufsize,int *outfld)
*
*   INPUT ARGUMENTS:
*      injpc - Input JPEG2000 code stream.
*    bufsize - Length (in bytes) of the input JPEG2000 code stream.
*
*   OUTPUT ARGUMENTS:
*     outfld - Output matrix of grayscale image values.
*
*   RETURN VALUES :
*          0 = Successful decode
*         -3 = Error decode jpeg2000 code stream.
*         -5 = decoded image had multiple color components.
*              Only grayscale is expected.
*
* REMARKS:
*
*      Requires JasPer Software version 1.500.4 or 1.700.2
*
* ATTRIBUTES:
*   LANGUAGE: C
*   MACHINE:  IBM SP
*
*$$$*/

{
    int ier;
    int i,j,k;
    jas_image_t *image=0;
    jas_stream_t *jpcstream;
    jas_image_cmpt_t *pcmpt;
    char *opts=0;
    jas_matrix_t *data;

//    jas_init();

    ier=0;
//   
//     Create jas_stream_t containing input JPEG200 codestream in memory.
//       

    jpcstream=jas_stream_memopen(injpc,bufsize);

//   
//     Decode JPEG200 codestream into jas_image_t structure.
//       
    image=jpc_decode(jpcstream,opts);
    if ( image == 0 ) {
       printf(" jpc_decode return = %d \n",ier);
       return -3;
    }
    
    pcmpt=image->cmpts_[0];

//   Expecting jpeg2000 image to be grayscale only.
//   No color components.
//
    if (image->numcmpts_ != 1 ) {
       printf("dec_jpeg2000: Found color image.  Grayscale expected.\n");
       return (-5);
    }

// 
//    Create a data matrix of grayscale image values decoded from
//    the jpeg2000 codestream.
//
    data=jas_matrix_create(jas_image_height(image), jas_image_width(image));
    jas_image_readcmpt(image,0,0,0,jas_image_width(image),
                       jas_image_height(image),data);
//
//    Copy data matrix to output integer array.
//
    k=0;
    for (i=0;i<pcmpt->height_;i++) 
      for (j=0;j<pcmpt->width_;j++) 
        outfld[k++]=data->rows_[i][j];
//
//     Clean up JasPer work structures.
//
    jas_matrix_destroy(data);
    ier=jas_stream_close(jpcstream);
    jas_image_destroy(image);

    return 0;

}
#endif

const double GRIB_MISSING_VALUE=1.e30;

typedef struct {
  int gds_templ_num;
  int earth_shape;
  int nx,ny;
  double slat,slon,latin1,latin2,splat,splon;
  union {
    double elat;
    double lad;
  } lats;
  union {
    double elon;
    double lov;
  } lons;
  union {
    double loinc;
    double dxinc;
  } xinc;
  union {
    double lainc;
    double dyinc;
  } yinc;
  int rescomp,scan_mode,proj_flag;
  int pds_templ_num;
  int param_cat,param_num,gen_proc,time_unit,fcst_time;
  int ens_type,perturb_num,derived_fcst_code,nfcst_in_ensemble;
  int lvl1_type,lvl2_type;
  double lvl1,lvl2;
  struct {
    int eyr,emo,edy,etime;
    int num_ranges,nmiss;
    int *proc_code,*incr_type,*time_unit,*time_length,*incr_unit,*incr_length;
  } stat_proc;
  struct {
    int stat_proc,type,num_points;
  } spatial_proc;
  struct {
    int split_method,miss_val_mgmt;
    int num_groups;
    float primary_miss_sub,secondary_miss_sub;
    struct {
	int ref,pack_width;
    } width;
    struct {
	int ref,incr,last,pack_width;
    } length;
    struct {
	int order,order_vals_width;
    } spatial_diff;
  } complex_pack;
  int drs_templ_num;
  float R;
  int E,D,num_packed,pack_width,orig_val_type;
  int bms_ind;
  unsigned char *bitmap;
} GRIB2Metadata;

typedef struct {
  GRIB2Metadata md;
  double *gridpoints;
  size_t gcapacity;
} GRIB2Grid;

typedef struct {
  unsigned char *buffer;
  size_t buffer_capacity;
  int offset;  /* offset in bytes to next GRIB2 section */
  int total_len,disc,ed_num;
  int center_id,sub_center_id,table_ver,local_table_ver,ref_time_type;
  int yr,mo,dy,time;
  int prod_status,data_type;
  GRIB2Metadata md;
  int num_grids;
  GRIB2Grid *grids;
  size_t grid_capacity;
} GRIB2Message;

/* get_bits gets the contents of the various GRIB octets
**   buf is the GRIB2 buffer as a stream of bytes
**   loc is the variable to hold the octet contents
**   off is the offset in BITS from the beginning of the buffer to the beginning
**       of the octet(s) to be unpacked
**   bits is the number of BITS to unpack - will be a multiple of 8 since GRIB2
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

void initialize(GRIB2Message *grib2_msg)
{
  grib2_msg->buffer=NULL;
  grib2_msg->buffer_capacity=0;
  grib2_msg->grids=NULL;
  grib2_msg->grid_capacity=0;
  grib2_msg->md.stat_proc.proc_code=NULL;
}

int unpack_IS(FILE *fp,GRIB2Message *grib2_msg)
{
  grib2_msg->num_grids=0;
  unsigned char temp[16];
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
  if ( (status=fread(&temp[4],1,12,fp)) == 0) {
    return 1;
  }
  get_bits(temp,&grib2_msg->disc,48,8);
  get_bits(temp,&grib2_msg->ed_num,56,8);
  get_bits(temp,&grib2_msg->total_len,96,32);
  grib2_msg->md.nx=grib2_msg->md.ny=0;
  size_t required_size=grib2_msg->total_len+4;
  if (required_size > grib2_msg->buffer_capacity) {
    if (grib2_msg->buffer != NULL) {
	free(grib2_msg->buffer);
    }
    grib2_msg->buffer_capacity=required_size;
    grib2_msg->buffer=(unsigned char *)malloc(grib2_msg->buffer_capacity*sizeof(unsigned char));
  }
  memcpy(grib2_msg->buffer,temp,16);
  size_t num=grib2_msg->total_len-16;
  status=fread(&grib2_msg->buffer[16],1,num,fp);
  if (status != num) {
    return 1;
  }
  else {
    if (strncmp(&((char *)grib2_msg->buffer)[grib2_msg->total_len-4],"7777",4) != 0) {
	fprintf(stderr,"Warning: no end section found\n");
    }
    grib2_msg->offset=128;
    return 0;
  }
}

void unpack_IDS(GRIB2Message *grib2_msg)
{
  int length;
/* length of the IDS */
  get_bits(grib2_msg->buffer,&length,grib2_msg->offset,32);
/* center ID */
  get_bits(grib2_msg->buffer,&grib2_msg->center_id,grib2_msg->offset+40,16);
/* sub-center ID */
  get_bits(grib2_msg->buffer,&grib2_msg->sub_center_id,grib2_msg->offset+56,16);
/* table version */
  get_bits(grib2_msg->buffer,&grib2_msg->table_ver,grib2_msg->offset+72,8);
/* local table version */
  get_bits(grib2_msg->buffer,&grib2_msg->local_table_ver,grib2_msg->offset+80,8);
/* significance of reference time */
  get_bits(grib2_msg->buffer,&grib2_msg->ref_time_type,grib2_msg->offset+88,8);
/* year */
  get_bits(grib2_msg->buffer,&grib2_msg->yr,grib2_msg->offset+96,16);
/* month */
  get_bits(grib2_msg->buffer,&grib2_msg->mo,grib2_msg->offset+112,8);
/* day */
  get_bits(grib2_msg->buffer,&grib2_msg->dy,grib2_msg->offset+120,8);
/* hours */
  int hh;
  get_bits(grib2_msg->buffer,&hh,grib2_msg->offset+128,8);
/* minutes */
  int mm;
  get_bits(grib2_msg->buffer,&mm,grib2_msg->offset+136,8);
/* seconds */
  int ss;
  get_bits(grib2_msg->buffer,&ss,grib2_msg->offset+144,8);
  grib2_msg->time=hh*10000+mm*100+ss;
/* production status */
  get_bits(grib2_msg->buffer,&grib2_msg->prod_status,grib2_msg->offset+152,8);
/* type of data */
  get_bits(grib2_msg->buffer,&grib2_msg->data_type,grib2_msg->offset+160,8);
  grib2_msg->offset+=length*8;
}

void unpack_LUS(GRIB2Message *grib2_msg)
{
}

void unpack_GDS(GRIB2Message *grib2_msg)
{
/* source of grid definition */
  int src;
  get_bits(grib2_msg->buffer,&src,grib2_msg->offset+40,8);
  if (src != 0) {
    fprintf(stderr,"Don't recognize predetermined grid definitions");
    exit(1);
  }
/* quasi-regular grid indication */
  int num_in_list;
  get_bits(grib2_msg->buffer,&num_in_list,grib2_msg->offset+80,8);
  if (num_in_list > 0) {
    fprintf(stderr,"Unable to unpack quasi-regular grids");
    exit(1);
  }
/* grid definition template number */
  get_bits(grib2_msg->buffer,&grib2_msg->md.gds_templ_num,grib2_msg->offset+96,16);
  switch (grib2_msg->md.gds_templ_num) {
/* Latitude/longitude grid */
    case 0:
    case 40:
    {
/* shape of the earth */
	get_bits(grib2_msg->buffer,&grib2_msg->md.earth_shape,grib2_msg->offset+112,8);
/* number of latitudes */
	get_bits(grib2_msg->buffer,&grib2_msg->md.nx,grib2_msg->offset+240,32);
/* number of longitudes */
	get_bits(grib2_msg->buffer,&grib2_msg->md.ny,grib2_msg->offset+272,32);
/* latitude of first gridpoint */
	int sign;
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+368,1);
	int value;
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+369,31);
	grib2_msg->md.slat=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.slat=-grib2_msg->md.slat;
	}
/* longitude of first gridpoint */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+400,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+401,31);
	grib2_msg->md.slon=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.slon=-grib2_msg->md.slon;
	}
/* resolution and component flags */
	get_bits(grib2_msg->buffer,&grib2_msg->md.rescomp,grib2_msg->offset+432,8);
/* latitude of last gridpoint */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+440,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+441,31);
	grib2_msg->md.lats.elat=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.lats.elat=-grib2_msg->md.lats.elat;
	}
/* longitude of last gridpoint */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+472,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+473,31);
	grib2_msg->md.lons.elon=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.lons.elon=-grib2_msg->md.lons.elon;
	}
/* longitude increment */
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+504,32);
	grib2_msg->md.xinc.loinc=value/1000000.;
/* latitude increment */
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+536,32);
	if (grib2_msg->md.gds_templ_num == 0) {
	  grib2_msg->md.yinc.lainc=value/1000000.;
	}
/* scanning mode flag */
	get_bits(grib2_msg->buffer,&grib2_msg->md.scan_mode,grib2_msg->offset+568,8);
	break;
    }
/* Mercator grid */
    case 10:
    {
	get_bits(grib2_msg->buffer,&grib2_msg->md.earth_shape,grib2_msg->offset+112,8);
/* number of points along a parallel */
	get_bits(grib2_msg->buffer,&grib2_msg->md.nx,grib2_msg->offset+240,32);
/* number of points along a meridian */
	get_bits(grib2_msg->buffer,&grib2_msg->md.ny,grib2_msg->offset+272,32);
/* latitude of first gridpoint */
	int sign;
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+304,1);
	int value;
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+305,31);
	grib2_msg->md.slat=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.slat=-grib2_msg->md.slat;
	}
/* longitude of first gridpoint */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+336,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+337,31);
	grib2_msg->md.slon=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.slon=-grib2_msg->md.slon;
	}
/* resolution and component flags */
	get_bits(grib2_msg->buffer,&grib2_msg->md.rescomp,grib2_msg->offset+368,8);
/* latin1 */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+376,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+377,31);
	grib2_msg->md.latin1=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.latin1=-grib2_msg->md.latin1;
	}
/* latitude of last gridpoint */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+408,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+409,31);
	grib2_msg->md.lats.elat=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.lats.elat=-grib2_msg->md.lats.elat;
	}
/* longitude of last gridpoint */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+440,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+441,31);
	grib2_msg->md.lons.elon=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.lons.elon=-grib2_msg->md.lons.elon;
	}
/* scanning mode flag */
	get_bits(grib2_msg->buffer,&grib2_msg->md.scan_mode,grib2_msg->offset+472,8);
/* x-direction increment */
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+512,32);
	grib2_msg->md.xinc.dxinc=value/1000.;
/* y-direction increment */
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+544,32);
	grib2_msg->md.yinc.dyinc=value/1000.;
	break;
    }
/* Lambert conformal grid */
    case 30:
    {
	get_bits(grib2_msg->buffer,&grib2_msg->md.earth_shape,grib2_msg->offset+112,8);
/* number of points along a parallel */
	get_bits(grib2_msg->buffer,&grib2_msg->md.nx,grib2_msg->offset+240,32);
/* number of points along a meridian */
	get_bits(grib2_msg->buffer,&grib2_msg->md.ny,grib2_msg->offset+272,32);
/* latitude of first gridpoint */
	int sign;
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+304,1);
	int value;
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+305,31);
	grib2_msg->md.slat=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.slat=-grib2_msg->md.slat;
	}
/* longitude of first gridpoint */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+336,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+337,31);
	grib2_msg->md.slon=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.slon=-grib2_msg->md.slon;
	}
/* resolution and component flags */
	get_bits(grib2_msg->buffer,&grib2_msg->md.rescomp,grib2_msg->offset+368,8);
/* LaD */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+376,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+377,31);
	grib2_msg->md.lats.lad=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.lats.lad=-grib2_msg->md.lats.lad;
	}
/* LoV */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+408,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+409,31);
	grib2_msg->md.lons.lov=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.lons.lov=-grib2_msg->md.lons.lov;
	}
/* x-direction increment */
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+440,32);
	grib2_msg->md.xinc.dxinc=value/1000.;
/* y-direction increment */
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+472,32);
	grib2_msg->md.yinc.dyinc=value/1000.;
/* projection center flag */
	get_bits(grib2_msg->buffer,&grib2_msg->md.proj_flag,grib2_msg->offset+504,8);
/* scanning mode flag */
	get_bits(grib2_msg->buffer,&grib2_msg->md.scan_mode,grib2_msg->offset+512,8);
/* latin1 */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+520,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+521,31);
	grib2_msg->md.latin1=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.latin1=-grib2_msg->md.latin1;
	}
/* latin2 */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+552,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+553,31);
	grib2_msg->md.latin2=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.latin2=-grib2_msg->md.latin2;
	}
/* latitude of southern pole of projection */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+584,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+585,31);
	grib2_msg->md.splat=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.splat=-grib2_msg->md.splat;
	}
/* longitude of southern pole of projection */
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+616,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+617,31);
	grib2_msg->md.splon=value/1000000.;
	if (sign == 1) {
	  grib2_msg->md.splon=-grib2_msg->md.splon;
	}
	break;
    }
    default:
    {
	fprintf(stderr,"Grid template %d is not understood\n",grib2_msg->md.gds_templ_num);
	exit(1);
    }
  }
}

void unpack_PDS(GRIB2Message *grib2_msg)
{
/* indication of hybrid coordinate system */
  int num_coords;
  get_bits(grib2_msg->buffer,&num_coords,grib2_msg->offset+40,16);
  if (num_coords > 0) {
    fprintf(stderr,"Unable to decode hybrid coordinates");
    exit(1);
  }
/* product definition template number */
  get_bits(grib2_msg->buffer,&grib2_msg->md.pds_templ_num,grib2_msg->offset+56,16);
  grib2_msg->md.stat_proc.num_ranges=0;
  switch (grib2_msg->md.pds_templ_num) {
    case 0:
    case 1:
    case 2:
    case 8:
    case 11:
    case 12:
    case 15:
    {
	grib2_msg->md.ens_type=-1;
	grib2_msg->md.derived_fcst_code=-1;
	grib2_msg->md.spatial_proc.type=-1;
/* parameter category */
	get_bits(grib2_msg->buffer,&grib2_msg->md.param_cat,grib2_msg->offset+72,8);
/* parameter number */
	get_bits(grib2_msg->buffer,&grib2_msg->md.param_num,grib2_msg->offset+80,8);
/* generating process */
	get_bits(grib2_msg->buffer,&grib2_msg->md.gen_proc,grib2_msg->offset+88,8);
/* time range indicator*/
	get_bits(grib2_msg->buffer,&grib2_msg->md.time_unit,grib2_msg->offset+136,8);
/* forecast time */
	get_bits(grib2_msg->buffer,&grib2_msg->md.fcst_time,grib2_msg->offset+144,32);
/* type of first level */
	get_bits(grib2_msg->buffer,&grib2_msg->md.lvl1_type,grib2_msg->offset+176,8);
/* value of first level */
	int factor;
	get_bits(grib2_msg->buffer,&factor,grib2_msg->offset+184,8);
	int sign;
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+192,1);
	int value;
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+193,31);
	if (sign == 1) {
	  value=-value;
	}
	grib2_msg->md.lvl1=(double)value/pow(10.,(double)factor);
/* type of second level */
	get_bits(grib2_msg->buffer,&grib2_msg->md.lvl2_type,grib2_msg->offset+224,8);
/* value of second level */
	get_bits(grib2_msg->buffer,&factor,grib2_msg->offset+232,8);
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+240,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+241,31);
	if (sign == 1) {
	  value=-value;
	}
	grib2_msg->md.lvl2=(double)value/pow(10.,(double)factor);
	switch (grib2_msg->md.pds_templ_num) {
	  case 1:
	  case 11:
	  {
	    get_bits(grib2_msg->buffer,&grib2_msg->md.ens_type,grib2_msg->offset+272,8);
	    get_bits(grib2_msg->buffer,&grib2_msg->md.perturb_num,grib2_msg->offset+280,8);
	    get_bits(grib2_msg->buffer,&grib2_msg->md.nfcst_in_ensemble,grib2_msg->offset+288,8);
	    switch (grib2_msg->md.pds_templ_num) {
		case 11:
		{
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.eyr,grib2_msg->offset+296,16);
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.emo,grib2_msg->offset+312,8);
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.edy,grib2_msg->offset+320,8);
		  int hh;
		  get_bits(grib2_msg->buffer,&hh,grib2_msg->offset+328,8);
		  int mm;
		  get_bits(grib2_msg->buffer,&mm,grib2_msg->offset+336,8);
		  int ss;
		  get_bits(grib2_msg->buffer,&ss,grib2_msg->offset+344,8);
		  grib2_msg->md.stat_proc.etime=hh*10000+mm*100+ss;
/* number of time range specifications */
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.num_ranges,grib2_msg->offset+352,8);
/* number of values missing from process */
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.nmiss,grib2_msg->offset+360,32);
		  if (grib2_msg->md.stat_proc.proc_code != NULL) {
		    free(grib2_msg->md.stat_proc.proc_code);
		    free(grib2_msg->md.stat_proc.incr_type);
		    free(grib2_msg->md.stat_proc.time_unit);
		    free(grib2_msg->md.stat_proc.time_length);
		    free(grib2_msg->md.stat_proc.incr_unit);
		    free(grib2_msg->md.stat_proc.incr_length);
		    grib2_msg->md.stat_proc.proc_code=NULL;
		  }
		  grib2_msg->md.stat_proc.proc_code=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.incr_type=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.time_unit=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.time_length=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.incr_unit=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.incr_length=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  size_t off=392;
		  for (size_t n=0; n < grib2_msg->md.stat_proc.num_ranges; ++n) {
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.proc_code[n],grib2_msg->offset+off,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_type[n],grib2_msg->offset+off+8,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.time_unit[n],grib2_msg->offset+off+16,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.time_length[n],grib2_msg->offset+off+24,32);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_unit[n],grib2_msg->offset+off+56,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_length[n],grib2_msg->offset+off+64,32);
		    off+=96;
		  }
		  break;
		}
	    }
	    break;
	  }
	  case 2:
	  case 12:
	  {
	    get_bits(grib2_msg->buffer,&grib2_msg->md.derived_fcst_code,grib2_msg->offset+272,8);
	    get_bits(grib2_msg->buffer,&grib2_msg->md.nfcst_in_ensemble,grib2_msg->offset+280,8);
	    switch (grib2_msg->md.pds_templ_num) {
		case 12:
		{
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.eyr,grib2_msg->offset+288,16);
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.emo,grib2_msg->offset+304,8);
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.edy,grib2_msg->offset+312,8);
		  int hh;
		  get_bits(grib2_msg->buffer,&hh,grib2_msg->offset+320,8);
		  int mm;
		  get_bits(grib2_msg->buffer,&mm,grib2_msg->offset+328,8);
		  int ss;
		  get_bits(grib2_msg->buffer,&ss,grib2_msg->offset+336,8);
		  grib2_msg->md.stat_proc.etime=hh*10000+mm*100+ss;
/* number of time range specifications */
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.num_ranges,grib2_msg->offset+344,8);
/* number of values missing from process */
		  get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.nmiss,grib2_msg->offset+352,32);
		  if (grib2_msg->md.stat_proc.proc_code != NULL) {
		    free(grib2_msg->md.stat_proc.proc_code);
		    free(grib2_msg->md.stat_proc.incr_type);
		    free(grib2_msg->md.stat_proc.time_unit);
		    free(grib2_msg->md.stat_proc.time_length);
		    free(grib2_msg->md.stat_proc.incr_unit);
		    free(grib2_msg->md.stat_proc.incr_length);
		    grib2_msg->md.stat_proc.proc_code=NULL;
		  }
		  grib2_msg->md.stat_proc.proc_code=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.incr_type=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.time_unit=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.time_length=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.incr_unit=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  grib2_msg->md.stat_proc.incr_length=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
		  size_t off=384;
		  for (size_t n=0; n < grib2_msg->md.stat_proc.num_ranges; ++n) {
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.proc_code[n],grib2_msg->offset+off,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_type[n],grib2_msg->offset+off+8,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.time_unit[n],grib2_msg->offset+off+16,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.time_length[n],grib2_msg->offset+off+24,32);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_unit[n],grib2_msg->offset+off+56,8);
		    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_length[n],grib2_msg->offset+off+64,32);
		    off+=96;
		  }
		  break;
		}
	    }
	    break;
	  }
	  case 8:
	  {
	    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.eyr,grib2_msg->offset+272,16);
	    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.emo,grib2_msg->offset+288,8);
	    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.edy,grib2_msg->offset+296,8);
	    int hh;
	    get_bits(grib2_msg->buffer,&hh,grib2_msg->offset+304,8);
	    int mm;
	    get_bits(grib2_msg->buffer,&mm,grib2_msg->offset+312,8);
	    int ss;
	    get_bits(grib2_msg->buffer,&ss,grib2_msg->offset+320,8);
	    grib2_msg->md.stat_proc.etime=hh*10000+mm*100+ss;
/* number of time range specifications */
	    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.num_ranges,grib2_msg->offset+328,8);
/* number of values missing from process */
	    get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.nmiss,grib2_msg->offset+336,32);
	    if (grib2_msg->md.stat_proc.proc_code != NULL) {
		free(grib2_msg->md.stat_proc.proc_code);
		free(grib2_msg->md.stat_proc.incr_type);
		free(grib2_msg->md.stat_proc.time_unit);
		free(grib2_msg->md.stat_proc.time_length);
		free(grib2_msg->md.stat_proc.incr_unit);
		free(grib2_msg->md.stat_proc.incr_length);
		grib2_msg->md.stat_proc.proc_code=NULL;
	    }
	    grib2_msg->md.stat_proc.proc_code=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
	    grib2_msg->md.stat_proc.incr_type=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
	    grib2_msg->md.stat_proc.time_unit=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
	    grib2_msg->md.stat_proc.time_length=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
	    grib2_msg->md.stat_proc.incr_unit=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
	    grib2_msg->md.stat_proc.incr_length=(int *)malloc(grib2_msg->md.stat_proc.num_ranges*sizeof(int));
	    size_t off=368;
	    for (size_t n=0; n < grib2_msg->md.stat_proc.num_ranges; ++n) {
		get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.proc_code[n],grib2_msg->offset+off,8);
		get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_type[n],grib2_msg->offset+off+8,8);
		get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.time_unit[n],grib2_msg->offset+off+16,8);
		get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.time_length[n],grib2_msg->offset+off+24,32);
		get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_unit[n],grib2_msg->offset+off+56,8);
		get_bits(grib2_msg->buffer,&grib2_msg->md.stat_proc.incr_length[n],grib2_msg->offset+off+64,32);
		off+=96;
	    }
	    break;
	  }
	  case 15:
	  {
	    get_bits(grib2_msg->buffer,&grib2_msg->md.spatial_proc.stat_proc,grib2_msg->offset+272,8);
	    get_bits(grib2_msg->buffer,&grib2_msg->md.spatial_proc.type,grib2_msg->offset+280,8);
	    get_bits(grib2_msg->buffer,&grib2_msg->md.spatial_proc.num_points,grib2_msg->offset+288,8);
	    break;
	  }
	}
	break;
    }
    default:
    {
	fprintf(stderr,"Product Definition Template %d is not understood\n",grib2_msg->md.pds_templ_num);
	exit(1);
    }
  }
}

void unpack_DRS(GRIB2Message *grib2_msg)
{
  union {
    float fval;
    int ival;
  } u;

/* number of packed values */
  get_bits(grib2_msg->buffer,&grib2_msg->md.num_packed,grib2_msg->offset+40,32);
/* data representation template number */
  get_bits(grib2_msg->buffer,&grib2_msg->md.drs_templ_num,grib2_msg->offset+72,16);
  switch (grib2_msg->md.drs_templ_num) {
    case 0:
    case 3:
#ifdef JASPER
    case 40:
    case 40000:
#endif
    {
	get_bits(grib2_msg->buffer,(int *)&grib2_msg->md.R,grib2_msg->offset+88,32);
	int sign;
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+120,1);
	int value;
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+121,15);
	if (sign == 1) {
	  value=-value;
	}
	grib2_msg->md.E=value;
	get_bits(grib2_msg->buffer,&sign,grib2_msg->offset+136,1);
	get_bits(grib2_msg->buffer,&value,grib2_msg->offset+137,15);
	if (sign == 1) {
	  value=-value;
	}
	grib2_msg->md.D=value;
	grib2_msg->md.R/=pow(10.,grib2_msg->md.D);
	get_bits(grib2_msg->buffer,&grib2_msg->md.pack_width,grib2_msg->offset+152,8);
	get_bits(grib2_msg->buffer,&grib2_msg->md.orig_val_type,grib2_msg->offset+160,8);
	if (grib2_msg->md.drs_templ_num == 3) {
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.split_method,grib2_msg->offset+168,8);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.miss_val_mgmt,grib2_msg->offset+176,8);
	  if (grib2_msg->md.orig_val_type == 0) {
	    get_bits(grib2_msg->buffer,&u.ival,grib2_msg->offset+184,32);
	    grib2_msg->md.complex_pack.primary_miss_sub=u.fval;
	    get_bits(grib2_msg->buffer,&u.ival,grib2_msg->offset+216,32);
	    grib2_msg->md.complex_pack.secondary_miss_sub=u.fval;
	  }
	  else if (grib2_msg->md.orig_val_type == 1) {
	    get_bits(grib2_msg->buffer,&u.ival,grib2_msg->offset+184,32);
	    grib2_msg->md.complex_pack.primary_miss_sub=u.ival;
	    get_bits(grib2_msg->buffer,&u.ival,grib2_msg->offset+216,32);
	    grib2_msg->md.complex_pack.secondary_miss_sub=u.ival;
	  }
	  else {
	    fprintf(stderr,"Unable to decode missing value substitutes for original value type %d\n",grib2_msg->md.orig_val_type);
	    exit(1);
	  }
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.num_groups,grib2_msg->offset+248,32);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.width.ref,grib2_msg->offset+280,8);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.width.pack_width,grib2_msg->offset+288,8);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.length.ref,grib2_msg->offset+296,32);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.length.incr,grib2_msg->offset+328,8);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.length.last,grib2_msg->offset+336,32);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.length.pack_width,grib2_msg->offset+368,8);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.spatial_diff.order,grib2_msg->offset+376,8);
	  get_bits(grib2_msg->buffer,&grib2_msg->md.complex_pack.spatial_diff.order_vals_width,grib2_msg->offset+384,8);
	}
	break;
    }
    default:
    {
	fprintf(stderr,"Data template %d is not understood\n",grib2_msg->md.drs_templ_num);
	exit(1);
    }
  }
}

void unpack_BMS(GRIB2Message *grib2_msg)
{
/* bit map indicator */
  int ind;
  get_bits(grib2_msg->buffer,&ind,grib2_msg->offset+40,8);
  switch (ind) {
    case 0:
    {
	int len;
	get_bits(grib2_msg->buffer,&len,grib2_msg->offset,32);
	len=(len-6)*8;
	grib2_msg->md.bitmap=(unsigned char *)malloc(len*sizeof(unsigned char));
	for (size_t n=0; n < len; ++n) {
	  int bit;
	  get_bits(grib2_msg->buffer,&bit,grib2_msg->offset+48+n,1);
	  grib2_msg->md.bitmap[n]=bit;
	}
	break;
    }
    case 254:
    {
	break;
    }
    case 255:
    {
	grib2_msg->md.bitmap=NULL;
	break;
    }
    default:
    {
	fprintf(stderr,"This code is not currently set up to deal with predefined bit-maps\n");
	exit(1);
    }
  }
}

void unpack_DS(GRIB2Message *grib2_msg,int grid_num)
{
  float D=pow(10.,grib2_msg->md.D),E=pow(2.,grib2_msg->md.E);
  size_t off=grib2_msg->offset+40;
  switch (grib2_msg->md.drs_templ_num) {
    case 0:
    {
	size_t required_size=grib2_msg->md.ny*grib2_msg->md.nx;
	if (required_size > grib2_msg->grids[grid_num].gcapacity) {
	  if (grib2_msg->grids[grid_num].gridpoints != NULL) {
	    free(grib2_msg->grids[grid_num].gridpoints);
	  }
	  grib2_msg->grids[grid_num].gcapacity=required_size;
	  grib2_msg->grids[grid_num].gridpoints=(double *)malloc(grib2_msg->grids[grid_num].gcapacity*sizeof(double));
	}
	for (size_t n=0; n < grib2_msg->md.ny*grib2_msg->md.nx; ++n) {
	  if (grib2_msg->md.bitmap == NULL || grib2_msg->md.bitmap[n] == 1) {
	    int pval;
	    get_bits(grib2_msg->buffer,&pval,off,grib2_msg->md.pack_width);
	    grib2_msg->grids[grid_num].gridpoints[n]=grib2_msg->md.R+pval*E/D;
	    off+=grib2_msg->md.pack_width;
	  }
	  else {
	    grib2_msg->grids[grid_num].gridpoints[n]=GRIB_MISSING_VALUE;
	  }
	}
	break;
    }
    case 3:
    {
	struct {
	  int *ref_vals,*widths,*lengths;
	  int *first_vals,sign,omin;
	  long long miss_val,group_miss_val;
	  int max_length;
	} groups;
	size_t required_size=grib2_msg->md.ny*grib2_msg->md.nx;
	if (required_size > grib2_msg->grids[grid_num].gcapacity) {
	  if (grib2_msg->grids[grid_num].gridpoints != NULL) {
	    free(grib2_msg->grids[grid_num].gridpoints);
	  }
	  grib2_msg->grids[grid_num].gcapacity=required_size;
	  grib2_msg->grids[grid_num].gridpoints=(double *)malloc(grib2_msg->grids[grid_num].gcapacity*sizeof(double));
	}
	if (grib2_msg->md.complex_pack.num_groups > 0) {
	  if (grib2_msg->md.complex_pack.miss_val_mgmt > 0) {
	    groups.miss_val=pow(2.,grib2_msg->md.pack_width)-1;
	  }
	  else {
	    groups.miss_val=GRIB_MISSING_VALUE;
	  }
	  groups.first_vals=(int *)malloc(grib2_msg->md.complex_pack.spatial_diff.order*sizeof(int));
	  for (int n=0; n < grib2_msg->md.complex_pack.spatial_diff.order; ++n) {
	    get_bits(grib2_msg->buffer,&groups.first_vals[n],off,grib2_msg->md.complex_pack.spatial_diff.order_vals_width*8);
	    off+=grib2_msg->md.complex_pack.spatial_diff.order_vals_width*8;
	  }
	  get_bits(grib2_msg->buffer,&groups.sign,off,1);
	  get_bits(grib2_msg->buffer,&groups.omin,off+1,grib2_msg->md.complex_pack.spatial_diff.order_vals_width*8-1);
	  if (groups.sign == 1) {
	    groups.omin=-groups.omin;
	  }
	  off+=grib2_msg->md.complex_pack.spatial_diff.order_vals_width*8;
	  groups.ref_vals=(int *)malloc(grib2_msg->md.complex_pack.num_groups*sizeof(int));
	  for (int n=0; n < grib2_msg->md.complex_pack.num_groups; ++n) {
	    get_bits(grib2_msg->buffer,&groups.ref_vals[n],off,grib2_msg->md.pack_width);
	    off+=grib2_msg->md.pack_width;
	  }
	  int pad;
	  if ( (pad=(off % 8)) > 0) {
	    off+=8-pad;
	  }
	  groups.widths=(int *)malloc(grib2_msg->md.complex_pack.num_groups*sizeof(int));
	  for (int n=0; n < grib2_msg->md.complex_pack.num_groups; ++n) {
	    get_bits(grib2_msg->buffer,&groups.widths[n],off,grib2_msg->md.complex_pack.width.pack_width);
	    off+=grib2_msg->md.complex_pack.width.pack_width;
	  }
	  if ( (pad=(off % 8)) > 0) {
	    off+=8-pad;
	  }
	  groups.lengths=(int *)malloc(grib2_msg->md.complex_pack.num_groups*sizeof(int));
	  for (int n=0; n < grib2_msg->md.complex_pack.num_groups; ++n) {
	    get_bits(grib2_msg->buffer,&groups.lengths[n],off,grib2_msg->md.complex_pack.length.pack_width);
	    off+=grib2_msg->md.complex_pack.length.pack_width;
	  }
	  if ( (pad=(off % 8)) > 0) {
	    off+=8-pad;
	  }
	  groups.max_length=0;
	  size_t end=grib2_msg->md.complex_pack.num_groups-1;
	  for (int n=0; n < end; ++n) {
	    groups.lengths[n]=grib2_msg->md.complex_pack.length.ref+groups.lengths[n]*grib2_msg->md.complex_pack.length.incr;
	    if (groups.lengths[n] > groups.max_length) {
		groups.max_length=groups.lengths[n];
	    }
	  }
	  groups.lengths[end]=grib2_msg->md.complex_pack.length.last;
	  if (groups.lengths[end] > groups.max_length) {
	    groups.max_length=groups.lengths[end];
	  }
// unpack the field of differences
	  size_t gridpoint_index=0;
	  for (int n=0; n < grib2_msg->md.complex_pack.num_groups; ++n) {
	    if (groups.widths[n] > 0) {
		if (grib2_msg->md.complex_pack.miss_val_mgmt > 0) {
		  groups.group_miss_val=pow(2.,groups.widths[n])-1;
		}
		else {
		  groups.group_miss_val=GRIB_MISSING_VALUE;
		}
		for (size_t m=0; m < groups.lengths[n]; ++m) {
		  int pval;
		  get_bits(grib2_msg->buffer,&pval,off,groups.widths[n]);
		  off+=groups.widths[n];
		  if ((grib2_msg->md.bitmap != NULL && grib2_msg->md.bitmap[gridpoint_index] == 0) || pval == groups.group_miss_val) {
		    grib2_msg->grids[grid_num].gridpoints[gridpoint_index]=GRIB_MISSING_VALUE;
		  }
		  else {
		    grib2_msg->grids[grid_num].gridpoints[gridpoint_index]=pval+groups.ref_vals[n]+groups.omin;
		  }
		  ++gridpoint_index;
		}
	    }
	    else {
// constant group
		for (int m=0; m < groups.lengths[n]; ++m) {
		  if ((grib2_msg->md.bitmap != NULL && grib2_msg->md.bitmap[gridpoint_index] == 0) || groups.ref_vals[n] == groups.miss_val) {
		    grib2_msg->grids[grid_num].gridpoints[gridpoint_index]=GRIB_MISSING_VALUE;
		  }
		  else {
		    grib2_msg->grids[grid_num].gridpoints[gridpoint_index]=groups.ref_vals[n]+groups.omin;
		  }
		  ++gridpoint_index;
		}
	    }
	  }
	  for (; gridpoint_index < grib2_msg->md.nx*grib2_msg->md.ny; ++gridpoint_index) {
	    grib2_msg->grids[grid_num].gridpoints[gridpoint_index]=GRIB_MISSING_VALUE;
	  }
	  float lastgp;
	  for (int n=grib2_msg->md.complex_pack.spatial_diff.order-1; n > 0; --n) {
	    lastgp=groups.first_vals[n]-groups.first_vals[n-1];
	    int num_not_missing=0;
	    for (size_t l=0; l < grib2_msg->md.nx*grib2_msg->md.ny; ++l) {
		if (grib2_msg->grids[grid_num].gridpoints[l] != GRIB_MISSING_VALUE) {
		  if (num_not_missing >= grib2_msg->md.complex_pack.spatial_diff.order) {
		    grib2_msg->grids[grid_num].gridpoints[l]+=lastgp;
		    lastgp=grib2_msg->grids[grid_num].gridpoints[l];
		  }
		  ++num_not_missing;
		}
	    }
	  }
	  lastgp=0.;
	  size_t num_not_missing=0;
	  for (size_t l=0; l < grib2_msg->md.nx*grib2_msg->md.ny; ++l) {
	    if (grib2_msg->grids[grid_num].gridpoints[l] != GRIB_MISSING_VALUE) {
		if (num_not_missing < grib2_msg->md.complex_pack.spatial_diff.order) {
		  grib2_msg->grids[grid_num].gridpoints[l]=grib2_msg->md.R+groups.first_vals[num_not_missing]*E/D;
		  lastgp=grib2_msg->md.R*D/E+groups.first_vals[num_not_missing];
		}
		else {
		  lastgp+=grib2_msg->grids[grid_num].gridpoints[l];
		  grib2_msg->grids[grid_num].gridpoints[l]=lastgp*E/D;
		}
		++num_not_missing;
	    }
	  }
	  if (grib2_msg->md.complex_pack.spatial_diff.order > 0) {
	    free(groups.first_vals);
	  }
	  if (grib2_msg->md.complex_pack.num_groups > 0) {
	    free(groups.ref_vals);
	    free(groups.widths);
	    free(groups.lengths);
	  }
	}
	else {
	  for (int n=0; n < grib2_msg->md.ny*grib2_msg->md.nx; ++n) {
	    grib2_msg->grids[grid_num].gridpoints[n]=GRIB_MISSING_VALUE;
	  }
	}
	break;
    }
#ifdef JASPER
    case 40:
    case 40000:
    {
	int len;
	get_bits(grib2_msg->buffer,&len,grib2_msg->offset,32);
	len=len-5;
	int *jvals=(int *)malloc(grib2_msg->md.ny*grib2_msg->md.nx*sizeof(int));
	size_t required_size=grib2_msg->md.ny*grib2_msg->md.nx;
	if (required_size > grib2_msg->grids[grid_num].gcapacity) {
	  if (grib2_msg->grids[grid_num].gridpoints != NULL) {
	    free(grib2_msg->grids[grid_num].gridpoints);
	  }
	  grib2_msg->grids[grid_num].gcapacity=required_size;
	  grib2_msg->grids[grid_num].gridpoints=(double *)malloc(grib2_msg->grids[grid_num].gcapacity*sizeof(double));
	}
	if (len > 0) {
	  dec_jpeg2000((char *)&grib2_msg->buffer[grib2_msg->offset/8+5],len,jvals);
	}
	size_t cnt=0;
	for (size_t n=0; n < grib2_msg->md.ny*grib2_msg->md.nx; ++n) {
	  if (grib2_msg->md.bitmap == NULL || grib2_msg->md.bitmap[n] == 1) {
	    if (len == 0) {
		jvals[cnt]=0;
	    }
	    grib2_msg->grids[grid_num].gridpoints[n]=grib2_msg->md.R+jvals[cnt++]*E/D;
	  }
	  else {
	    grib2_msg->grids[grid_num].gridpoints[n]=GRIB_MISSING_VALUE;
	  }
	}
	free(jvals);
	break;
    }
#endif
  }
}

int unpackgrib2(FILE *fp,GRIB2Message *grib2_msg)
{
  int status;
  if ( (status=unpack_IS(fp,grib2_msg)) != 0) {
    return status;
  }
  unpack_IDS(grib2_msg);
/* find out how many grids are in this message */
  size_t off=grib2_msg->offset;
  while (strncmp(&((char *)grib2_msg->buffer)[off/8],"7777",4) != 0) {
    int len;
    get_bits(grib2_msg->buffer,&len,off,32);
    int sec_num;
    get_bits(grib2_msg->buffer,&sec_num,off+32,8);
    if (sec_num == 7) {
	++grib2_msg->num_grids;
    }
    off+=len*8;
  }
/* allocate space for the grids */
  if (grib2_msg->num_grids > grib2_msg->grid_capacity) {
    if (grib2_msg->grids != NULL) {
	for (size_t n=0; n < grib2_msg->grid_capacity; ++n) {
	  if (grib2_msg->grids[n].md.bitmap != NULL) {
	    free(grib2_msg->grids[n].md.bitmap);
	    grib2_msg->grids[n].md.bitmap=NULL;
	  }
	  free(grib2_msg->grids[n].gridpoints);
	}
	free(grib2_msg->grids);
	grib2_msg->grids=NULL;
    }
    grib2_msg->grid_capacity=grib2_msg->num_grids;
    grib2_msg->grids=(GRIB2Grid *)malloc(grib2_msg->grid_capacity*sizeof(GRIB2Grid));
    for (size_t n=0; n < grib2_msg->grid_capacity; ++n) {
	grib2_msg->grids[n].gridpoints=NULL;
	grib2_msg->grids[n].gcapacity=0;
    }
  }
/* now decode the message */
  int grid_num=0;
  while (strncmp(&((char *)grib2_msg->buffer)[grib2_msg->offset/8],"7777",4) != 0) {
    int len;
    get_bits(grib2_msg->buffer,&len,grib2_msg->offset,32);
    int sec_num;
    get_bits(grib2_msg->buffer,&sec_num,grib2_msg->offset+32,8);
    switch (sec_num) {
	case 2:
	{
	  unpack_LUS(grib2_msg);
	  break;
	}
	case 3:
	{
	  unpack_GDS(grib2_msg);
	  break;
	}
	case 4:
	{
	  unpack_PDS(grib2_msg);
	  break;
	}
	case 5:
	{
	  unpack_DRS(grib2_msg);
	  break;
	}
	case 6:
	{
	  unpack_BMS(grib2_msg);
	  break;
	}
	case 7:
	{
	  grib2_msg->grids[grid_num].md=grib2_msg->md;
	  unpack_DS(grib2_msg,grid_num);
	  ++grid_num;
	  break;
	}
    }
    grib2_msg->offset+=len*8;
  }
  return 0;
}
