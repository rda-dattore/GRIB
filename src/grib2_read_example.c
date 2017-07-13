#include <stdio.h>
#include <stdlib.h>
#include "unpackgrib2.c"

/*  If the compiler complains about the "pow" function being an undefined
**  symbol, recompile like this:
**    % cc -std=c99 grib2_read_example.c -lm
*/

int main(int argc,char **argv)
{
  if (argc != 2) {
    fprintf(stderr,"usage: %s GRIB2_file_name\n",argv[0]);
    exit(1);
  }
  FILE *fp;
  if ( (fp=fopen(argv[1],"rb")) == NULL) {
    fprintf(stderr,"Error opening %s\n",argv[1]);
    exit(1);
  }
  GRIB2Message grib2_msg;
  initialize(&grib2_msg);
  int status;
  size_t nmsg=0;
  while ( (status=unpackgrib2(fp,&grib2_msg)) == 0) {
    ++nmsg;
/* print some header information */
    int hr=grib2_msg.time/10000;
    int min=(grib2_msg.time/100) % 100;
    int sec=grib2_msg.time % 10000;
    printf("Message Number: %d  GRIB Edition: %d  Discipline: %d  Table Version: %d-%d  Source ID: %d-%d  Date: %4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d  Number of Grids: %d\n",nmsg,grib2_msg.ed_num,grib2_msg.disc,grib2_msg.table_ver,grib2_msg.local_table_ver,grib2_msg.center_id,grib2_msg.sub_center_id,grib2_msg.yr,grib2_msg.mo,grib2_msg.dy,hr,min,sec,grib2_msg.num_grids);
    for (size_t l=0; l < grib2_msg.num_grids; ++l) {
	printf("  Grid #: %d  Type: %d  Parameter: %d/%d  Level: %d/%f  R: %f\n",l+1,grib2_msg.grids[l].md.gds_templ_num,grib2_msg.grids[l].md.param_cat,grib2_msg.grids[l].md.param_num,grib2_msg.grids[l].md.lvl1_type,grib2_msg.grids[l].md.lvl1,grib2_msg.grids[l].md.R);
	printf("    Definition:  Dimensions: %d x %d  ",grib2_msg.md.nx,grib2_msg.md.ny);
	if (grib2_msg.md.scan_mode == 0) {
	  printf("NW Corner (Lat/Lon): %f,%f",grib2_msg.md.slat,grib2_msg.md.slon);
	}
	else if (grib2_msg.md.scan_mode == 0x40) {
	  printf("NW Corner (Lat/Lon): %f,%f",grib2_msg.md.lats.elat,grib2_msg.md.slon);
	}
	switch (grib2_msg.md.gds_templ_num) {
	  case 0:
	  case 40:
	  {
	    if (grib2_msg.md.scan_mode == 0) {
		printf("  SE Corner (Lat/Lon): %f,%f",grib2_msg.md.lats.elat,grib2_msg.md.lons.elon);
	    }
	    else if (grib2_msg.md.scan_mode == 0x40) {
		printf("  SE Corner (Lat/Lon): %f,%f",grib2_msg.md.slat,grib2_msg.md.lons.elon);
	    }
	    if (grib2_msg.md.gds_templ_num == 0) {
	      printf("  Lat/Lon Resolution: %f,%f\n",grib2_msg.md.yinc.lainc,grib2_msg.md.xinc.loinc);
	    }
	    else if (grib2_msg.md.gds_templ_num == 40) {
	      printf("  Lat Circles %d, Lon Resolution: %f\n",grib2_msg.md.yinc.lainc,grib2_msg.md.xinc.loinc);
	    }
	    break;
	  }
	}
/* print out the gridpoints for the grids in the first message */
	if (nmsg == 1) {
	  size_t x=0;
	  for (size_t n=0; n < grib2_msg.md.ny; ++n) {
	    for (size_t m=0; m < grib2_msg.md.nx; ++m) {
		printf("(i,j)=(%d,%d)",m,n);
		if (grib2_msg.grids[l].gridpoints[x] == GRIB_MISSING_VALUE) {
		  printf(" value=MISSING\n");
		}
		else {
		  printf(" value=%f\n",grib2_msg.grids[l].gridpoints[x]);
		}
		++x;
	    }
	  }
	}
    }
  }
  if (status == -1) {
    printf("EOF - end of file found\n");
  }
  else {
    printf("Read error after %d messages\n",nmsg);
  }
}
