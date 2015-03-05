/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:    saa output_file
 *
 * Description: Writes a FITS ASCII table with the SAA contours.
 *              
 *
 * History:     07/21/98  E. Murphy   Begin work.
 *
 * References:  SAA contours from Scott Friedman 07/16/98, same contours as
 *              used in SPIKE.  Note that the contours do not extend below
 *              a latitude of 30 degrees.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "calfuse.h"
 
#define VERSION 2
#define EFFMJD 50000.0

int main(int argc, char *argv[])

{

/* Define variables. */
int i, retcode, status=0, anynull, hdutype;
int numlimits=10;
double maxlat, maxlon, minlon, delta;
int nrows=10;
int tfields=3;
float saa_lat[10]=   {-30.0,-25.0, -20.0, -15.0,-10.0,  0.0,  5.0,  10.0, 12.0, 13.0};
float saa_minlon[10]={-91.0,-98.0,-102.0,-100.0,-93.0,-68.3,-56.0, -42.0,-32.0,-27.0};
float saa_maxlon[10]={ 24.1, 47.0,  41.0,  33.0, 25.0,  8.0, -3.0, -14.0,-20.0,-27.0};
float floatnull;
fitsfile *outfits, *infits;
long dumar[2]={0,0};
char extname[]="SAA_CONTOURS";
char *ttype[]={"LATITUDE","MIN_LONGITUDE","MAX_LONGITUDE"};
char *tform[]={"F8.3","F8.3","F8.3"};
char *tunit[]={"deg","deg","deg"};
int  vers;
float effmjd;

FITS_create_file(&outfits, "saac002.fit", &status);

FITS_create_img(outfits, SHORT_IMG, 0, dumar, &status);

     FITS_write_key(outfits, TSTRING, "CALFTYPE", "SAAC", 
                    "Calibration file type", &status);

     vers=VERSION;
     FITS_write_key(outfits,TINT,"CALFVERS",&vers,
                    "Calibration file version", &status);

     effmjd=EFFMJD;
     FITS_write_key(outfits,TFLOAT,"EFFMJD",&effmjd,
                    "Date on which file should be applied (MJD)", &status);

     FITS_write_date(outfits, &status);

     FITS_write_key(outfits,TSTRING,"AUTHOR","EDWARD MURPHY",
                    "Author of file", &status);

FITS_create_tbl(outfits, ASCII_TBL, nrows, tfields, ttype, tform, tunit, 
extname, &status);

FITS_write_col(outfits, TFLOAT, 1, 1, 1, nrows, saa_lat, &status);
FITS_write_col(outfits, TFLOAT, 2, 1, 1, nrows, saa_minlon, &status);
FITS_write_col(outfits, TFLOAT, 3, 1, 1, nrows, saa_maxlon, &status);

FITS_close_file(outfits, &status);

return 0;

}
