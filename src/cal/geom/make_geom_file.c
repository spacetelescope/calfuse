/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:    make_geom_file
 *
 * Description: Make the geom(Y distortion) file:  copy the distortion
 *		image from a FITS file produced by Rich Robinson, and
 *		write it out in the format required by CALFUSE.
 *
 * Usage:	make_geom_file version
 *
 *		The program will read the files ydist_1a_2d.fit, 
 *		ydist_1b_2d.fit, etc and write the files geom1axxx.fit
 *		(where xxx is the version number specified on the command line).
 *
 *		The output file has an empty primary extension (containing
 *		just header keywords), an empty extension for X distortions
 *		(not presently used as such), and a full-size image containing
 *		the Y distortion as a function of position (FLOAT's).
 *
 * History:     10/19/00        jwk     Begin work.
 ******************************************************************************/

#include <stdio.h>
#include "calfuse.h"
#define EFFMJD 50000.0
#define NXMAX 16384
#define NYMAX 1024

static char CF_PRGM_ID[] = "make_geom_file";
static char CF_VER_NUM[] = "1.1";


int main(int argc, char *argv[])
{
    char *segment[]={"1A","1B","2A","2B"};
    char *segments[]={"1a","1b","2a","2b"};
    char extname[]="GEOMETRIC DISTORTION";
    char infile[80];
    char outfile[80];
    int	version;
    int   i, j, hdutype, status=0, nullval=0, anynull=0;
    int   ix, iy;
    int	ival;
    int   naxis;
    long  naxes[2];
    int   nxpix, nypix;
    float effmjd;
    float bscale, bzero, datamin, datamax;
    float *inbuf;
    fitsfile *geomfits, *infits;

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    if (argc != 2)
	{
	fprintf (stderr, "Usage: make_geom_file version\n");
	exit (-1);
	}

    version = (int) strtol (argv[1], (char **) NULL, 10);
    if ((version <0) || (version > 999))
	{
	fprintf (stderr, "Illegal version value: %s\n", argv[1]);
	exit (-1);
	}


for (i=0; i<4; i++) {  /* 1A,1B,2A,2B */
    /* open the y-distortion file, move to first extension (the one with
       the Y distortion image in it) */
    sprintf(infile, "ydist_%2.2s_2d.fit", segments[i]);
    FITS_open_file (&infits, infile, READONLY, &status);
    FITS_movabs_hdu (infits, 2, &hdutype, &status);
    FITS_read_key (infits, TINT, "NAXIS1", &nxpix, NULL, &status);
    FITS_read_key (infits, TINT, "NAXIS2", &nypix, NULL, &status);
    if (nxpix != NXMAX)
	{
	fprintf (stderr, "NAXIS1 = %d, not %d\n", nxpix, NXMAX);
	exit (-1);
	}

    if (nypix != NYMAX)
	{
	fprintf (stderr, "NAXIS2 = %d, not %d\n", nypix, NYMAX);
	exit (-1);
	}

    /* buffers to hold one row of data at a time */
    inbuf = (float *) calloc (nxpix, sizeof(float));

    /* create the output file */
    naxis = 0;
    naxes[0] = 0;
    naxes[1] = 0;

    sprintf(outfile,"geom%2.2s%03d.fit", segments[i], version);
    printf("%20.20s\n",outfile); 

    FITS_create_file (&geomfits,outfile,&status);
    FITS_create_img (geomfits, SHORT_IMG, naxis, naxes, &status);

    /* populate basic keywords */
    FITS_write_key (geomfits, TSTRING, "CALFTYPE", "GEOM", 
		   "Calibration file type", &status);

    FITS_write_key (geomfits,TINT,"CALFVERS",&version,
		   "Calibration file version", &status);

    FITS_write_key (geomfits,TSTRING,"DETECTOR",segment[i],
                   "detector (1A, 1B, 2A, 2B", &status);

    effmjd=EFFMJD;
    FITS_write_key (geomfits,TFLOAT,"EFFMJD",&effmjd,
		   "Date on which file should be applied (MJD)", &status);

    ival = 1;
    FITS_write_key (geomfits, TINT, "SPECBINX", &ival,
		"Binning in Detector X coordinate for HIST mode", &status);

    FITS_write_key (geomfits, TINT, "SPECBINY", &ival,
		"Binning in Detector Y coordinate for HIST mode", &status);

    FITS_write_date (geomfits, &status);

    FITS_write_key (geomfits,TSTRING,"AUTHOR","RICH ROBINSON",
		   "Author of file", &status);

    /* create empty first extension (HDU 2) for X distortions;
     * not used presently, but is kept for possible later use.
     */
    FITS_create_img (geomfits, SHORT_IMG, naxis, naxes, &status);
    FITS_write_comment (geomfits, "No x correction for geometric distortion",
	&status);

    /* now create extension to hold Y distortions */
    naxis = 2;
    naxes[0] = nxpix;
    naxes[1] = nypix;
    bzero = 0.;
    bscale = 0.01;
    datamin = 0.;
    datamax = 0.;

    /* use SHORT_IMG to force bitpix=16 */
    FITS_create_img (geomfits, SHORT_IMG, naxis, naxes, &status);
    FITS_write_key (geomfits, TFLOAT, "BSCALE", &bscale, "Scale factor",
	&status);
    FITS_write_key (geomfits, TFLOAT, "BZERO", &bzero, "Zero level",
	&status);

    /* now loop over rows of distortion image */
    for (iy=0; iy<nypix; iy++)
	{
	FITS_read_img (infits, TFLOAT, nxpix*iy+1, nxpix, &nullval,
		inbuf, &anynull, &status);
	for (ix=0; ix<nxpix; ix++)
		{
		if (inbuf[ix] < datamin)
			datamin = inbuf[ix];
		if (inbuf[ix] > datamax)
			datamax = inbuf[ix];
		}
	/* cfitsio will do conversion to SHORT using bscale, bzero */
	FITS_write_img (geomfits, TFLOAT, nxpix*iy+1, nxpix, inbuf, &status);
	}

    FITS_write_key (geomfits, TFLOAT, "DATAMIN", &datamin, "Minimum value",
	&status);
    FITS_write_key (geomfits, TFLOAT, "DATAMAX", &datamax, "Maximum value",
	&status);


    /* done with this segment */
    FITS_close_file(geomfits, &status);
    FITS_close_file(infits, &status);

} /* end of loop over segments */

    return 0;
}
