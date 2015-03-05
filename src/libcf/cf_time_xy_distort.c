/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_time_xy_distort(fitsfile *header, long nevents,
 *		float *xfarf, float *yfarf, unsigned char *locflags)
 *
 * Description: Corrects each photon event in the active area for long-term
 *		drifts in the X and Y coordinate scales.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nevents       The number of events
 *              float     *xfarf        An array of event X positions
 *              float     *yfarf        An array of event Y positions
 *              unsigned char *locflags Location flags for each event
 *
 * Calls:
 *
 * Return:      0  on success
 *
 * History:     11/01/06   1.1   wvd    Based on cf_count_rate_y_distort.c
 *		11/22/06   1.2   wvd	Read BIN_FACT from file header.
 *					Make X correction a function of X.
 *		01/19/07   1.3   wvd	Clean up i/o.
 *		04/07/07   1.4   wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_time_xy_distort(fitsfile *header, long nevents,
	float *xfarf, float *yfarf, unsigned char *locflags)
{
    char CF_PRGM_ID[] = "cf_time_xy_distort";
    char CF_VER_NUM[] = "1.4";
    
    char  tmxyfile[FLEN_VALUE];
    int   errflg=0, hdutype, status=0, anynull=0;
    int   binfact, i, ii, i_max, tndx, xndx, yndx;
    long  j, xlen, ylen;
    double expstart;
    float *xstretch2d, *ystretch2d, xstretch1d[NXMAX], ystretch1d[NYMAX];
    fitsfile *tmxyfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /*
     *  Read MJD of exposure start.
     */
    FITS_read_key(header, TDOUBLE, "EXPSTART", &expstart, NULL, &status);

    /*
     *  Read the X correction array.
     *  Note: the first column of the calibration array contains MJD.
     */
    FITS_read_key(header, TSTRING, "TMXY_CAL", tmxyfile, NULL, &status);
    cf_verbose(3, "Reading TMXY_CAL file %s", tmxyfile);
    FITS_open_file(&tmxyfits, cf_cal_file(tmxyfile), READONLY, &status);
    FITS_movabs_hdu(tmxyfits, 2, &hdutype, &status);
    FITS_read_key(tmxyfits, TLONG, "NAXIS1", &xlen, NULL, &status);
    FITS_read_key(tmxyfits, TLONG, "NAXIS2", &ylen, NULL, &status);
    FITS_read_key(tmxyfits, TINT, "BIN_FACT", &binfact, NULL, &status);
    xstretch2d = (float *) cf_malloc(sizeof(float)*xlen*ylen);
    FITS_read_img(tmxyfits, TFLOAT, 1L, xlen*ylen, NULL, xstretch2d,
    	&anynull, &status);

    /* 
     * Find the row in the calibration array that corresponds to the 
     * exposure start time.
     */
    i = 0;
    while (expstart > *(xstretch2d+i*xlen) && i < ylen) i++;
    tndx = i-1;
    if (tndx < 0) tndx = 0;
    cf_verbose(2,"EXPSTART = %g.", expstart);
    cf_verbose(2,"Using X corrections for MJD >= %g", *(xstretch2d+tndx*xlen));
    if (tndx < ylen-1)
	cf_verbose(2,"                    and MJD <  %g",
	*(xstretch2d+(tndx+1)*xlen));

    /* 
     *  Now compute a 1-D array containing corrections to X
     *  coordinate as a function of X position on the detector.
     *  Shift information binned by BIN_FACT pixels.
     *  We interpolate to the nearest X pixel.  Otherwise,
     *  counts pile up at the boundaries.
     */
    i_max = (xlen-2) * binfact;
    for (i = 0; i < i_max; i++) {
	ii = (i/binfact) * binfact;
	xstretch1d[i] = ((ii+binfact-i) * (xstretch2d+tndx*xlen)[1+i/binfact] +
	                (i-ii) * (xstretch2d+tndx*xlen)[2+i/binfact]) / binfact;
    }
    for ( ; i < NXMAX; i++) 
	xstretch1d[i] = (xstretch2d+tndx*xlen)[xlen-1];

    /*
     *  Read the Y correction array.
     */
    FITS_movabs_hdu(tmxyfits, 3, &hdutype, &status);
    FITS_read_key(tmxyfits, TLONG, "NAXIS1", &xlen, NULL, &status);
    FITS_read_key(tmxyfits, TLONG, "NAXIS2", &ylen, NULL, &status);
    FITS_read_key(tmxyfits, TINT, "BIN_FACT", &binfact, NULL, &status);
    ystretch2d = (float *) cf_malloc(sizeof(float)*xlen*ylen);
    FITS_read_img(tmxyfits, TFLOAT, 1L, xlen*ylen, NULL, ystretch2d,
    	&anynull, &status);
    FITS_close_file(tmxyfits, &status);

    /* 
     * Find the row in the calibration array that corresponds to the 
     * exposure start time.
     */
    i = 0;
    while (expstart > *(ystretch2d+i*xlen) && i < ylen) i++;
    tndx = i-1;
    if (tndx < 0) tndx = 0;
    cf_verbose(3,"Using Y corrections for MJD >= %g", *(ystretch2d+tndx*xlen));
    if (tndx < ylen-1)
	cf_verbose(3,"                    and MJD <  %g",
	*(ystretch2d+(tndx+1)*xlen));

    /* 
     *  Now compute a 1-D array containing corrections to Y 
     *  coordinate as a function of Y position on the detector.
     */
    i_max = (xlen-2) * binfact;
    for (i = 0; i < i_max; i++) {
	ii = (i/binfact) * binfact;
	ystretch1d[i] = ((ii+binfact-i) * (ystretch2d+tndx*xlen)[1+i/binfact] +
	                (i-ii) * (ystretch2d+tndx*xlen)[2+i/binfact]) / binfact;
    }
    for ( ; i < NYMAX; i++) 
	ystretch1d[i] = (ystretch2d+tndx*xlen)[xlen-1];
    cf_verbose(3, "For y = 430, ystretch = %f\n", ystretch1d[430]);
    cf_verbose(3, "For y = 810, ystretch = %f\n", ystretch1d[810]);

    /* 
     * Loop through all events. Move only events in active region.
     */
    for (j = 0; j < nevents; j++) {
	if (!(locflags[j] & LOCATION_SHLD)) {
	    xndx = (int) xfarf[j];
	    if (xndx < 0) xndx = 0;
	    else if (xndx >= NXMAX) xndx = NXMAX-1;

	    yndx = (int) yfarf[j];
	    if (yndx < 0) yndx = 0;
	    else if (yndx >= NYMAX) yndx = NYMAX-1;

	    xfarf[j] += xstretch1d[xndx];
	    yfarf[j] += ystretch1d[yndx];
	}
    }

    free(xstretch2d);
    free(ystretch2d);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
