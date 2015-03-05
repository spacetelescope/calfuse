/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_count_rate_y_distort(fitsfile *header, long nevents,
 *                                      float *time, float *yfarf,
 *                                      unsigned char *locflags, long nseconds,
 *                                      float *timeline, float *fec_rate)
 *
 * Description: Corrects the Y position of each event based on the
 *              instantaneous count rate.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nevents       The number of events
 *              float     *time         An array of event times
 *              float     *yfarf        An array of event Y positions
 *              unsigned char *locflags Location flags for each event
 *              long      nseconds      The number of seconds in the timeline
 *              float     *timeline     The timeline in seconds
 *              float     *fec_rate     Front End Counter rate
 *
 * Calls:
 *
 * Return:      0  on success
 *
 * History:     08/06/02   1.1   peb    Begin work
 *              10/27/02   1.2   peb    Added a time-dependent correction
 *                                      which uses the Front End Counter
 *                                      (FEC) rates.
 *              11/11/02   1.3   peb    Removed non-time-dependent correction.
 *                                      (Code now expects timeline extension
 *                                      data.)
 *              11/11/02   1.4   peb    Fixed compile error - added buffer
 *                                      variable.
 *              11/12/02   1.5   peb    Added check to move only events in
 *                                      active region.
 *		03/11/03   1.6   wvd    Changed locflags to unsigned char
 *              05/20/03   1.7   rdr    Added call to cf_proc_check
 *		07/29/03   1.8   wvd    If cf_proc_check fails, return errflg.
 *		08/04/03   1.9   wvd    Convert fec_rate to type short.
 *		11/26/03   1.10  wvd    Convert fec_rate to type float.
 *		08/18/04   1.11  wvd    Interpolate stretch correction among
 *					tabulated Y values.
 *		10/06/05   1.12  wvd    For HIST data, set each element of
 *					fec_rate to weighted mean of array.
 *		04/07/07   1.13  wvd    Clean up compiler warnings.
 *
 ****************************************************************************/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_count_rate_y_distort(fitsfile *header, long nevents, float *time,
			float *yfarf, unsigned char *locflags, long nseconds,
			float *timeline, float *fec_rate)
{
    char CF_PRGM_ID[] = "cf_count_rate_y_distort";
    char CF_VER_NUM[] = "1.13";
    
    char  instmode[FLEN_VALUE], ystrfile[FLEN_VALUE];
    short xlen, ylen;
    int   errflg=0, hdutype, status=0, anynull=0;
    int   i, ii, i_max, rndx, yndx;
    long  j, k;
    float *ystretch, ystr_1d[NYMAX], flt=0.;
    fitsfile *ystrfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /*
     *  For HIST data, replace fec_rate with weighted mean value.
     */
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    if (!strncmp(instmode, "HIST", 4)) {
	float *mean_rate, ratio, numerator=0, denominator=0;
	mean_rate = (float *) cf_malloc(sizeof(float) * nseconds);
	for (i = 0; i < nseconds; i++) {
	    numerator += fec_rate[i] * fec_rate[i];
	    denominator += fec_rate[i];
	}
	ratio = numerator / denominator;
	for (i = 0; i < nseconds; i++) mean_rate[i] = ratio;
	fec_rate = mean_rate;
	cf_verbose(2, "Weighted mean count rate = %f", ratio);
    }

    /*
     *  Read the rate calibration file.
     */
    FITS_read_key(header, TSTRING, "RATE_CAL", ystrfile, NULL, &status);
    FITS_open_file(&ystrfits, cf_cal_file(ystrfile), READONLY, &status);
    FITS_movabs_hdu(ystrfits, 2, &hdutype, &status);
    FITS_read_key(ystrfits, TSHORT, "NAXIS1", &xlen, 0, &status);
    FITS_read_key(ystrfits, TSHORT, "NAXIS2", &ylen, 0, &status);
    ystretch = (float *) cf_malloc(sizeof(float)*xlen*ylen);
    FITS_read_img(ystrfits, TFLOAT, 1L, xlen*ylen, &flt, ystretch, &anynull,
		  &status);
    FITS_close_file(ystrfits, &status);

    i_max = (xlen-1) * 10;
    for(j=k=0; j<nevents; j++) {	/* Loop through all events. */
	/*
	 *  Move only events in active region.
	 */
	if (!(locflags[j] & LOCATION_SHLD)) {

	    /* If time has changed, compute a new ystretch array. */
	    while(timeline[k]-FRAME_TOLERANCE < time[j] && k < nseconds) {
		k++;

		/* This index reflects the count rate. */
	        rndx = (int)((log10(fec_rate[k])-1.)*10.);
	        if (rndx < 0)
		    rndx = 0;
	        else if (rndx >= ylen)
		    rndx = ylen-1;

		/* 
		 *  This index reflects Y position on the detector.
		 *  Shift information is provided for every 10 Y pixels.
		 *  We interpolate to the nearest Y pixel.  Otherwise,
		 *  counts pile up at the 10-pixel boundaries.
		 */
		for (i = 0; i < i_max; i++) {
		    ii = (i/10) * 10;
		    ystr_1d[i] = ((ii+10-i) * (ystretch+rndx*xlen)[i/10] +
				  (i - ii)  * (ystretch+rndx*xlen)[i/10+1]) / 10.;
		}
		for ( ; i < NYMAX; i++) 
		    ystr_1d[i] = (ystretch+rndx*xlen)[xlen-1];
	    }

	    /* Now apply the shift. */
	    yndx = (int) yfarf[j];
	    if (yndx < 0)
		yndx = 0;
	    else if (yndx >= NYMAX)
		yndx = NYMAX-1;
	    yfarf[j] -= ystr_1d[yndx];
	}
    }
    free(ystretch);
    if (!strncmp(instmode, "HIST", 4)) free(fec_rate);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
