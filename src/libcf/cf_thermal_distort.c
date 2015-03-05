/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_thermal_distort(fitsfile *header, long nevents, 
 *				   float *xfarf, float *yfarf, float *weight,
 *				   unsigned char *locflags)
 *
 * Description: Corrects the X position of the event for thermal distortion
 *              using the locations of the stim pulses.  The correction is
 *              made only to stim pulses and events in the detector active 
 *		region.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nevents       The number of events
 *              float     *xfarf        An array of event X positions
 *              float     *yfarf        An array of event Y positions
 *              float     *weight       Number of photons per pixel
 *              unsigned char *locflags Location flags for each event
 *
 * Calls:       cf_estimate_drift_coefficients
 *              cf_get_stim_positions   These are static functions.
 *
 * Return:      0  on success
 *
 * History:     08/06/02   1.1    peb   Begin work
 *              11/11/02   1.2    peb   Corrected function description and
 *                                      added cf_timestamp after ELEC_COR
 *                                      check.
 *              11/12/02   1.3    peb   Sets left and right stim pulse flags.
 *              11/12/02   1.4    peb   Added check to move only events in
 *                                      stim-pulse and active region.
 *		03/11/03   1.5    wvd   Changed locflags to unsigned char
 *              05/20/03   1.6    rdr   Added cf_proc_check call
 *		07/29/03   1.8    wvd	If cf_proc_check fails, return errflg.
 *              02/27/04   1.9    rdr   Added weights to the calculation - 
 *                                      needed for histogram mode
 *              03/01/04   1.10   rdr   Corrected bug in procedure for finding
 *                                       stim pulses
 *              07/21/04   1.13   wvd	Delete variable i; unused.
 *              03/03/05   1.14   wvd	Read thermal coefficients from STIM_CAL
 *					file if STIM pulses are unavailable.
 *              04/06/05   1.15   wvd	Require at least 500 counts in a stim
 *					pulse to compute a centroid.
 *		10/05/07   1.16   bot	Corrected jmax in estimate_drift
 *
 ****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

#define STIMS (LOCATION_STIML|LOCATION_STIMR)

static int
cf_estimate_drift_coefficients (fitsfile *header, float *x0, float *x1,
	float *y0, float *y1)
{
    char	stimfile[FLEN_FILENAME];
    int		status=0;
    long	j, jmax;
    float	*x0_array=NULL, *x1_array=NULL, *y0_array=NULL, *y1_array=NULL;
    double	expstart, *mjd=NULL;
    fitsfile	*stimfits;

    /* Read header keywords. */
    FITS_read_key(header, TDOUBLE, "EXPSTART", &expstart, NULL, &status);
    FITS_read_key(header, TSTRING, "STIM_CAL", stimfile, NULL, &status);
    cf_verbose(3, "Drift-coefficient calibration file = %s", stimfile);

    /* Read drift-coefficients from the calibration file. */
    FITS_open_file(&stimfits, cf_cal_file(stimfile), READONLY, &status);
    FITS_movabs_hdu(stimfits, 2, NULL, &status);
    jmax = cf_read_col(stimfits, TDOUBLE, "MJD", (void *) &mjd);
    jmax = cf_read_col(stimfits, TFLOAT,  "X0", (void *) &x0_array);
    jmax = cf_read_col(stimfits, TFLOAT,  "X1", (void *) &x1_array);
    jmax = cf_read_col(stimfits, TFLOAT,  "Y0", (void *) &y0_array);
    jmax = cf_read_col(stimfits, TFLOAT,  "Y1", (void *) &y1_array);
    FITS_close_file(stimfits, &status);

    /* Find the coefficients appropriate for this exposure. */
    j = 0;
    while (j < jmax-1 && mjd[j+1] < expstart) j++;
    *x0 = x0_array[j];
    *x1 = x1_array[j];
    *y0 = y0_array[j];
    *y1 = y1_array[j];

    free(mjd);
    free(x0_array);
    free(x1_array);
    free(y0_array);
    free(y1_array);

    return status;
}


static int
cf_get_stim_positions(long nevents, float *xfarf, float *yfarf, float *weight, 
                   unsigned char *loc, float width, float *stimlx, float *stimly,
		   float *stimrx, float *stimry)
{
    int   stim_status=0;
    long  j;
    double lx=0., ly=0., rx=0., ry=0., nl=0., nr=0.;

    cf_verbose(3,"initial: stimlx = %0.1f, stimly=%0.1f, width=%0.1f ",
               *stimlx, *stimly, width) ;
    cf_verbose(3,"initial: stimrx = %0.1f, stimry=%0.1f, width=%0.1f ",
		*stimrx, *stimry, width) ;

    for (j=0; j<nevents; j++) {
	if      (xfarf[j] >= *stimlx-width && xfarf[j] <= *stimlx+width &&
		 yfarf[j] >= *stimly-width && yfarf[j] <= *stimly+width) {
	    lx += xfarf[j] * weight[j];
	    ly += yfarf[j] * weight[j];
	    nl += weight[j];
	    loc[j] |= LOCATION_STIML;
	}
	else if (xfarf[j] >= *stimrx-width && xfarf[j] <= *stimrx+width &&
		 yfarf[j] >= *stimry-width && yfarf[j] <= *stimry+width) {
	    rx += xfarf[j] * weight[j];
	    ry += yfarf[j] * weight[j];
	    nr +=  weight[j];
	    loc[j] |= LOCATION_STIMR;
	}
	else {
	    loc[j] &= ~STIMS;
	}
    }

    cf_verbose(3,"nl=%0.1f,  nr=%0.1f ",nl,nr) ;

    if (nl > 500) {
	*stimlx = lx/nl;
	*stimly = ly/nl;
    } else {
	cf_if_warning("Cannot determine left STIM position.  Estimating.");
	stim_status = 1;
    }
    if (nr > 500) {
	*stimrx = rx/nr;
	*stimry = ry/nr;
    } else {
	cf_if_warning("Cannot determine right STIM position.  Estimating.");
	stim_status = 1;
    }

    if (!stim_status) {
	cf_verbose(3,"after iteration: stimlx = %0.1f, stimly=%0.1f ", *stimlx, *stimly);
	cf_verbose(3,"after iteration: stimrx = %0.1f, stimry=%0.1f ", *stimrx, *stimry);
    }

    return stim_status;
}

int
cf_thermal_distort(fitsfile *header, long nevents, float *xfarf, 
		   float *yfarf, float *weight, unsigned char *locflags)
{
    char CF_PRGM_ID[] = "cf_thermal_distort";
    char CF_VER_NUM[] = "1.16";

    char  elecfile[FLEN_VALUE], detector[FLEN_VALUE];
    char  keyword[FLEN_KEYWORD];
    int   errflg=0, status=0, stim_status=0;
    long  j;
    float farflx, farfly, farfrx, farfry, x0, x1, y0, y1;
    float stimlx, stimly, stimrx, stimry, oldlx, oldly, oldrx, oldry;
    float width=128., tol=0.01;
    fitsfile *elecfits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    FITS_read_key(header, TSTRING, "DETECTOR", detector, NULL, &status);
    /*
     *  Get the name of the calibration file and open it
     */
    FITS_read_key(header, TSTRING, "ELEC_CAL", elecfile, NULL, &status);
    FITS_open_file(&elecfits, cf_cal_file(elecfile), READONLY, &status);
    /*
     *  Read stim locations in FARF frame.
     */
    FITS_read_key(elecfits, TFLOAT, "STIMWDTH", &width, NULL, &status);
    FITS_read_key(elecfits, TFLOAT, "STIMTOLR", &tol, NULL, &status);
    sprintf(keyword, "STIMLX%s", detector);
    FITS_read_key(elecfits, TFLOAT, keyword, &farflx, NULL, &status);
    sprintf(keyword, "STIMLY%s", detector);
    FITS_read_key(elecfits, TFLOAT, keyword, &farfly, NULL, &status);
    sprintf(keyword, "STIMRX%s", detector);
    FITS_read_key(elecfits, TFLOAT, keyword, &farfrx, NULL, &status);
    sprintf(keyword, "STIMRY%s", detector);
    FITS_read_key(elecfits, TFLOAT, keyword, &farfry, NULL, &status);
    FITS_close_file(elecfits, &status);
    /*
     *  Get stim positions (iteratively).
     */
    stimlx = farflx; stimly = farfly;
    stimrx = farfrx; stimry = farfry;
    do {
	oldlx = stimlx;	oldly = stimly;
	oldrx = stimrx;	oldry = stimry;
	width /= 2.;
	stim_status = cf_get_stim_positions(nevents, xfarf, yfarf, weight,
		locflags, width, &stimlx, &stimly, &stimrx, &stimry);
	if (!stim_status)
	    cf_verbose(3,"dstimlx=%.1f,  dstimly=%.1f, dstimrx=%.1f,"
	    "dstimry=%.1f ", fabs(oldlx-stimlx), fabs(oldly-stimly),
	    fabs(oldrx-stimrx), fabs(oldry-stimry) ) ;
    } while (fabs(oldlx-stimlx) > tol && fabs(oldly-stimly) > tol && 
	     fabs(oldrx-stimrx) > tol && fabs(oldry-stimry) > tol &&
	     !stim_status);
    /*
     *  Write stim pulse positions to file header.
     */
    if (!stim_status) {
	FITS_update_key(header, TFLOAT, "STIM_L_X", &stimlx, NULL, &status);
	FITS_update_key(header, TFLOAT, "STIM_L_Y", &stimly, NULL, &status);
	FITS_update_key(header, TFLOAT, "STIM_R_X", &stimrx, NULL, &status);
	FITS_update_key(header, TFLOAT, "STIM_R_Y", &stimry, NULL, &status);
    }
    /*
     *  Calculate drift coefficients from stim pulse positions
     */
    if (!stim_status) {
	x0 = (stimlx*farfrx - stimrx*farflx)/(stimlx - stimrx);
	x1 = (farflx - farfrx)/(stimlx - stimrx);
	y0 = 0.5*((farfly + farfry) - (stimly + stimry));
	y1 = 1.0;
    }
    /* 
     *  ...or estimate them from date of exposure.
     */
    else
	cf_estimate_drift_coefficients (header, &x0, &x1, &y0, &y1);
    cf_verbose(2, "Drift coefficients: %f %f, %f %f", x0, x1, y0, y1);
    /*
     *  Apply thermal correction to events in active area and
     *  stim-pulse regions.
     */
    for (j=0; j<nevents; j++) {
	if (!(locflags[j] & LOCATION_SHLD) || (locflags[j] & STIMS)) {
	    xfarf[j] = x1*xfarf[j] + x0;
	    yfarf[j] = y1*yfarf[j] + y0;
	}
    }
    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
