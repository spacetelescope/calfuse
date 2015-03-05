/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences  
 *              FUSE
 *****************************************************************************
 *
 * Synopsis: cf_scale_bkgd(infits, nevents, x_in, y_in, w_in, channel,
 *		timeflags, locflags, ngood, good_index, dtime, ntime, 
 *		&binx, &biny, &nx_lif, &ny_lif,
 *              &ymin_lif, img_lif, &nx_sic, &ny_sic, &ymin_sic, &img_sic)
 *
 * Arguments:  *fitsfile     infits      Input IDF file
 *             long          nevents     Number of photon events recorded
 *             float         x_in, y_in  Coordinates of each photon
 *	       float	     w_in	 Weighting factor for each photon
 *             unsigned char channel     Aperture ID for each photon
 *             unsigned char timeflags   Photon status flags regarding time
 *             unsigned char locflags    Photon status flags regarding location
 *             long          ngood       Number of good photons
 *             long          good_index  List of indices for the good photons
 *             long          dtime       Length of the daytime exposure
 *             long          ntime       Length of the night time exposure
 *             int           binx, biny  Binning factors of the background
 *                                         array along the X and Y axes
 *             int         nx_lif, nx_sic X Size of background array
 *             int         ny_lif, ny_sic Y Size of background array
 *             int     ymin_lif, ymin_sic Y value corresponding to the
 *                                          first row of the bkgd image
 *             float   *img_lif, *img_sic 2D background image (binned)
 *
 * Description:
 *  (a) Determine the total counts in a background sample region for data
 *      taken during the night and during the day. Use the data in the
 *      IDF event list, after screening out the geocoronal lines.
 *  (b) Read estimated intrinsic count rates from the background characterization
 *      file and estimate the intrinsic counts in the background region from
 *      the exposure times and number of pixels in the region.
 *  (c) Subtract the estimated intrinsic counts from the measured total to 
 *      get an estimate of the scattered light contribution for both day
 *      and night portions of the observation
 *  (d) Use the scattered light counts to scale daytime and night background 
 *      images obtained from the background calibration files
 *  (e) Combine the two background images to obtain the final background model
 *
 * A previous assumption that the intrinsic background is constant was
 * found to be incorrect. The program now iterates the assumed intrinsic 
 * background count rate (assumed constant in day and night) and checks the
 * calculated model against the observations by comparing the intensity 
 * variations along the x direction. In this comparison, we mask out regions
 * that are affected by geocoronal contamination.
 *
 * Caveats: The "intrinsic" component of the background is neither intrinsic
 * nor constant in day and night.  It is better called the "spatially uniform"
 * component and should be computed independently for day and night portions
 * of the exposure.  wvd 01/24/2006
 *
 * History:     01/11/03   v1.1   rdr   Started work - adapted program
 *                                              from cf_make_ttag_bkgd
 *              03/01/03   v1.3   wvd   Correct use of pointer in FITS_read_key
 *		03/11/03   v1.4   wvd   Update documentation regarding use of 
 *					unsigned char.
 *              03/28/03   v1.5   rdr   corrected error in freeing storage
 *              03/31/03   v1.6   rdr   corrected error in assigning subarrays
 *                                      for lif and sic apertures
 *              05/08/03   v1.7   rdr   read limits to the background region
 *                                      from the header of the IDF. Also 
 *                                      changed print statements to cf_verbose.
 *              05/20/03   v1.8   rdr   - Make a model using only exposure time
 *                                      if no background sample regions have
 *                                      been specified (generally for HIST data)
 *                                      - Look into parm file to determine
*                                       whether to generate a background model
 *              05/22/03   v1.9    rdr  Correct problems occuring when exptime=0
 *              06/04/03   v1.10   rdr  Change method for determining limits of 
 *                                      output arrays
 *              09/08/03   v1.12   wvd  Use cf_read_col in get_limits()
 *              09/29/03   v1.13   wvd  Fix bug in population of 2-D bkgd
 *					array; clean up i/o.
 *              10/07/03   v1.14   wvd  Reject airglow photons using a
 *					bit-wise comparison of locflags.
 *              10/31/03   v1.15   wvd  Increase aperture pad from 10 to 15.
 *					Change channel to unsigned char.
 *              10/31/03   v1.16   wvd  Use cf_read_col in get_geocoronal_mask()
 *              03/22/04   v1.17   wvd  Change get_limits so that sizes of
 *					bkgd and probability arrays are equal.
 *              04/08/04   v1.18   bjg  Replace NXMAX by nx
 *                                      Remove unused variables
 *                                      Change format in printf to match 
 *                                      arg type.
 *		04/27/04   v1.19   wvd	Replace RUN_MKBK with RUN_BKGD
 *					throughout.
 *		07/01/04   v1.20   wvd	Fix bug: nxb was undefined for
 *					RUN_BKGD = NO.
 *		03/08/05   v1.21   wvd	Convert to cf_scale_bkgd.
 *					Treat data and models the same way.
 *					Modify get_xvar_img.  Add
 *					make_bad_pixel_mask and sum_good_pixels.
 *		06/15/05   v1.22   wvd	BUG FIX: get_limits always read the
 *					probability arrays for point-source 
 *					targets from the WGTS_CAL file.  Now
 *					it distinguishes between point-source
 *					and extended targets.
 *		01/20/06   v1.23   wvd	BUG FIX: Apply uniform X limits to
 *					data and background models.
 *					BUG FIX: Allow output bkgd images
 *					to fall off the detector if target 
 *					centroid requires it.
 *		01/24/06   v1.24   wvd	BUG FIX: Test y coordinate of photon
 *					before adding it to the data image.
 *		02/01/06   v1.25   wvd	BUG FIX: Make npix a long in 
 *					sum_good_pixels.  Correct error in
 *					extraction of sub-arrays.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calfuse.h"


/*********************************************************************
 *
 *    GET_XVAR_IMG
 *
 *   Subroutine to compute the intensity variation along the X direction
 *   of a background image.  Pixels flagged as bad are excluded.
 *     bkgimg = image from which to extract the data
 *     maskimg = mask for data image (1 = good pixel)
 *     nxb = size of bkgimg in x
 *     nxsam = number of elements in the xvars array
 *     xvars = array of x variations
 *
 ***********************************************************************/


int
get_xvar_img(float *bkgimg, float *maskimg, int nxb, int nxsam, float *xvars)
{
    int     nbin;
    long    i, j;

    nbin = nxb/nxsam ;
    for (i = 0; i < NYMAX; i++) 
	for (j = 0; j < nxb; j++)
	    xvars[j/nbin] += bkgimg[i*nxb+j] * maskimg[i*nxb+j];

    return 0 ;
}

/*********************************************************************
 *
 *       MAKE_MODEL
 *
 *  Subroutine to combine a day and night calibration image to generate the 
 *    final background model:
 *       npix = number of pixels in the full background image
 *       nsam = number of pixels in the background regions of the detector
 *       nxb = size of the x dimension of the background images
 *       intcr = assumed intrinsic count rate
 *       expnight, expday = exposure times during the night and day
 *       data_sum_night, data_sum_day = sum of weights within the background
 *             sample regions for night and day exposures
 *       model_sum_night, model_sum_day = number of counts in the unscaled
 *              night and day background calibration images
 *		- background regions only
 *       bkgimgn = image containing the night calibration image.
 *       bkgimgd = image containing the daytime calibration image 
 *       bkgimg = final background image 
 *	 maskimg = mask defining good background regions
 *
 *****************************************************************************/

int
make_model(long npix, long nsam, int nxb, float intcr, float expnight,
    float expday, double data_sum_night, double data_sum_day,
    double model_sum_night, double model_sum_day,
    float *bkgimgn, float *bkgimgd, float *bkgimg, float *maskimg)
{
    int   bin_size ;
    long  i;
    float int_counts;
    float sfd, sfn, intnight, intday, scatnight, scatday ;

    /* bin_size = binning factor of the background image */
    bin_size = NXMAX / nxb ;

    cf_verbose(3, "\n\tConstructing the background model\n") ;
    cf_verbose(3, "total counts in the data: night= %.1f, day = %.1f",
        data_sum_night, data_sum_day);
    cf_verbose(3, "intrinsic count rate = %.2f x 10^-7 ",intcr ) ;
    cf_verbose(3, "exptime night = %.1f, exptime day = %.1f ",
        expnight, expday );

    /* Sum the intrinsic background counts */
    intnight = nsam * bin_size * intcr * expnight * 1.e-7 ;
    intday = nsam * bin_size * intcr * expday * 1.e-7 ;
    cf_verbose(3, "total intrinsic counts: night= %.1f , day= %.1f ",
        intnight, intday) ;

    /* Sum of scattered light in data (measured - intrinsic) */
    scatnight = data_sum_night - intnight ;
        if (scatnight <= 0.) scatnight = 0. ;
    scatday = data_sum_day - intday ;
        if (scatday <= 0.) scatday = 0. ;
    cf_verbose(3, "total scattered light count: night = %.1f, day = %.1f",
        scatnight, scatday);
    cf_verbose(3, "total counts in unscaled background model: "
        "night = %.2f, day = %.2f", model_sum_night, model_sum_day);

    /* Calculate the scale factor for the night scattered light image */
    sfn = 0. ;
    if (model_sum_night > 0)  sfn = scatnight / model_sum_night ;

    /* Calculate the scale factor for the day scattered light image */
    sfd = 0. ;
    if (model_sum_day > 0) sfd = scatday / model_sum_day ;
    cf_verbose(3, "scale factors for background images: "
        "night = %.1f, day = %.1f ", sfn, sfd) ;

    /* Calculate the intrinsic counts per binned pixel */
    int_counts = bin_size * intcr * (expnight + expday) * 1.e-7 ;
    cf_verbose(3, "intrinsic counts per binned pixel = %.6f ", int_counts ) ;

    /* Use the scattered light scaling factors and the tabulated intrinsic
       count rate to model the final background image */

    for (i = 0; i < npix; i++)
        bkgimg[i] = (sfn * bkgimgn[i]) + (sfd * bkgimgd[i]) + int_counts;

    return 0 ;
}


/*********************************************************************
 *
 *     MAKE_MODEL_EXPTIME  : 
 *
 *  Subroutine to combine a day and night calibration image based on
 *    exposure times only to generate the final background model:
 *      npix = total number of pixels in the image
 *      nxb = length of the x axis in the background images
 *     intcr = assumed intrinsic count rate
 *     expnight, expday = exposure times during the night and day
 *     bkging = image containing the night calibration image.
 *     bkgimgd = image containing the daytime calibration image 
 *     bkgimg = final background model    
 *
 *********************************************************************/
 

int
make_model_exptime(long npix, int nxb, float intcr, float expnight,
    float expday, float *bkgimgn, float *bkgimgd, float *bkgimg)
{

    int i, bin_size ;
    float intnight, intday ;

    /* Compute the binning factor in x */
    bin_size = NXMAX/nxb ;

    /* Calculate the expected intrinsic counts in the background */ 
    intnight = bin_size * intcr * expnight * 1.e-7 ;
    intday = bin_size * intcr * expday * 1.e-7 ;

    /* Combine the day and night images based on exposure times */
    for (i=0 ; i<npix ; i++)
    bkgimg[i] = intnight + intday + bkgimgn[i]*expnight + bkgimgd[i]*expday ;

   return 0 ;
}

/*****************************************************************************
 *
 *    MAKE_BAD_PIXEL_MASK :
 *
 *  Produces a mask the same size as the background image. 
 *  Good pixels are flagged with 1 and bad pixels with 0.  Bad pixels are
 *  - all pixels that are NOT in the user-selected background regions;
 *  - all pixels containing airglow lines (from AIRG_CAL); and
 *  - all pixels lying outside of the active area of the detector.
 *	bkgd_num = number of background regions
 *	ylim - limits of background regions
 *      npix = number of pixels in the background array
 *      nxb = dimension of the x axis of the array
 *      maskimg = array to contain the mask (must be initialized to zero)
 *      nsam = number of good pixels in the mask
 *
 ****************************************************************************/

int
make_bad_pixel_mask (fitsfile *infits, int bkgd_num, short *ylim,
	long npix, int nxb, float *maskimg, long *nsam)
{

    int nbin;
    int xmin,xmax,ymin,ymax ;
    long nrow, ndx, i, j, k, n_air, n_good;
    short *gxmin, *gxmax, *gymin, *gymax ;

    /* Compute the binning factor in x */
    nbin = NXMAX/nxb;

    /* Ignore regions near the detector edges. */
    xmin =  2048 / nbin;
    xmax = 14336 / nbin;

    /* Flag all pixels in the user-defined background regions as good (1). */
    for (i = 0; i < bkgd_num; i++)
	for (j = ylim[2*i]; j <= ylim[2*i+1]; j++)
	    for (k = xmin; k < xmax; k++) {
		ndx = (long) (j*nxb+k) ;
		maskimg[ndx] = 1.;
	    }

    /* Read the positions of the geocoronal lines. */  
    nrow = cf_get_geocorona(infits, &gxmin, &gxmax, &gymin, &gymax);
    cf_verbose(3, "Number of geocoronal lines = %ld", nrow) ;

    /* Flag pixels affected by geocoronal lines as bad (0). */
    n_air = 0;
    for (i=0; i<nrow; i++) {
	xmin = gxmin[i] / nbin ;
	xmax = gxmax[i] / nbin ;
	ymin = gymin[i] ;
	ymax = gymax[i] ;
	for (j=ymin; j<ymax; j++) 
	    for (k=xmin; k<xmax; k++) {
		ndx = j * nxb + k ;
		if (maskimg[ndx] > 0) {
		    maskimg[ndx] = 0. ;
		    n_air++ ;
		}
	    }
	}

    free(gxmin);
    free(gxmax);
    free(gymin);
    free(gymax);

    cf_verbose(3, "Number of bkgd pixels affected by geocoronal lines = %ld",
	n_air) ;

    /* Count up and return the number of good pixels. */
    n_good = 0;
    for (i = 0; i < npix; i++) if(maskimg[i] > 0) n_good++;
    cf_verbose(3, "Number of bkgd pixels flagged as good = %ld", n_good) ;

    *nsam = n_good;
    return 0 ;
}


/******************************************************************************
 *
 *             GET_LIMITS
 *
 *  Set limits of background subimage to match those of probability array.
 *
 *
 ******************************************************************************/

int get_limits(fitsfile *infits, int extended, int aperture,
	int *ymin, int *ymax)
{

	char wgtsfile[FLEN_VALUE], ycentname[FLEN_KEYWORD];
	int hdu, nywt, status=0;
	float ycentv, wgt_cent;
	fitsfile *wgtsfits;

	cf_verbose(3, "Determining Y limits for aperture %d background array.",
		aperture);

	/* Read the measured centroid of the target spectrum. */
	sprintf(ycentname,"YCENT%1d", aperture);
	FITS_movabs_hdu(infits, 1, NULL, &status);
	FITS_read_key(infits, TFLOAT, ycentname, &ycentv, NULL, &status);

	/* Read the centroid and size of the probability array. */ 
	FITS_read_key(infits, TSTRING, "WGTS_CAL", wgtsfile, NULL, &status);
	FITS_open_file(&wgtsfits, cf_cal_file(wgtsfile), READONLY, &status);
	hdu = extended * 8 + aperture + 1;
	cf_verbose(3, "WGTS_CAL: extended = %d, aperture = %d, hdu = %d",
		extended, aperture, hdu);
	FITS_movabs_hdu(wgtsfits, hdu, NULL, &status);
	FITS_read_key(wgtsfits, TINT, "NAXIS2", &nywt, NULL, &status);
	FITS_read_key(wgtsfits, TFLOAT, "WGT_CENT", &wgt_cent, NULL, &status);
	FITS_close_file(wgtsfits, &status);

        cf_verbose(3, "Spectrum centroid=%6.2f, weights centroid=%5.2f, "
	   "weights ny=%4d", ycentv, wgt_cent, nywt) ;

	/* Note that we allow these values to exceed the limits of the 
	 * detector.  This is required if the bkgd arrays are to cover
	 * the same area as the 2-D probability and data arrays. */

	*ymin = cf_nint (ycentv - wgt_cent);
        *ymax = *ymin + nywt - 1;

	return status ; 
}


/*********************************************************************
 *
 *    SUM_GOOD_PIXELS
 *
 *    Subroutine multiplies image and mask arrays, then sums the result.
 *     image = image from which to extract the data
 *     mask = mask for data image (1 = good pixels)
 *     npix = size of image and mask arrays
 *
 ***********************************************************************/


int
sum_good_pixels(float *image, float *mask, long npix, double *image_sum)
{
    
    long    i;

    *image_sum = 0;

    for (i = 0; i < npix; i++) 
            *image_sum += image[i] * mask[i];

    return 0 ;
}


/******************************************************************************
 *
 *             MAIN PROGRAM BEGINS HERE
 *
 ******************************************************************************/


int cf_scale_bkgd(fitsfile *infits, long nevents, float *x_in, float *y_in,
    float *w_in, unsigned char *channel, unsigned char *timeflags,
    unsigned char *locflags, long ngood, long *good_index, long dtime,
    long ntime, int *binx, int *biny, int *nx_lif, int *ny_lif,
    int *ymin_lif, float **img_lif, int *nx_sic, int *ny_sic, int *ymin_sic,
    float **img_sic)
{

   char CF_PRGM_ID[] = "cf_scale_bkgd";
   char CF_VER_NUM[] = "1.25";

    fitsfile  *bkgdfits, *scrnfits, *bchrfits, *parmfits ;
    char     buffer[FLEN_CARD], bkgdfile[FLEN_CARD], scrnfile[FLEN_CARD];
    char     bchrfile[FLEN_CARD], keyname[FLEN_CARD], parmfile[FLEN_CARD] ;
    char     run_bkgd[FLEN_VALUE]; 
    char     comment[FLEN_CARD], datestr[FLEN_CARD];
    int      nx, nxb, ny, status, hdutype, frow, felem, niter;
    int      src_ap[2], nullval, anynull, bkgd_num;
    int      ymin, ymax, phalow, nbinx, nxsam, bkgdtype ;
    int      bin_size, extended, keyval, timeref;
    double   data_sum_day, data_sum_night, data_sum;
    double   model_sum_night, model_sum_day, model_sum ;
    float    *xvarg, *bkgarr, *xvars0, *xvarm0 ;
    float    *bkgimg, *bkgimgd, *bkgimgn, *xvars, *xvarm;
    float    floatnull, exptime, expnight, expday, xvars_ave ;
    float    intcrarray[9], intcr, scattot, rms, rms0, dcr ;
    float    intcnt, scat_int_ratio, zero, mf, numtot, xvars_tot, xvarm_tot ;
    float    *dataimg, *maskimg;
    short    *ylim ;
    long     fpixel, npix, i, j, k, nsam;
    long     nbgood, npts, ndx, ndx1, ndx2, bkgd_ndx;

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");

    /* Initialize CFITSIO variables */
    anynull = 0 ;
    hdutype = 0 ;
    frow = 1 ;
    felem = 1 ;
    floatnull = 0. ;
    fpixel = 1 ;
    status = 0 ;

    /* Initialize other variables */
    data_sum_day=0.;
    data_sum_night=0.;
    zero = 0.;

    /* Open parameter file and check whether user wants a background model. */
    FITS_movabs_hdu(infits, 1, &hdutype, &status) ;
    FITS_read_key(infits, TSTRING, "PARM_CAL", parmfile, NULL, &status) ;
    cf_verbose(3, "parameter file is %s ", parmfile) ;
    FITS_open_file(&parmfits,cf_parm_file(parmfile), READONLY, &status) ;
    FITS_read_key(parmfits, TSTRING, "RUN_BKGD", run_bkgd, NULL, &status) ;
    FITS_close_file(parmfits, &status) ;
    cf_verbose(3, "RUN_BKGD = %s", run_bkgd) ;

    /* If the RUN_BKGD = NO, then fill background array with zeros. */ 
    if (!strncmp(run_bkgd, "N", 1) || !strncmp(run_bkgd, "n", 1)) {
        cf_verbose(1, "RUN_BKGD = NO.  No background subtraction.") ;
	FITS_write_comment(infits, " ", &status);
	FITS_write_comment(infits, "RUN_BKGD = NO.  No background subtraction.",
	     &status);
	fits_get_system_time(datestr, &timeref, &status);
        sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
        FITS_write_comment(infits, comment, &status);
	FITS_write_comment(infits, " ", &status);

        bin_size = 1 ;
        nx=nxb=NXMAX ;
        ny=NYMAX ;
        npix=nx * ny ;
        bkgimg = (float *) cf_calloc(npix, sizeof(float)) ;
        goto fin ;
    }

/***********************************************************************
 *
 *      Set up the parameters for the analysis 
 *
 ***********************************************************************/

    /* Set day and night exposure times */
    expday = (float) dtime ;
    expnight = (float) ntime ;
    exptime = expday + expnight ;

    /* Get name of background characterization file */
    FITS_read_key(infits, TSTRING, "BCHR_CAL",bchrfile, NULL, &status) ;
    cf_verbose(3, "BCHR parameter file is %s ",bchrfile) ;

    /* Get name of background calibration file */
    FITS_read_key(infits, TSTRING, "BKGD_CAL", bkgdfile, buffer, &status);
    cf_verbose(3,"Background calibration file = %s ",bkgdfile ) ;

    /* Get name of screening file */
    FITS_read_key(infits, TSTRING, "SCRN_CAL",scrnfile, NULL, &status) ;
    cf_verbose(4, "screening file is %s ",scrnfile) ; 

    /* Read pha screening limit from the input file header */
    FITS_read_key(infits,TINT,"PHALOW",&phalow, NULL, &status) ; 
    if (phalow > 8) phalow=8 ;
    cf_verbose(2, "pha screening limit = %d", phalow) ;

    /* Determine the type of background model to make by checking the
       BKGDTYPE flag in the screening file. If BKGDTYPE > 0, do a full
       model. If BKGDTYPE < 0, calculate a model based on exposure time. */
    FITS_open_file(&scrnfits,cf_parm_file(scrnfile), READONLY, &status) ;
    FITS_read_key(scrnfits, TINT, "BKGDTYPE",&bkgdtype, NULL, &status) ;
    FITS_close_file(scrnfits, &status) ;

    /* Read the intrinsic count rate from the BCHR file */
    FITS_open_file(&bchrfits,cf_parm_file(bchrfile), READONLY, &status) ; 
    FITS_movabs_hdu(bchrfits, 3, &hdutype, &status) ;
    FITS_read_img(bchrfits, TFLOAT, fpixel, 9, &nullval, intcrarray,
      &anynull, &status) ;
    FITS_close_file(bchrfits, &status) ;
    intcr=intcrarray[phalow] ;
    cf_verbose(2,"Initial guess for the intrinsic count rate = "
	"%.1f x 10^-7 c/s/pixel ", intcr) ;

    /* Read the limits of the background sample regions from the header
	of the input file */
    FITS_read_key(infits, TINT, "BKGD_NUM", &bkgd_num, NULL, &status) ;
    cf_verbose(3, "Number of background sample regions = %d ", bkgd_num) ;

    ylim = (short *) cf_calloc(bkgd_num*2, sizeof(short) ) ;
    for (i=0; i<bkgd_num; i++) {
	sprintf(keyname,"BKG_MIN%ld",i) ;
	FITS_read_key(infits, TINT, keyname,&keyval, NULL, &status) ; 
	ylim[2*i] = keyval ;
	sprintf(keyname,"BKG_MAX%ld",i) ;
	FITS_read_key(infits, TINT, keyname,&keyval, NULL, &status) ; 
	ylim[2*i+1] = keyval ;
    }
    cf_verbose(3, "Limits of background regions") ;
    for (i=0; i<bkgd_num; i++) cf_verbose(3,
	"region %d = %d to %d ", i, ylim[2*i], ylim[2*i+1] ) ;

    /* If no background sample regions have been specified (as is usual for
       HIST data), then generate a model based only on exposure time */
    if (bkgd_num == 0) {
	bkgdtype = -1 ;
	cf_verbose(1, "No background sample regions specified "
	    "- probably HIST") ;
    }


/*************************************************************************
 *    
 *   Read in the background calibration information
 *
 *************************************************************************/
 
    /* Open background image file and go to the night image hdu */
    FITS_open_file(&bkgdfits, cf_cal_file(bkgdfile), READONLY, &status);
    FITS_movabs_hdu(bkgdfits, 2, &hdutype, &status);
    FITS_read_key(bkgdfits, TINT, "NAXIS1", &nx, buffer, &status);
    FITS_read_key(bkgdfits, TINT, "NAXIS2", &ny, buffer, &status);
    cf_verbose(3, "Reading background calibration file:",
	"naxis1=%d,  naxis2 = %d ",nx,ny) ;
       
    /* Background image files are binned in X, but not Y.
	Determine the degree of compression. */
    bin_size= NXMAX/nx ; 
    cf_verbose(3, "Background images are binned by %d pixels", bin_size) ;
            
    /* Allocate space to hold the background, data, and mask images */
    npix = nx * ny ;
    nxb = nx ;
    bkgimgn = (float *) cf_calloc(npix, sizeof(float)) ;
    bkgimgd = (float *) cf_calloc(npix, sizeof(float)) ;
    bkgimg  = (float *) cf_calloc(npix, sizeof(float)) ;
    dataimg = (float *) cf_calloc(npix, sizeof(float)) ;
    maskimg = (float *) cf_calloc(npix, sizeof(float)) ;
    
    /* Read the night and day scattered-light images */
    FITS_read_img(bkgdfits, TFLOAT, fpixel, npix, &nullval, bkgimgn,
	&anynull, &status) ;
    FITS_movabs_hdu(bkgdfits, 3, &hdutype, &status);
    FITS_read_img(bkgdfits, TFLOAT, fpixel, npix, &nullval, bkgimgd,
	&anynull, &status) ;
    FITS_close_file(bkgdfits, &status) ;

    /* Rescale the background images if they have been compressed. */
    if (bin_size > 1) for (i=0; i<npix; i++) {
	bkgimgn[i] *= bin_size ;
	bkgimgd[i] *= bin_size ;
    }

    /* Generate a mask of the same dimensions as the background image.
	It is 1 in the user-defined background regions, but 0 around
	airglow lines and near the edges of the detector. */
    make_bad_pixel_mask (infits, bkgd_num, ylim, npix, nxb, maskimg, &nsam);


/****************************************************************************
 *
 *      Read the photon list.
 *	Store in an array the same size as the the background image. 
 *	Populate only regions flagged as good in the background mask.
 *	Omit photons flagged as airglow.
 *	This scheme rejects non-airglow photons that drift into airglow
 *	regions due to mirror or spacecraft motion.  We thus under-estimate
 *	the background, but assume that the resulting error is small.
 *
 ****************************************************************************/

    /* If the background model is to be calculated only by exposure time,
       then skip the IDF file analysis */
    if (bkgdtype < 0) {
      cf_verbose(3, "Skipping the analysis of the IDF file") ;
    }
    else {	/* Analysis of the data begins here. */
      cf_verbose(3, "Beginning the analysis of the IDF file") ;

    nbgood=0;
    for (i=0; i<ngood; i++) {
	ndx=good_index[i] ;
	if (locflags[ndx] & LOCATION_AIR) continue;
	if ((j = cf_nint(y_in[ndx])) < 0 || j >= NYMAX) continue;
	bkgd_ndx = j*nxb + cf_nint(x_in[ndx])/bin_size;
	if (maskimg[bkgd_ndx] > 0) {
		dataimg[bkgd_ndx] += w_in[ndx];
      		if (timeflags[ndx] & TEMPORAL_DAY) data_sum_day += w_in[ndx] ;
       		else data_sum_night += w_in[ndx]; 
	  	nbgood++ ;
	}
    }

    cf_verbose(3, "Number of photons events in the background region = %ld",
	nbgood); 
    cf_verbose(3, "Sum of daytime weights  = %.1lf ", data_sum_day) ;
    cf_verbose(3, "Sum of nighttime weights  = %.1lf ", data_sum_night) ;

    /* If the background scattered light contribution is too large, the
	background may be contaminated, perhaps by a source in one of the
	other apertures.  This is common in crowded fields or with nebulocity.
	In this case, we calculate the background model based on the 
	exposure times, as done for the histogram data */

    /* Check the level of the scattered-light background */
    intcnt = intcr * nsam * bin_size * exptime * 1.e-7 ;
    scattot = data_sum_day + data_sum_night - intcnt ;
    scat_int_ratio = 0. ;
    if (intcnt > 0) scat_int_ratio = scattot / intcnt ;
    cf_verbose(3, "int cnt = %.1f, scat cnt = %.1f , scat/int = %f ",
	intcnt, scattot, scat_int_ratio) ;

    if (scat_int_ratio > 10) {
	cf_verbose(1, "The scattered light background is too large.");
	bkgdtype = -1;
    }
    else {
	cf_verbose(2,
	    "Scattered light appears reasonable. Proceeding with analysis.") ;
    }

    }	/* Analysis of the data ends here. */

    /* If the background model is to be calculated using only the exposure
	time, then do it now. */
    if (bkgdtype < 0) {
        cf_verbose(1, "Constructing background based only on exposure time.") ;
	make_model_exptime( npix, nxb, intcr, expnight, expday,
	    bkgimgn, bkgimgd, bkgimg) ;
        goto fin ; 
    }

/****************************************************************************
 *
 *      In this section, we attempt to determine the intrinsic count rate.
 *      The scattered light shows little structure in the Y dimension, so
 *	we project all of the arrays onto the X axis.  Only good pixels
 *	in the user-specified background regions are included in the sum.
 *
 ****************************************************************************/

    /* Make a scattered light reference image assuming an intrinsic count
       rate of zero. */
    make_model_exptime(npix, nxb, zero, expnight, expday,
	bkgimgn, bkgimgd, bkgimg) ;

    /* We bin the data by another factor of 8 in X      */
    /*   nbinx = number of detector pixels in one bin	*/
    /*   nxsam = number of samples in the array		*/
    nbinx=128 ;
    nxsam = NXMAX/nbinx ;
    cf_verbose(3, "X variation array information: "
                  "nbinx = %d, nxsam= %d ", nbinx, nxsam) ; 

    xvars  = (float *) cf_calloc(nxsam,sizeof(float)) ;
    xvarm  = (float *) cf_calloc(nxsam,sizeof(float)) ;
    xvars0 = (float *) cf_calloc(nxsam,sizeof(float)) ;
    xvarm0 = (float *) cf_calloc(nxsam,sizeof(float)) ;
    xvarg  = (float *) cf_calloc(nxsam,sizeof(float)) ;

    /* Project each image onto the X axis */
    get_xvar_img(bkgimg,  maskimg, nxb, nxsam, xvarm0) ;
    get_xvar_img(dataimg, maskimg, nxb, nxsam, xvars0) ;
    get_xvar_img(maskimg, maskimg, nxb, nxsam, xvarg) ;

    /* Calculate the average counts per "good" bin in the data array. */
    xvars_ave = 0;
    for (i = j = 0; i < nxsam; i++)
	if (xvarg[i] > 0) {
	    xvars_ave += xvars0[i];
	    j++;
	}
    xvars_ave /= j;
    cf_verbose (3, "xvars_ave = %.1f (should be > 100)", xvars_ave);

    /* If the average number of counts in the xvars0 array (which holds 
    the observed intensity variation along the x direction) is small,
    then we can't determine an intrinsic count rate from the data and
    we use the tabulated rate.  If the count is large enough, then we
    can do a more detailed fit of the model to the observations. */

    if (xvars_ave > 100) {	/* Begin fit to intrinsic background */

    cf_verbose(1, "Calculating a full background model.") ;
    cf_verbose(3, "\n\tInterpolating to derive intrinsic count rate\n") ;

    /* Sum the counts in the data and model arrays. */
    sum_good_pixels(dataimg, maskimg, npix, &data_sum);
    sum_good_pixels(bkgimg, maskimg, npix, &model_sum);

    cf_verbose(3, "Total counts in initial scattered-light model = %7.0lf",
	model_sum);
    cf_verbose(3, "Intrinsic counts summed over good regions =     %7.0f",
	(intcnt = intcr * nsam * bin_size * exptime * 1.e-7));
    cf_verbose(3, "Sum of above                                  = %7.0lf",
	model_sum + intcnt);
    cf_verbose(3, "Total counts in good regions of data array =    %7.0lf",
	data_sum);

    /* Sum the counts in the data and model x variation arrays */
    xvarm_tot = 0.;
    xvars_tot = 0.;
    for (i = j = 0; i < nxsam; i++) if (xvarg[i] > 0) {
	xvarm_tot += xvarm0[i];
	xvars_tot += xvars0[i];
    }
    cf_verbose(3, "Total counts in scattered-light xvar array =    %7.0f",
	xvarm_tot);
    cf_verbose(3, "Total counts in data x variation array =        %7.0f",
	xvars_tot);

    /* Convert the x variation arrays from total counts to units of
	1.e-7 c/s/pix, the units of the intrinsic count rate. */
    mf = bin_size * exptime * 1.e-7 ;
    for (i=0 ; i<nxsam ; i++) if (xvarg[i] > 0) {
	xvars0[i] /= xvarg[i] * mf ;
	xvarm0[i] /= xvarg[i] * mf ;
    }

    /* Sum the counts in the rescaled xvar arrays */
    xvarm_tot = 0.;
    xvars_tot = 0.;
    for (i = j = 0; i < nxsam; i++) if (xvarg[i] > 0) {
	xvars_tot += xvars0[i];
	xvarm_tot += xvarm0[i];
	j++;
    }
    cf_verbose(3, "Mean xvar count rates in units of 10^-7 c/s/pixel:");
    cf_verbose(3, "     Initial model = %.1f", xvarm_tot/j);
    cf_verbose(3, "  Intrinsic counts = %.1f", intcr);
    cf_verbose(3, "        Data array = %.1f", xvars_tot/j);

    /* Determine the value of the intrinsic background by comparing the X
	variations of the data (after removing a guess for the intrinsic
	background) with the model scattered light variations, making sure
	that the final average count rate is the same for both.  Iterate
	until the intrinsic background changes by less than 0.1 units. */

    dcr = 1. ;
    niter = 0 ;
    rms0 = 1.e10 ;
    while ( fabs(dcr) > 0.1) { 

        /* Subtract intrinsic counts from the data array and sum it. */
	xvars_tot = 0. ;
	for (i=0; i<nxsam; i++) if (xvarg[i] > 0) {
	    xvars[i] = xvars0[i] - intcr ;
	    xvars_tot += xvars[i] ;
	}

    	/* Scale the model array to match the data.  Sum it. */
	numtot = 0. ;
	for (i=0 ; i<nxsam ; i++) if (xvarg[i] > 0) {
	    xvarm[i] = xvarm0[i] * (xvars_tot / xvarm_tot) ;
	    numtot += xvarm[i] ;
	}
        
	/* Sum the absolute difference between the data and model arrays. */
	rms = 0. ;
	for (i=0 ; i<nxsam ; i++) if (xvarg[i] > 0)  
	    rms += (fabs(xvars[i]-xvarm[i])) ; 
	cf_verbose(3, "intr cr= %5.2f, data-model = %6.2f, "
	    "dat cnts = %5.0f, mod cnts=%5.0f", intcr,rms,xvars_tot,numtot) ;

	niter++ ;
	if (niter > 20) break ;

	/* if the current rsm value is greater than the preceeding one,
	   reverse direction and reduce the step size by 50% */
	if ( rms > rms0 )  dcr = -0.5*dcr ;
	    intcr += dcr ;
 
	/* initialize the next iteration */
	rms0 = rms ;
    }
    }		/* End of fit to intrinsic background */

    else { 
    cf_verbose(1, "Background count too small for a detailed fit.");
    cf_verbose(1, "Assuming an intrinsic count rate of %.1f cts/s", intcr ) ;
    }

    /* Sum the good pixels in the unscaled day and night background images. */
    sum_good_pixels(bkgimgn, maskimg, npix, &model_sum_night);
    sum_good_pixels(bkgimgd, maskimg, npix, &model_sum_day);

    /* Make the final background model */
    make_model(npix, nsam, nxb, intcr, expnight, expday,
	data_sum_night, data_sum_day, model_sum_night, model_sum_day, 
	bkgimgn, bkgimgd, bkgimg, maskimg) ;

    /* Sum the counts in the final background image. */
    sum_good_pixels(bkgimg, maskimg, npix, &model_sum);

    cf_verbose(3, "Total counts in final background model = %.0lf", model_sum);
    cf_verbose(3, "                       measured counts = %.0lf",
	data_sum_day+data_sum_night);

    free(bkgimgn) ;
    free(bkgimgd) ;
    free(dataimg) ;
    free(maskimg) ;
    free(xvars) ;
    free(xvars0) ;
    free(xvarm) ;
    free(xvarm0) ;
    free(xvarg) ;

    fin:

  /**************************************************************************
   *
   *     Extract the subarrays covering the active LiF and SiC apertures 
   *
   ***************************************************************************/

    extended = cf_source_aper(infits, src_ap) ;
    cf_verbose(3, "\n\tSelecting subarrays\n") ;
    cf_verbose(3, "source apertures: lif=%d,  sic=%d ",src_ap[0], src_ap[1]) ; 

    /* Populate some output values */
    *binx=bin_size ;
    *biny=1 ;

    for (k=0; k<2 ; k++) {	/* Loop over LiF and SiC channels */

    /* Read the y limits for the output array */
    get_limits(infits, extended, src_ap[k], &ymin, &ymax) ;
  
    /* Set up array to contain the image */ 
    npts = nxb * (ymax-ymin+1) ;
    bkgarr = (float *) cf_calloc(npts, sizeof(float)) ;
				   
    /* Fill the array */
    ndx1 = -1 ;
    model_sum = 0 ;
    for (i=ymin; i<=ymax ; i++) {
	if (i < 0 || i >= ny)
	    ndx1 += nxb ;
	else {
	    for (j=0; j<nxb; j++) {
	        ndx2 = i*nxb + j ;
	        if (ndx2 >= npix) break ;
	        ndx1++;
	        if (ndx1 >= npts) break ;
	        bkgarr[ndx1]=bkgimg[ndx2] ;
		model_sum += bkgarr[ndx1] ;
	    }
	}
    }

    /* Populate the output arrays and values */
    if (k == 0) {
	*nx_lif = nxb ;
	*ny_lif = (ymax - ymin + 1) ;
	*ymin_lif = ymin ;
	*img_lif = bkgarr ;
    }
    else {
	*nx_sic = nxb ;
	*ny_sic = (ymax - ymin + 1) ;
	*ymin_sic = ymin ;
	*img_sic = bkgarr ;
    }

    cf_verbose(2, "for aperture %d, nx=%d, ny=%d, ymin=%d, model_sum=%.0f",
	src_ap[k], nxb, (ymax-ymin+1), ymin, model_sum) ; 

    }	/* End loop over LiF and SiC channels */

    /* release memory for the background image */
    free(bkgimg) ;

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
    return 0;
}
