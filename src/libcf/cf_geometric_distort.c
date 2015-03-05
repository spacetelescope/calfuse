/****************************************************************************
 *
 * Synopsis:    cf_geometric_distort(fitsfile *header, long nevents,
 *                                   float *xfarf, float *yfarf,
 *                                   unsigned char *locflags)
 *
 * Description: Corrects the X and Y positions of the event for geometric
 *              distortion.
 *
 * Arguments:   fitsfile  *header       Pointer to the location of the FITS
 *                                      header of the Intermediate Data File
 *              long      nevents       number of photons in the file
 *              float     *xfarf        X positions for each photon event
 *              float     *yfarf        Y positions for each photon event
 *              unsigned char *locflags Location flags for each event
 *
 * Calls :      None
 *
 * Return:      0  on success
 *
 * History:     08/09/02   1.1   rdr    Begin work 
 *              11/12/02   1.2   peb    Added check to move only events in
 *                                      active region.  Added locflags to
 *                                      function cal.
 *		03/11/03   1.3   wvd    Changed locflags to unsigned char.
 *              05/20/03   1.4   rdr    Added call to cf_proc_check
 *              07/29/03   1.5   wvd    If cf_proc_check fails, return errflg.
 *              03/26/04   1.6   rdr    Modify specbiny
 *              04/05/04   1.7   wvd	Modify specbiny only if specbiny > 1.
 *              04/07/07   1.8   wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_geometric_distort(fitsfile *header, long nevents, float *xfarf,
		     float *yfarf, unsigned char *locflags)
{
    char CF_PRGM_ID[] = "cf_geometric_distort";
    char CF_VER_NUM[] = "1.8";

    char geomfile[FLEN_VALUE]={'\0'};
    short *geomx_img, *geomy_img;
    int   errflg=0, status=0, hdutype, anynull, nullval;
    int   specbiny ;
    long  fpixel=1, j;
    long  xdim_gx, ydim_gx, xdim_gy, ydim_gy, npix_gx, npix_gy;
    float scale_gx=1., zero_gx=0., scale_gy=1., zero_gy=0.;
    float yexp ;
    fitsfile *geomfits=NULL;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /*
     *  Get the name of the calibration file and open it
     */
    FITS_read_key(header, TINT, "SPECBINY", &specbiny, NULL, &status);
    cf_verbose(5,"Initial value of specbiny=%d",specbiny) ;
    FITS_read_key(header, TSTRING, "GEOM_CAL", geomfile, NULL, &status);
    cf_verbose(4, "Geometric correction file = %s", geomfile);
    FITS_open_file(&geomfits, cf_cal_file(geomfile), READONLY, &status);

    /* Read the y expansion factor from the calibration header */
    FITS_movabs_hdu(geomfits, 1, &hdutype, &status);
    FITS_read_key(geomfits, TFLOAT,  "YEXPAND", &yexp,  NULL, &status);
    cf_verbose(3,"Y scale expansion factor = %f",yexp) ;

    /* If SPECBINY > 1, scale by Y expansion factor */
    if (specbiny > 1) {
       specbiny = (int) (0.5 + (float) specbiny * yexp) ;
       FITS_update_key(header, TINT, "SPECBINY", &specbiny, NULL, &status) ;
       cf_verbose(5,"Corrected value of specbiny = %d",specbiny) ;
    }

    /*
     *  Read in the X distortion calibration image
     */
    FITS_movabs_hdu(geomfits, 2, &hdutype, &status);
    FITS_read_key(geomfits, TLONG,  "NAXIS1", &xdim_gx,  NULL, &status);
    FITS_read_key(geomfits, TLONG,  "NAXIS2", &ydim_gx,  NULL, &status);
    npix_gx = xdim_gx * ydim_gx;

    FITS_read_key(geomfits, TFLOAT, "BSCALE", &scale_gx, NULL, &status);
    FITS_read_key(geomfits, TFLOAT, "BZERO",  &zero_gx,  NULL, &status);
    geomx_img = (short *) cf_malloc(sizeof(short) * npix_gx);
    fits_set_bscale(geomfits, 1., 0., &status);
    FITS_read_img(geomfits, TSHORT, fpixel, npix_gx, &nullval, geomx_img,
		  &anynull, &status);
    /*
     *  Read the Y distortion calibration image
     */
    FITS_movabs_hdu(geomfits, 3, &hdutype, &status);
    FITS_read_key(geomfits, TLONG,  "NAXIS1", &xdim_gy,  NULL, &status);
    FITS_read_key(geomfits, TLONG,  "NAXIS2", &ydim_gy,  NULL, &status);
    npix_gy = xdim_gy * ydim_gy;

    FITS_read_key(geomfits, TFLOAT, "BSCALE", &scale_gy, NULL, &status);
    FITS_read_key(geomfits, TFLOAT, "BZERO",  &zero_gy,  NULL, &status);
    geomy_img = (short *) cf_malloc(sizeof(short) * npix_gy);
    fits_set_bscale(geomfits, 1., 0., &status);
    FITS_read_img(geomfits, TSHORT, fpixel, npix_gy, &nullval, geomy_img,
		  &anynull, &status);
    FITS_close_file(geomfits, &status);
    /*
     *  Make sure that the X and Y distortion images are the same size
     */
    if (xdim_gx != xdim_gy || ydim_gx != ydim_gy) 
	cf_if_warning("X and Y distortion images not the same size.");
    /*
     *  Use the data in the geometric distortion calibration images to
     *  correct XFARF and YFARF
     */
    for (j=0; j<nevents; j++) {
	/*
	 *  Move only events in active region.
	 */
	if (!(locflags[j] & LOCATION_SHLD)) {
	    int xndx = (int) xfarf[j], yndx = (int) yfarf[j], ndx;
	    /*
	     *  Check that the y indices are within bounds
	     */
	    if (yndx < 0)
		yndx = 0;
	    else if (yndx > ydim_gy-1)
		yndx = ydim_gy-1;

	    ndx = yndx*xdim_gx + xndx;
	    xfarf[j] += geomx_img[ndx]*scale_gx + zero_gx;
	    yfarf[j] += geomy_img[ndx]*scale_gy + zero_gy;
	}
    }
    free(geomx_img);
    free(geomy_img);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");

    return 0;
}
