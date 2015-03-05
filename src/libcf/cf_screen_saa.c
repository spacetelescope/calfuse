/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_screen_saa(fitsfile *infits, long nseconds,
 *                            unsigned char *timeline_status, 
 *			      float *timeline_long, float *timeline_lat);
 *
 * Description: Set the screening bit due to SAA boundary.
 *
 * Arguments:   fitsfile  *infits       Input FITS file pointer
 *              long      nseconds      Number of points in the timeline table
 *              unsigned char  *timeline_status The status flag array
 *              float     *timeline_long The longitude array
 *              float     *timeline_lat  The latitude array
 *
 * Calls:	saa
 *
 * Returns:	0 on success
 *
 * History:	10/29/02   1.1   jch    Initial coding
 * 		12/20/02   1.3   jch    Include reading SAA table here
 * 		12/20/02   1.4   wvd    Delete call to calfuse.h 
 *              05/20/03   1.5   rdr    Add call to cf_proc_check
 *              08/28/03   1.6   bjg    New definition of structure saareg
 * 		09/02/03   1.7   wvd    Change calfusettag.h to calfuse.h
 * 		09/10/03   1.8   wvd    Tidy up code.
 *              02/17/05   1.9   wvd    Place parentheses around assignment
 *                                      used as truth value.
 *
 ****************************************************************************/

#include <stdlib.h>

#include "calfuse.h"

int
cf_screen_saa(fitsfile *infits, long nseconds, unsigned char *timeline_status,
	      float *timeline_long, float *timeline_lat)
{
    char CF_PRGM_ID[] = "cf_screen_saa";
    char CF_VER_NUM[] = "1.9";

    char  saa_cal[FLEN_VALUE];
    int   errflg=0, status=0, hdutype;
    long  i;
    saareg saa_bounds;
    fitsfile *saafits;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    /*
     *  Read SAA boundaries from the specified calibration file.
     */
    FITS_read_key(infits, TSTRING, "SAAC_CAL", saa_cal, NULL, &status);
    FITS_open_file(&saafits, cf_cal_file(saa_cal), READONLY, &status);
    FITS_movabs_hdu(saafits, 2, &hdutype, &status);

    saa_bounds.n_points =
	   cf_read_col(saafits, TFLOAT, "LATITUDE",  (void **) &(saa_bounds.lat));
    (void) cf_read_col(saafits, TFLOAT, "LONGITUDE", (void **) &(saa_bounds.lon));
    
    fits_close_file(saafits, &status);

    /*
     *  Set the SAA boundary bit if in SAA
     */
    for (i = 0; i < nseconds; i++) {
	if (saa(&saa_bounds, (double)timeline_long[i], (double)timeline_lat[i])) {
	    timeline_status[i] |= TEMPORAL_SAA; 
	}
    }
    free(saa_bounds.lon);
    free(saa_bounds.lat);

    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
