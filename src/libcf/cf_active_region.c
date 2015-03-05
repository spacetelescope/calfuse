/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_active_region(fitsfile *header, long nevents, float *xfarf,
 *			float *yfarf, unsigned char *locflags)
 *
 * Description: Flags events outside of detector active area.
 *
 *		Note: Both cf_ttag_init and cf_hist_init set LOCATION_SHLD
 *		when initializing IDF, but they look only at XRAW and use
 *		limits that are less restrictive than those used here.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nevents       The number of events
 *              float     *xfarf        An array of event X positions
 *              float     *yfarf        An array of event Y positions
 *              unsigned char *locflags Location flags of each event
 *
 * Calls:
 *
 * Return:      0  on success
 *
 * History:     11/12/02   1.1    peb   Begin and finish work
 *		03/11/03   1.2    wvd   Changed locflags to unsigned char
 *              05/20/03   1.3    rdr   Added call to cf_proc_check 
 *		07/29/03   1.4    wvd   If cf_proc_check fails, return errflg.
 *              04/07/07   1.5    wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_active_region(fitsfile *header, long nevents, float *xfarf, float *yfarf,
		 unsigned char *locflags)
{
    char CF_PRGM_ID[] = "cf_active_region";
    char CF_VER_NUM[] = "1.5";

    char  elecfile[FLEN_VALUE]={'\0'};
    char  keyword[FLEN_KEYWORD]={'\0'}, detector[FLEN_VALUE]={'\0'};
    int   errflg=0, status=0, active_l, active_r, active_b, active_t;
    long  j;
    fitsfile *elecfits=NULL;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    FITS_read_key(header, TSTRING, "DETECTOR", detector, NULL, &status);
    /*
     *  Read keywords from ELEC_CAL file
     */
    FITS_read_key(header, TSTRING, "ELEC_CAL", elecfile, NULL, &status);
    FITS_open_file(&elecfits, cf_cal_file(elecfile), READONLY, &status);
    sprintf(keyword, "ACTVL_%s", detector);
    FITS_read_key(elecfits, TINT, keyword, &active_l, NULL, &status);
    sprintf(keyword, "ACTVR_%s", detector);
    FITS_read_key(elecfits, TINT, keyword, &active_r, NULL, &status);
    sprintf(keyword, "ACTVB_%s", detector);
    FITS_read_key(elecfits, TINT, keyword, &active_b, NULL, &status);
    sprintf(keyword, "ACTVT_%s", detector);
    FITS_read_key(elecfits, TINT, keyword, &active_t, NULL, &status);
    FITS_close_file(elecfits, &status);
    cf_verbose(3, "Active area limits: X from %d to %d, Y from %d to %d",
	active_l, active_r, active_b, active_t);
    /*
     *  Apply active region flags to each event.
     */
    for(j=0; j<nevents; j++) {
	/*
	 *  Flag events that fall outside of detector active region.
	 */
	if (xfarf[j] < active_l || xfarf[j] > active_r ||
	    yfarf[j] < active_b || yfarf[j] > active_t) {
	    locflags[j] |= LOCATION_SHLD;
	}
    }

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
