/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_pha_x_distort(fitsfile *header, long nevents, 
 *			unsigned char *pha, float *xfarf,
 *			unsigned char *locflags)
 *
 * Description: Corrects the X position of low-pulse-height events.
 *
 *		Note: To conform with the CalFUSE standard, the PHAX
 *		correction is now ADDED to the input photon X coordinate.
 *		This change from previous versions of the code requires
 *		a new set of calibration files.  The program now checks the
 *		version number of the PHAX_CAL file and returns an error if
 *		it is less than 4.
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              long      nevents       The number of events
 *              unsigned  char *pha     An array of PHA values
 *              float     *xfarf        An array of event X positions
 *              unsigned char *locflags Location flags of each event
 *
 * Calls:
 *
 * Return:      0  on success
 *
 * History:     11/11/02   1.1    peb   Begin work
 *              11/12/02   1.3    peb   Added check to move only events in
 *                                      active region.
 *		03/11/03   1.4    wvd   Changed pha and locflags to unsigned
 *					char
 *              05/20/03   1.5    rdr   Added call to cf_proc_check
 *		07/29/03   1.6    wvd	If cf_proc_check fails, return errflg.
 *		10/14/03   1.7    wvd	Add, rather than subtract, PHAX
 *					correction.  Check cal file version.
 *		10/17/03   1.8    wvd	Don't crash if CALFVERS keyword 
 *					is missing.
 *              04/07/07   1.9   wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_pha_x_distort(fitsfile *header, long nevents, unsigned char *pha,
		float *xfarf, unsigned char *locflags)
{
    char CF_PRGM_ID[] = "cf_pha_x_distort";
    char CF_VER_NUM[] = "1.9";

    char  phaxfile[FLEN_VALUE]={'\0'};
    int   anynull=0, calfvers, errflg=0, status=0, xlen, ylen;
    long  j;
    float *phax, flt=0.;
    fitsfile *phaxfits=NULL;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /*
     *  Read PHA correction (Walk) file
     */
    FITS_read_key(header, TSTRING, "PHAX_CAL", phaxfile, NULL, &status);
    FITS_open_file(&phaxfits, cf_cal_file(phaxfile), READONLY, &status);
    FITS_read_key(phaxfits, TINT, "NAXIS1", &xlen, NULL, &status);
    FITS_read_key(phaxfits, TINT, "NAXIS2", &ylen, NULL, &status);

    /* If CALFVERS keyword missing or less than 4, abort pipeline. */
    fits_read_key(phaxfits, TINT, "CALFVERS", &calfvers, NULL, &status);
    if (status || calfvers < 4) {
	status = 0;
	FITS_close_file(phaxfits, &status);
	cf_if_error("Old versions of PHAX_CAL have wrong sign.  Aborting.");
    }

    /* Otherwise, read walk correction as an image. */
    phax = (float *) cf_malloc(sizeof(float)*xlen*ylen);
    FITS_read_img(phaxfits, TFLOAT, 1L, xlen*ylen, &flt, phax, &anynull, &status);
    FITS_close_file(phaxfits, &status);

    if (xlen != NXMAX)
	cf_if_warning("NAXIS1 of %s != 16384", phaxfile);
    /*
     *  Apply PHA X correction to each event.
     */
    for(j=0; j<nevents; j++) {
	/*
	 *  Move only events in active region.
	 */
	if (!(locflags[j] & LOCATION_SHLD)) {
	    xfarf[j] += phax[(int)(pha[j]*xlen + xfarf[j])];
	}
    }
    free(phax);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
