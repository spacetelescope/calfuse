/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_screen_bad_pixels (fitsfile *header, long nevents,
 *			      float *x, float *y, unsigned char *location)
 *
 * Description: Flag photon events in bad-pixel regions.
 *
 * Arguments:   fitsfile  *header	Input FITS file pointer
 *              long      nevents	Number of photon events
 *              float     *x		X coordinate array
 *              float     *y		Y coordinate array
 *              unsigned char *location Photon location flag array
 *
 * Calls:	cf_get_potholes
 *
 * Returns:	0 on success
 *
 * History:	11/22/05   1.1   wvd    Initial coding
 *
 ****************************************************************************/

#include <stdlib.h>
#include "calfuse.h"

/***********************************************************************
 *
 * Procedure to read pothole locations from the QUAL_CAL file
 *
 ***********************************************************************/

int
cf_get_potholes(fitsfile *header, float **xbad, float **ybad,
        float **rxbad, float **rybad)
{
    char qualfile[FLEN_VALUE];
    int nbad, status=0;
    fitsfile *qualfits;

    cf_verbose(3, "Entering cf_get_potholes.");

    /*
     *  Read limits of bad-pixel regions from the first extension 
     *  of the QUAL_CAL file.
     */
    FITS_movabs_hdu(header, 1, NULL, &status);
    FITS_read_key(header, TSTRING, "QUAL_CAL", qualfile, NULL, &status);
    cf_verbose(3, "QUAL_CAL = %s", qualfile);

    FITS_open_file(&qualfits, cf_cal_file(qualfile), READONLY, &status);
    FITS_movabs_hdu(qualfits, 2, NULL, &status);
    nbad = cf_read_col(qualfits, TFLOAT, "X",  (void **) xbad);
    nbad = cf_read_col(qualfits, TFLOAT, "Y",  (void **) ybad);
    nbad = cf_read_col(qualfits, TFLOAT, "RX", (void **) rxbad);
    nbad = cf_read_col(qualfits, TFLOAT, "RY", (void **) rybad);
    FITS_close_file(qualfits, &status);

    cf_verbose(3, "Exiting cf_get_potholes.");
    return nbad;
}


int
cf_screen_bad_pixels (fitsfile *header, long nevents, float *x, float *y,
	unsigned char *location)
{
    char CF_PRGM_ID[] = "cf_screen_bad_pixels";
    char CF_VER_NUM[] = "1.1";

    int   errflg=0, nbad;
    long  counter=0, i, j;
    float *xbad, *ybad, *rxbad, *rybad, *rx2, *ry2;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /*
     *  Read limits of bad-pixel regions from QUAL_CAL file.
     */
    nbad = cf_get_potholes(header, &xbad, &ybad, &rxbad, &rybad);
    cf_verbose(2, "For first pothole: x = %.0f, y = %.0f, rx = %.0f, ry = %.0f",
	xbad[0], ybad[0], rxbad[0], rybad[0]);

    /*
     *  Square the semi-major axes.
     */
    for (j = 0; j < nbad; j++) {
	rxbad[j] *= rxbad[j];
	rybad[j] *= rybad[j];
    }
    rx2 = rxbad;
    ry2 = rybad;
    cf_verbose(2, "For first pothole: rx2 = %.0f, ry2 = %.0f", rx2[0], ry2[0]);

    /*
     *  Flag photons in the bad-pixel regions.
     */
    for (i = 0; i < nevents; i++) {
	for (j = 0; j < nbad; j++) {
            if ((x[i] - xbad[j]) * (x[i] - xbad[j]) / rx2[j] +
		(y[i] - ybad[j]) * (y[i] - ybad[j]) / ry2[j] < 1.0) {
		location[i] = location[i] | LOCATION_BADPX;
		counter++;
		break;
	    }
        }
    }
    cf_verbose(2, "Flagged %d pixels as bad", counter);

    free(xbad);
    free(ybad);
    free(rxbad);
    free(rybad);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
