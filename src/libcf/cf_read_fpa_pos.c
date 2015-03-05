/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_read_fpa_pos (fitsfile *infits, float *fpaxsic, float *fpaxlif)
 *
 * Description: Returns SiC and LiF FPA positions in microns.  Can read both
 *              old and new FPA position keywords.
 *
 * Arguments:   fitsfile    *infile       Input data file name
 *              float       *fpaxsic      FPA SiC position (returned)
 *              float       *fpaxlif      FPA LiF position (returned)
 *
 * Returns:     int	     status	  CFITSIO status value
 *
 * History:     08/21/01   v1.1   wvd     Subroutine taken from cf_wfits.c
 * 		01/14/03   v1.3   wvd     Change name to cf_read_fpa_pos.
 *              12/18/03   v1.4   bjg     Change calfusettag.h to calfuse.h
 *
 ****************************************************************************/

#include <string.h>
#include "calfuse.h"

/* for all FPA's an increase in "data number" moves the FPA towards +Xipcs */
#define FPALIF1DN2X0 -69.949	/* DN to microns conversion, zero-point */
#define FPALIF1DN2X1 .14144	/* DN to microns conversion, linear term  */
#define FPALIF2DN2X0 -73.189	/* DN to microns conversion, zero-point */
#define FPALIF2DN2X1 .14487	/* DN to microns conversion, linear term  */
#define FPASIC1DN2X0 -132.507	/* DN to microns conversion, zero-point */
#define FPASIC1DN2X1 .14156	/* DN to microns conversion, linear term  */
#define FPASIC2DN2X0 -89.116	/* DN to microns conversion, zero-point */
#define FPASIC2DN2X1 .14206	/* DN to microns conversion, linear term  */

int cf_read_fpa_pos (fitsfile *infits, float *fpaxsic, float *fpaxlif)
{
    char  detector[FLEN_VALUE]; /* detector keyword from FITS header */
    int	status=0;
    float tmpxsic, tmpxlif;

    FITS_read_key(infits, TSTRING, "DETECTOR", detector, NULL, &status);

    /* Get FPA positions */
    /* Use ffgky instead of FITS_read_key so we can read either the old
     * or new FPA position keywords */
    ffgky (infits, TFLOAT, "FPASXPOS", fpaxsic, NULL, &status);
    if (status != 0) {
	/* New keywords are not present; must be the old ones.
	 * Determine which detector provided current spectrum
	 * and convert raw FPA position to microns.
	 */
	status = 0;
	if (strncmp(detector,"1", 1)==0) {
	    FITS_read_key (infits, TFLOAT, "FP1SXPOS", &tmpxsic, NULL, 
			   &status);
	    FITS_read_key (infits, TFLOAT, "FP1LXPOS", &tmpxlif, NULL, 
			   &status);
	    *fpaxsic = FPASIC1DN2X0 + FPASIC1DN2X1 * tmpxsic;
	    *fpaxlif = FPALIF1DN2X0 + FPALIF1DN2X1 * tmpxlif;
	}
	else if (strncmp(detector,"2", 1)==0) {
	    FITS_read_key (infits, TFLOAT, "FP2SXPOS", &tmpxsic, NULL, 
			   &status);
	    FITS_read_key (infits, TFLOAT, "FP2LXPOS", &tmpxlif, NULL, 
			   &status);
	    *fpaxsic = FPASIC2DN2X0 + FPASIC2DN2X1 * tmpxsic;
	    *fpaxlif = FPALIF2DN2X0 + FPALIF2DN2X1 * tmpxlif;
	}
	else {
	    cf_if_error("Could not parse DETECTOR keyword, cannot "
			"convert raw FPA positions.");
	}
    }
    else {
	/* New SiC FPA keyword present. Assume new LiF FPA keyword
	 * is also present */
	FITS_read_key(infits, TFLOAT, "FPALXPOS", fpaxlif, NULL, &status);
    }
    return status;
}
