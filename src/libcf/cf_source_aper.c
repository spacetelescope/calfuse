/**************************************************************************
 *           Johns Hopkins University
 *           Center for Astrophysical Sciences
 *           FUSE
 *************************************************************************
 *
 * Synopsis:    cf_source_aper(infits, ap)
 *
 * Description: Procedure to determine which aperture contains the source
 *
 * Arguments :  fitsfile     infits     Pointer to Integmediate Data File
 *              int          *ap        Array containing the two active
 *                                      apertures
 *
 * Returns      int          srctype    Type of source in aperture
 *
 * HISTORY      12/12/02  v1.1   rdr    First delivery
 *              12/13/02  v1.3   wvd    Changed file name
 *              02/24/03  v1.4   peb    Changed include file to calfitsio.h
 *		04/17/03  v1.5   wvd	Change to cf_source_aper throughout.
 *              07/11/03  v1.6   rdr    Change so that src_type returns
 *                                      1 rather than 3 for extended source
 *		08/21/03  v1.7   wvd	If APERTURE=RFPT, treat as LWRS.
 *              04/09/04  v1.8   bjg    Include string.h
 *                                      Remove unused variables
 *
 ************************************************************************/

#include "calfuse.h"
#include <string.h>

int
cf_source_aper(fitsfile *infits, int *active_ap )
{
     /* char CF_PRGM_ID[] = "cf_source_aper";         */
     /* char CF_VER_NUM[] = "1.8";                    */

    char aper[FLEN_CARD], src[FLEN_CARD];
    int status=0, hdutype, src_type;

    /* Determine the active aperture */
    FITS_movabs_hdu(infits, 1, &hdutype, &status); 
    FITS_read_key(infits, TSTRING, "APERTURE", aper, NULL, &status);
 
    /* Specify the channels of the active aperture */
    if (!strncmp(aper,"HIRS",4)) {
	active_ap[0] = 1;
	active_ap[1] = 5; }
    else if (!strncmp(aper,"MDRS",4)) {
	active_ap[0] = 2;
	active_ap[1] = 6; }
    else if (!strncmp(aper,"LWRS",4) || !strncmp(aper,"RFPT",4)) {
	active_ap[0] = 3;
	active_ap[1] = 7;
    }
    else
	cf_if_error("APERTURE keyword corrupted.");
    cf_verbose(3, "APERTURE = %s, channels are %d and %d",
		aper, active_ap[0], active_ap[1]);

    FITS_read_key(infits, TSTRING,"SRC_TYPE", src, NULL, &status);
    src_type = 0;
    if (!strncmp(src,"E",1) || !strncmp(src,"e",1) ) src_type = 1;

    cf_verbose(3, "SRC_TYPE = %s, source type = %d", src, src_type);
    return src_type;
}
