/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_screen_airglow (fitsfile *header, long nevents,
 *			      float *x, float *y, unsigned char *location)
 *
 * Description: Flag photon events in airglow regions.
 *
 * Arguments:   fitsfile  *header	Input FITS file pointer
 *              long      nevents	Number of photon events
 *              float     *x		X coordinate array
 *              float     *y		Y coordinate array
 *              unsigned char *location Photon location flag array
 *
 * Calls:	cf_get_geocorona
 *
 * Returns:	0 on success
 *
 * History:	03/10/05   1.1   wvd    Initial coding
 *
 ****************************************************************************/

#include <stdlib.h>
#include "calfuse.h"

int
cf_screen_airglow (fitsfile *header, long nevents, float *x, float *y,
	unsigned char *location)
{
    char CF_PRGM_ID[] = "cf_screen_airglow";
    char CF_VER_NUM[] = "1.1";

    int   errflg=0, ngeo;
    long  i, j;
    short *gxmin, *gxmax, *gymin, *gymax;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /*
     *  Read limits of airglow regions from AIRG_CAL file.
     */
    ngeo = cf_get_geocorona(header, &gxmin, &gxmax, &gymin, &gymax);

    /*
     *  Set the airglow bit of photons in the airglow regions.
     */
    for (i = 0; i < nevents; i++) {
	for (j=0; j<ngeo; j++) {
            if (x[i] > gxmin[j] && x[i] < gxmax[j] &&
		y[i] > gymin[j] && y[i] < gymax[j] ) {
		location[i] = location[i] | LOCATION_AIR;
		break;
	    }
        }
    }

    free(gxmin);
    free(gxmax);
    free(gymin);
    free(gymax);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
