/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_modify_hist_times (fitsfile *infits, long nevents,
 *                          float *time, GTI *gti);
 *
 * Description: For histogram data, set photon-arrival times to the midpoint
 *              of the longest good-time interval.
 *
 * Arguments:   fitsfile   *infits        	Input FITS file pointer
 *              long       nevents       	Number of photon events
 *              float      *time 		Time array for photon list
 *              GIT        *gti           	Good time interval(s)
 *
 * Calls:
 *
 * Returns:     TRUE if time array is changed.
 *              FALSE if time array is not changed.
 *
 * History:     06/02/04   1.1   wvd    Initial coding
 *		06/07/04   1.2   wvd	Use standard return values.
 *		02/17/05   1.3   wvd	Initalize jmax to 0.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_modify_hist_times (fitsfile *infits, long nevents, float *time, GTI *gti)
{
    char CF_PRGM_ID[] = "cf_modify_hist_times";
    char CF_VER_NUM[] = "1.3";

    int    errflg=0, status=0;
    long   j, jmax=0;
    float  exptime, gti_time, max_time, photon_time, rawtime;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* If data were not taken in HIST mode, exit now. */
    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    /* Read header keywords. */
    FITS_read_key(infits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    FITS_read_key(infits, TFLOAT, "RAWTIME", &rawtime, NULL, &status);

    /* 
     *  If the entire exposure was rejected by the screening routines,
     *	we don't bother changing photon-arrival times.  If no time was
     *	lost to screening, the default arrival times are OK.  Exit now.
     */
    if (exptime < 1. || exptime > rawtime - 1.) {
	cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
	return -1;
    }

    /* Determine which is the longest good-time interval. */
    max_time = 0.;
    for (j = 0; j < gti->ntimes; j++) {
	if ((gti_time = gti->stop[j] - gti->start[j]) > max_time) {
	    max_time = gti_time;
	    jmax = j;
	}
    }

    /* Set all photon-arrival times to midpoint of longest GTI. */
    photon_time = (gti->start[jmax] + gti->stop[jmax]) / 2.;
    for (j = 0; j < nevents; j++)
	time[j] = photon_time;
    cf_verbose(1, "Setting photon-arrival times to %.1f", photon_time);

    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
