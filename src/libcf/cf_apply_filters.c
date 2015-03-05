/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_apply_filters(fitsfile *header, int tscreen, long nevents,
 *                      unsigned char *timeflags, unsigned char *loc_flags,
 *			long ntimes, unsigned char *gtiflags,
 *			float *dtime, float *ntime, long *ngood, long **index)
 *
 * Description: Select photons that satisfy filter criteria
 *
 * Arguments:   fitsfile  *header       Pointer to FITS file containing the
 *                                      header of the intermediate data file
 *              int       tscreen       Flag indicating whether to screen on
 *                                       time flags (yes if > 0)
 *              long      nevents       The number of events
 *              unsigned char *timeflags Photon event time flags
 *              unsigned char *loc_flags Photon event location flags
 *              long      ntimes        Number of entries in timeline table
 *              unsigned char *gtiflags Status flags from timeline table
 *              float     *dtime        Day time exposure length
 *              float     *ntime        Night time exposure length
 *              long      *ngood        Number of good photons
 *              long      **index       List of good photon indexes
 *
 * Return:      status
 *
 * History:     02/19/03   1.1    peb   Begin work
 *              03/20/03   1.2    peb   Fixed daytime and nighttime
 *                                      calculation and fill header with
 *                                      correct EXPTIME and EXPNIGHT value.
 *              04/08/03   1.3    wvd   Get timeline-table arrays from 
 *					calling routine.
 *		05/10/03   1.4    wvd	Read DAYNIGHT from file header;
 *					test only first character.
 *              05/20/03   1.5    rdr   Move to top level HDU before reading
 *                                      header keywords
 *              06/09/03   1.6    rdr   Ignore timing flags in the screening
 *                                      if tscreen is 0.
 *              10/06/03   1.8    wvd	Properly deal with new format of
 *					timeline table.
 *              12/08/03   1.9    wvd	If tscreen is 0, ignore temporal
 *					flags when calculating dtime & ntime.
 *              02/11/04   1.10   wvd	Correct error in calculation of
 *					dtime and ntime.
 *              06/02/04   1.11   wvd	Mask BRST, SAA, and LIMB flags
 *					when calculating exposure time
 *					for HIST data.
 *              06/10/04   1.12   wvd	Delete timeline times array from 
 *					argument list.
 *		01/27/05   1.13   wvd	Use TEMPORAL_MASK rather than
 *					TEMPORAL_DAY to screen photon status
 *					flags.  They differ only in HIST mode.
 *		07/18/08   1.14   wvd	If EXP_STAT = 2, target is bright
 *					earth or airglow.  Mask limb-angle flag.
 *		08/22/08   1.15   wvd	Compare EXP_STAT with TEMPORAL_LIMB, 
 *					rather than a particular value.
 *
 ****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "calfuse.h"

int
cf_apply_filters(fitsfile *header, int tscreen, long nevents, 
	unsigned char *timeflags, unsigned char *loc_flags,
	long ntimes, unsigned char *gtiflags,
	long *dtime, long *ntime, long *ngood, long **index)
{
    char CF_PRGM_ID[] = "cf_apply_filters";
    char CF_VER_NUM[] = "1.15";

    char daynight[FLEN_VALUE], instmode[FLEN_VALUE];
    int  day=0, expstat, night=0, status=0;
    long j;
    float exptime;
    unsigned char TEMPORAL_MASK;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* 
     *  Read INSTMODE keyword.  If in HIST mode, mask out
     *  LIMB, SAA, and BRST flags.  TEMPORAL_DAY is the default.
     */
    TEMPORAL_MASK = TEMPORAL_DAY;
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    if (!strncmp(instmode, "HIST", 4)) {
	TEMPORAL_MASK |= TEMPORAL_LIMB;
	TEMPORAL_MASK |= TEMPORAL_SAA;
	TEMPORAL_MASK |= TEMPORAL_BRST;
    }
    /*
     * If EXP_STAT = TEMPORAL_LIMB, target is bright earth or airglow.  
     * Mask limb-angle flag.
     */
    FITS_read_key(header, TINT, "EXP_STAT", &expstat, NULL, &status);
    if (expstat == (int) TEMPORAL_LIMB) TEMPORAL_MASK |= TEMPORAL_LIMB;

    /*
     *  Calculate day/night exposure time.
     *  If tscreen is 0, ignore the temporal flags (other than day/night).
     */
    *dtime = *ntime = 0.;
    for (j=0; j<ntimes-1; j++) {
	if ( !(gtiflags[j] & ~TEMPORAL_MASK) || !tscreen ) {
	   if (gtiflags[j] & TEMPORAL_DAY)
		*dtime += 1;
	   else
		*ntime += 1;
	}
    }
    cf_verbose(2, "Exposure time: day = %ld, night = %ld", *dtime, *ntime);
    /*
     *  Read EXPTIME and DAYNIGHT keywords from file header.
     */
    FITS_movabs_hdu(header, 1, NULL, &status) ;
    FITS_read_key(header, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    FITS_read_key(header, TSTRING, "DAYNIGHT", daynight, NULL, &status);

    /* Decode DAYNIGHT keyword */
    if (!strncmp(daynight, "D", 1) || !strncmp(daynight, "d", 1)) {
	day = 1;
	*ntime = 0.;
    }
    else if (!strncmp(daynight, "N", 1) || !strncmp(daynight, "n", 1)) {
	night = 1;
	*dtime = 0.;
    }
    else if (!strncmp(daynight, "B", 1) || !strncmp(daynight, "b", 1)) {
	day = night = 1;
    }
    else
	cf_if_error("Unknown DAYNIGHT keyword value");
    cf_verbose(3, "DAYNIGHT flag set to %s", daynight);

    /* Quick test of internal consistency */
    if (fabs(*dtime + *ntime - exptime) > exptime/10.)
	cf_if_warning("Sum of day and night time differs from EXPTIME by > 10%%");

    /*
     *  Apply screening to each photon event.
     *  First mask out airglow and day flags and test for zero
     *  (good) values.  Next, see if the user wants daytime,
     *  night time, or both and confirm that appropriate flag is set.
     */

    *index = (long *) cf_calloc(nevents, sizeof(long));

    if (tscreen) {
        for(*ngood=j=0; j<nevents; j++) {
	    if (!(loc_flags[j] & ~LOCATION_AIR) && 
	        !(timeflags[j] & ~TEMPORAL_MASK) &&
	        ((day   && (timeflags[j] & TEMPORAL_DAY)) || 
	         (night && (timeflags[j] ^ TEMPORAL_DAY)))) {
	        (*index)[(*ngood)++] = j;
	    }
	}
    }

    /*  If tscreen is 0, then screen only on the location flags.  */
    else {
	cf_verbose(1, "Ignoring time flags in the screening") ;
	for(*ngood=j=0; j<nevents; j++) 
	    if (!(loc_flags[j] & ~LOCATION_AIR)) 
	        (*index)[(*ngood)++] = j;
    }

    cf_verbose(2, "%ld good events", *ngood);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
