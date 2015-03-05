/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_set_good_time_intervals(fitsfile *infits, long nseconds,
 *                          float *timeline_time, 
 *			    unsigned char *timeline_status, GTI *gti);
 *
 * Description: Compute the new good time intervals based on the timeline 
 *              status flag.
 *
 * Arguments:   fitsfile   *infits        	Input FITS file pointer
 *              long       nseconds       	Length of timeline table
 *              float      *timeline_time 	Time array of timeline list
 *     	   unsigned char   *timeline_status    	Status array of timeline list
 *              GIT        *gti           	Good time interval(s)
 *
 * Calls:
 *
 * Returns:     0 on success
 *              -1 on error
 *
 * History:     11/01/02   1.1   jch    Initial coding
 * 		12/20/02   1.3   jch    Make sure that flags are unsigned char
 *		05/09/03   1.4   wvd	Include D/N screening in GTI's.
 *					Calculate EXPTIME as SUM (GTI's).
 *					Update header keyword EXPTIME.
 *              05/20/03   1.5   rdr    Add call to cf_proc_check
 *              05/29/03   1.6   wvd    Allow lower-case value of DAYNIGHT.
 *					Correct bug in calculation of gti.stop 
 *					If in HIST mode, ignore LIMB, SAA,
 *					and BRST flags.
 *              09/10/03   1.7   wvd    Change calfusettag.h to calfuse.h
 *              10/13/03   1.8   wvd    Modify for new timeline table format.
 *					No time skips.
 *              04/20/04   1.9   bjg    Remove unused variables
 *              12/01/04   1.10  bjg    Set gti start and stop to zero 
 *                                      for the single GTI when there are no
 *                                      Good Times.
 *              08/23/07   1.11  wvd    Change malloc to cf_malloc.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_set_good_time_intervals(fitsfile *infits, long nseconds,
	float *timeline_time, unsigned char *timeline_status, GTI *gti)
{
    char CF_PRGM_ID[] = "cf_set_good_time_intervals";
    char CF_VER_NUM[] = "1.11";

    char   daynight[FLEN_VALUE], instmode[FLEN_VALUE];
    int	   in_bad;  		   /* True between good-time intervals */
    int    day=FALSE, night=FALSE, errflg=0, status=0;
    long   j, jj;
    float  exptime=0.;
    double last_good_time=-99;
    unsigned char TEMPORAL_MASK;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID) )) return errflg;

    /* Read and interpret DAYNIGHT keyword. */

    FITS_read_key(infits, TSTRING, "DAYNIGHT", daynight, NULL, &status);
    if (!strncmp(daynight, "D", 1) || !strncmp(daynight, "d", 1)) {
        day = TRUE;
    }
    else if (!strncmp(daynight, "N", 1) || !strncmp(daynight, "n", 1)) {
        night = TRUE;
    }
    else if (!strncmp(daynight, "B", 1) || !strncmp(daynight, "b", 1)) {
        day = night = TRUE;
    }
    else
        cf_if_error("Unknown DAYNIGHT keyword value: %s", daynight);
    cf_verbose(3, "DAYNIGHT flag set to %s", daynight);

    /* 
     *  Read INSTMODE keyword.  If in HIST mode, mask out
     *  LIMB, SAA, and BRST flags.  TEMPORAL_DAY is the default.
     */
    TEMPORAL_MASK = TEMPORAL_DAY;
    FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL, &status);
    if (!strncmp(instmode, "HIST", 4)) {
	TEMPORAL_MASK |= TEMPORAL_LIMB;
	TEMPORAL_MASK |= TEMPORAL_SAA;
	TEMPORAL_MASK |= TEMPORAL_BRST;
    }

    /*
     *  Determine limits of good-time intervals.
     */
    gti->start = (double *) cf_malloc(sizeof(double)*nseconds);
    gti->stop  = (double *) cf_malloc(sizeof(double)*nseconds);

    in_bad = TRUE;
    jj = -1;
    for (j = 0; j < nseconds; j++) {
        /*
         *  If we are in a good time interval and were in a bad one,
	 *  increment jj and set gti->start.
         */
        if (!(timeline_status[j] & ~TEMPORAL_MASK) &&
            ((day   && (timeline_status[j] & TEMPORAL_DAY)) ||
             (night && (timeline_status[j] ^ TEMPORAL_DAY)))) {
	    if (in_bad) {
		jj++;
		gti->start[jj] = timeline_time[j];
		in_bad = FALSE;
	    }
        /*
         *  If we are in a good time interval and were in a good one,
	 *  update last_good_time.
         */
	    last_good_time = timeline_time[j];
	}
        /*
         *  If we are in a bad time interval and were in a good one,
	 *  set in_bad to TRUE and update gti->stop.
         */
	else if (in_bad == FALSE) {
	    gti->stop[jj] = timeline_time[j];
	    in_bad = TRUE;
	}
    }
    gti->ntimes = jj + 1;

    /* If the entire exposure is bad, return a zero-length GTI. */
    if (gti->ntimes == 0) {
	gti->ntimes = 1;
	gti->start[0] = gti->stop[0] = 0;
    }

    /* Calculate total exposure time and write it to file header */
    for (j = 0; j < gti->ntimes; j++)
	exptime += gti->stop[j] - gti->start[j];
    FITS_update_key(infits, TFLOAT, "EXPTIME", &exptime, NULL, &status);

    /* Results to trailer file. */
    cf_verbose(2, "GTI's:\t START\t  STOP");
    for (j = 0; j < gti->ntimes; j++)
        cf_verbose(2, "%3d\t%6.0f\t%6.0f", j, gti->start[j], gti->stop[j]);
    cf_verbose(2, "Total EXPTIME = %6.0f", exptime);

    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
