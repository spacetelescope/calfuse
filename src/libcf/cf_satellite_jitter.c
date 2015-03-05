/*****************************************************************************
 *              Johns Hopkins University
 *              Center for Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_satellite_jitter infile, nevents, time, x, y, channel,
 *              	nsec, ttime, status_flag)
 *
 * Description: Correct photon coordinates for spacecraft motion.
 *
 *		Operates only on point-source observations made in TTAG mode.
 *		Only photon events in the target aperture are corrected.
 *		Minimum trustworthy TRKFLG value is read from parameter file.
 *
 *              Key to new TRKFLG values:
 *      	5 = dx, dy from known FPDs and q_cmd from telemetry
 *      	4 = good dx, dy from ACS q_est and q_cmd from telemetry
 *      	3 = maybe good dx, dy from ACS q_est
 *      	2 = dx, dy from known FPDs, but q_cmd from FITS header coords
 *      	1 = dx, dy from ACS q_est, but q_cmd from FITS header coordis
 *      	0 = no pointing information (missing telemetry)
 *             -1 = pointing is assumed to be bad (never achieved known track)
 *
 * Arguments:   fitsfile  infile        Pointer to the IDF
 *              long      nevents       The number of events
 *              float     *time         An array of event times
 *              float     *x            An array of event X positions
 *              float     *y            An array of event Y positions
 *     unsigned char      *channel      Channel assignment for each photon
 *              long      nsec          Number of seconds in the timeline
 *              long      *ttime        time (sec) of each point in the timeline
 *     unsigned char      *status_flag  Timeline status flags
 *
 * Returns:	0 upon success
 *
 * History:     09/12/02   v1.1   RDR   Adapted from cf_ttag_jitter
 *              10/02/02   v1.2   RDR   Modified to be compatable with
 *                                      new IDF file format
 *              12/10/02   v1.3   RDR   Adjusted method of populating the 
 *                                      timeline jitter status flag. 
 *              12/13/02   v1.4   RDR   Exit gracefully if there is no jitr file 
 *              03/01/03   v1.6   wvd   Correct use of pointer in FITS_read_key
 *              04/07/03   v1.7   wvd   Confirm that TRKFLG <= 5.
 *					Replace printf with cf_verbose.
 *              04/14/03   v1.8   rdr   Changed flag from char to unsigned char
 *              04/22/03   v1.9   wvd   Move only photons in target aperture.
 *					Do not assume that ttime is continuous.
 *              04/29/03   v1.10  wvd   Allow possiblity that jitter
 *                                      filename is in lower-case letters.
 *              05/20/03   v1.11  rdr   Added call to cf_proc_check
 *              05/22/03   v1.12  wvd   Direct cf_error_init to stderr
 *              07/11/03   v1.13  rdr   Set all points as bad if JIT_STAT=1
 *                                       (no reference position found)
 *              08/06/03   v1.14  wvd   Improve documentation, delete
 *					comparison with GTI's, change channel
 *					array to unsigned char.
 *              08/29/03   v1.15  wvd   Set all points as bad if JIT_STAT != 0
 *              09/18/03   v1.16  wvd   Print name of jitter file only if
 *					the file exists.
 *              10/03/03   v1.17  wvd   Check HKEXISTS before looking for 
 *					jitter file.
 *              02/12/04   v1.18  wvd   Make interpolation immune
 *					to errors in jitter time array.
 *              02/24/04   v1.19  rdr   Change limits of acceptable jitter
 *					values
 *              03/01/04   v1.20  rdr   Populate the PLATESC keyword in the
 *					input file header
 *              06/01/04   v1.21  bjg   Exit when keyword JIT_STAT not found.
 *              06/07/04   v1.22  wvd   Return EXP_JITR to calling routine.
 *		01/26/05   v1.23  wvd   Fix tiny bug in a verbose stmt.
 *		04/05/05   v1.24  wvd   If JIT_STAT != 0, APERTURE = RFPT,
 *					or target = airglow, exit without 
 *					flagging any times as bad.
 *					Use new TRKFLG values (see above).
 *					Delete adjust_trkflg subroutine and
 *					the J_TRKFLG and JTIMEGAP keywords.
 *					Use cf_read_col to read jitter file.
 *		07/06/05   v1.25  wvd   If JIT_VERS < 2.0, issue a warning.
 *		11/29/05   v1.26  wvd	Move screening functions into 
 *					cf_screen_jitter.
 *              08/21/06   v1.27  wvd   Use EXPSTART in jitter and data
 *                                      headers to correct jitter time array.
 *              08/25/06   v1.28  wvd   Change meaning of TRKFLG values.
 *                                      TRKFLG = 3 is only used internally.
 *              11/06/06   v1.29  wvd   Change meaning of TRKFLG values.
 *					Read APER_COR keyword in file header.
 *					Read min TRKFLG value from parmfile.
 *					Don't have to test observing mode.
 *              03/23/07   v1.30  wvd   Store offset between the jitter and 
 *					data values of EXPSTART in variable
 *					time_diff.
 *              04/07/07   v1.31  wvd   Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "calfuse.h"

int
cf_satellite_jitter(fitsfile *infits, long nevents, float *time,
	float *x, float *y, unsigned char *channel, long nsec,
	float *ttime, unsigned char *status_flag)
{
    char CF_PRGM_ID[] = "cf_satellite_jitter"; 
    char CF_VER_NUM[] = "1.31";

    fitsfile  *jitfits, *parmfits;
    char  det[FLEN_VALUE];
    char  jitr_cal[FLEN_VALUE], parmfile[FLEN_VALUE];
    char  comment[FLEN_COMMENT];
    short c, *trkflg, trkflg_min;
    int   active_ap[2], extended, status=0, hdutype;
    long  j, k, njitr, nrej, *time_jit, time_diff;
    double data_expstart, jitr_expstart;
    float *dx_jit, *dy_jit, x_factor, y_factor;
    float *dx_tt, *dy_tt;
    float factor[2][4] = {{1.743,1.743,-1.743, -1.743},
			  {1.154,1.163,-0.6335,-0.581}};

    /* Initialize error checking */ 
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Read exposure start time from IDF header */
    FITS_read_key(infits, TDOUBLE, "EXPSTART", &data_expstart, NULL, &status);

    /* If APER_COR is not set to "COMPLETE", exit now. */
    FITS_read_key(infits, TSTRING, "APER_COR", comment, NULL, &status);
    if (strncmp(comment, "COMPLETE", 8)) {
        cf_verbose(1, "APER_COR = %s.  Omitting photon shift.", comment);
	cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
	return (status);
    }

    /* If target is an extended source, exit now. */
    extended = cf_source_aper(infits, active_ap);
    if (extended) {
	cf_verbose(1, "Extended source.  Omitting photon shift.");
	cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
	return (status);
    }

    /* Determine the scale factor to convert arcseconds to pixels. */
    FITS_read_key(infits, TSTRING, "DETECTOR", det, NULL, &status);
    c = det[0] =='2'? 2: 0;
    if (det[1] == 'B')
	c++;
    x_factor=factor[0][c];
    y_factor=factor[1][c];

    /* Update the PLATESC (plate scale) keyword in the header */
    FITS_update_key(infits, TFLOAT, "PLATESC", &y_factor, NULL, &status);

    cf_verbose(3, "Arcsec to pixels: X = %5.3f\tY = %5.3f",
	x_factor, y_factor);

    /* Try to open the jitter file. */
    FITS_read_key(infits, TSTRING, "JITR_CAL", jitr_cal, NULL, &status);
    fits_open_file(&jitfits, jitr_cal, READONLY, &status);

    /* If the jitter file is not found, try a lower-case filename. */
    if (status != 0) {
        status = 0;
        jitr_cal[0] = (char) tolower(jitr_cal[0]);
        fits_open_file(&jitfits, jitr_cal, READONLY, &status);
    }

    /* Read exposure start time from jitter file header */
    FITS_read_key(jitfits, TDOUBLE, "EXPSTART", &jitr_expstart, NULL, &status);

    /* Allocate space for correction arrays */
    dx_tt = (float *) cf_calloc(nsec, sizeof(float));
    dy_tt = (float *) cf_calloc(nsec, sizeof(float));

    /* Read the jitter file. */
    FITS_movabs_hdu(jitfits, 2, &hdutype, &status);
    njitr=cf_read_col(jitfits, TLONG, "TIME", (void **) &time_jit);
    njitr=cf_read_col(jitfits, TFLOAT, "DX", (void **) &dx_jit);
    njitr=cf_read_col(jitfits, TFLOAT, "DY", (void **) &dy_jit);
    njitr=cf_read_col(jitfits, TSHORT, "TRKFLG", (void **) &trkflg);
    FITS_close_file(jitfits, &status);

    /* Read minimum acceptable TRKFLG value from PARM_CAL file */
    FITS_read_key(infits, TSTRING, "PARM_CAL", parmfile, NULL, &status);
    FITS_open_file(&parmfits, cf_parm_file(parmfile), READONLY, &status);
    FITS_read_key(parmfits, TSHORT, "TRKFLG", &trkflg_min, NULL, &status);
    FITS_close_file(parmfits, &status);

    /* Correct for difference in data and jitter EXPSTART times */
    time_diff = cf_nlong((data_expstart - jitr_expstart) * 86400.);
    for (j=0; j< njitr; j++) time_jit[j] -= time_diff;

    /*
     * Map times in the jitter file to times in the timeline table
     * and calculate the jitter offsets.
     */
    for (j=k=0; k < nsec; k++) {
	while ((float) time_jit[j+1] < ttime[k] && j+1 < njitr) j++;
	if (trkflg[j] >= trkflg_min) {
	    dx_tt[k] = x_factor * dx_jit[j];
	    dy_tt[k] = y_factor * dy_jit[j];
	}
    }
    /*
     *  Map times in the timeline table to the times associated
     *  with individual photons and apply the appropriate shifts.
     *  Move only photons in the active aperture and for which the
     *  timeline status jitter flag has not been set.
     *  Record number of rejected events.
     *  Note: active_ap[0] = lif, active_ap[1] = sic
     */
    nrej = 0;
    for (j=k=0; j<nevents; j++) {
	while (ttime[k+1] - FRAME_TOLERANCE < time[j] && k+1 < nsec) k++;
	    if (channel[j] == active_ap[0] || channel[j] == active_ap[1]) {
		if (status_flag[k] & TEMPORAL_JITR) nrej++;
		else {
		   x[j] += dx_tt[k];
		   y[j] += dy_tt[k];
		}
	    }
    }

    /* Write to trailer file */
    cf_verbose(2, "Rejected %ld photons from target aperture", nrej);

    /* Deallocate storage space */
    free(time_jit);
    free(dx_jit);
    free(dy_jit);
    free(trkflg);
    free(dx_tt);
    free(dy_tt);
   
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");

    return (status);
}
