/*****************************************************************************
 *              Johns Hopkins University
 *              Center for Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_screen_jitter (infile, nsec, ttime, status_flag)
 *
 * Description: Flag times when target is out of the aperture or pointing
 *		is known to be bad.  Read the minimum trustworthy TRKFLG 
 *		value and aperture limits from the parameter file.
 *		
 *              Key to new TRKFLG values:
 *      	5 = dx, dy from known FPDs and q_cmd from telemetry
 *      	4 = good dx, dy from ACS q_est and q_cmd from telemetry
 *      	3 = maybe good dx, dy from ACS q_est
 *      	2 = dx, dy from known FPDs, q_cmd from FITS header coordinates
 *      	1 = dx, dy from ACS q_est, q_cmd from FITS header coordinates
 *      	0 = no pointing information (missing telemetry)
 *             -1 = Pointing is assumed to be bad (never achieved known track)
 *
 * Arguments:   fitsfile  infile        Pointer to the IDF
 *              long      nsec          Number of seconds in the timeline
 *              long      *ttime        time (sec) of each point in the timeline
 *     unsigned char      *status_flag  Timeline status flags
 *
 * Returns:	0 upon success
 *
 * History:     11/29/05   v1.1   wvd   Adapted from cf_satellite_jitter
 *		03/28/06    1.2   wvd	Test both dx and dy.
 *		05/09/06    1.3   wvd	Apply jitter screening to RFPT obs.
 *		06/12/06    1.4   wvd	Call cf_verbose when jitter file not
 *					found.
 *		08/21/06    1.5   wvd	Use EXPSTART in jitter and data
 *					headers to correct jitter time array.
 *		08/25/06    1.6   wvd	Change meaning of TRKFLG values.
 *					TRKFLG = 3 is only used internally.
 *		11/06/06    1.7   wvd	Change meaning of TRKFLG values.
 *					Read minimum TRKFLG value and
 *					aperture limits from parmfile.
 *					Require jitter file version >= 3.0.
 *					Don't screen moving targets.
 *					Bail out only if JIT_STAT = 1.
 *		12/29/06    1.8   wvd	Exit if JIT_VERS < 3.0.
 *              03/23/07    1.9   wvd   Store offset between the jitter and 
 *					data values of EXPSTART in variable
 *					time_diff.
 *              04/07/07    1.10  wvd   Clean up compiler warnings.
 *              08/16/07    1.11  wvd   Don't apply jitter correction to
 *					observations taken after the 
 *					spacecraft lost pointing control.
 *              07/25/08    1.12  wvd   Use EXP_STAT keyword, rather than
 *					file name, to check for bright-earth
 *					or airglow observations.
 *              11/21/08    1.13  wvd   Don't apply jitter correction to
 *					ERO observations of Jupiter.
 *					
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "calfuse.h"

#define	JIT_VERS_MIN	3.0
#define	BAD_TRKFLG	-1
#define LAST_MJD	54300.

int
cf_screen_jitter(fitsfile *infits, long nsec, float *ttime,
	unsigned char *status_flag)
{
    char CF_PRGM_ID[] = "cf_screen_jitter"; 
    char CF_VER_NUM[] = "1.13";

    char  jitr_cal[FLEN_VALUE], parmfile[FLEN_VALUE];
    char  comment[FLEN_COMMENT], hkexists[FLEN_VALUE], jit_vers[FLEN_VALUE];
    char  mov_targ[FLEN_VALUE], program[FLEN_VALUE];
    int   errflg=0, status=0, expstat, jitter_status, hdutype;
    long  j, k, njitr, *time_jit, time_diff;
    short *trkflg, trkflg_min;
    double data_expstart, jitr_expstart;
    float *dx_jit, *dy_jit, dx_max, dy_max;
    fitsfile  *jitfits, *parmfits;

    /* Initialize error checking */ 
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Read exposure start time from IDF header */
    FITS_read_key(infits, TDOUBLE, "EXPSTART", &data_expstart, NULL, &status);

    /* Check whether program is appropriate for this data set. */
    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    /* Don't apply jitter screening to bright-earth or airglow observations. */
    FITS_read_key(infits, TINT, "EXP_STAT", &expstat, NULL, &status);
    if (expstat == (int) TEMPORAL_LIMB) {
        cf_verbose(1, "Bright-earth or airglow observation.  ",
	"No jitter correction.");
        cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
        return 0;
    }

    /* Don't apply jitter screening to moving targets. */
    FITS_read_key(infits, TSTRING, "MOV_TARG", mov_targ, NULL, &status);
    if (!strncmp(mov_targ, "M", 1)) {
        cf_verbose(1, "Moving-target observation.  No jitter correction.");
        cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
        return 0;
    }

    /* Don't apply jitter screening to ERO observations of Jupiter. */
    FITS_read_key(infits, TSTRING, "PRGRM_ID", program, NULL, &status);
    if (!strncmp(program, "X006", 4)) {
        cf_verbose(1, "Moving-target observation.  No jitter correction.");
        cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
        return 0;
    }

    /* Don't apply jitter screening after pointing control is lost. */
    if (data_expstart > LAST_MJD) {
        cf_verbose(1, "No pointing control.  No jitter correction.");
        cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
        return 0;
    }

    /* Check HKEXISTS before looking for jitter file. */
    FITS_read_key(infits, TSTRING, "HKEXISTS", hkexists, NULL, &status);
    if(!strncasecmp(hkexists, "N", 1) || !strncasecmp(hkexists, "n", 1)) {
        cf_verbose(1, "HKEXISTS = %s.  No jitter correction.", hkexists);
        cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
        return 0;
    }

    /* Try to open the jitter file. */
    FITS_read_key(infits, TSTRING, "JITR_CAL", jitr_cal, NULL, &status);
    fits_open_file(&jitfits, jitr_cal, READONLY, &status);

    /* If the jitter file is not found, try a lower-case filename. */
    if (status != 0) {
        status = 0;
        jitr_cal[0] = (char) tolower(jitr_cal[0]);
        fits_open_file(&jitfits, jitr_cal, READONLY, &status);
    }

    /* If the jitter file is still not found, exit the program */
    if (status != 0) {
        cf_verbose(1, "Jitter file not found.  No jitter correction.");
        cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
        return 0;
    }
    cf_verbose(3, "Jitter calibration file = %s", jitr_cal);

    /* Check status keyword in jitter file.  Exit if missing or non-zero. */
    fits_read_key(jitfits, TINT, "JIT_STAT", &jitter_status, comment, &status);
    if (status != 0) {
	status = 0;
	cf_verbose(1, "Keyword JIT_STAT not found. No jitter correction.");
	cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
	return 0;
    }
    if (jitter_status == 1) {
	cf_verbose(1, "Keyword JIT_STAT = %d.  No jitter correction.",
	jitter_status) ;
	cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
	return 0;
    }

    /* Read exposure start time from jitter file header */
    FITS_read_key(jitfits, TDOUBLE, "EXPSTART", &jitr_expstart, NULL, &status);

    /* Check JIT_VERS keyword.  Exit if missing or value too low. */
    fits_read_key(jitfits, TSTRING, "JIT_VERS", jit_vers, NULL, &status);
    if (status != 0 || atof(jit_vers) < JIT_VERS_MIN) {
	status = 0;
	cf_verbose(1, "Jitter file format is out of date. No jitter correction.");
	cf_proc_update(infits, CF_PRGM_ID, "SKIPPED");
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");
	return 0;
    }

    /* Read the jitter file. */
    FITS_movabs_hdu(jitfits, 2, &hdutype, &status);
    njitr=cf_read_col(jitfits, TLONG, "TIME", (void **) &time_jit);
    njitr=cf_read_col(jitfits, TFLOAT, "DX", (void **) &dx_jit);
    njitr=cf_read_col(jitfits, TFLOAT, "DY", (void **) &dy_jit);
    njitr=cf_read_col(jitfits, TSHORT, "TRKFLG", (void **) &trkflg);
    FITS_close_file(jitfits, &status);

    /* Correct for difference in data and jitter EXPSTART times */
    time_diff = cf_nlong((data_expstart - jitr_expstart) * 86400.);
    cf_verbose(2, "Offset between IDF and jitter EXPSTART values: %d sec",
	time_diff);
    for (j=0; j< njitr; j++) time_jit[j] -= time_diff;

    /* Read minimum TRKFLG value and aperture boundaries from PARM_CAL file */
    FITS_read_key(infits, TSTRING, "PARM_CAL", parmfile, NULL, &status);
    FITS_open_file(&parmfits, cf_parm_file(parmfile), READONLY, &status);
    FITS_read_key(parmfits, TSHORT, "TRKFLG", &trkflg_min, NULL, &status);
    FITS_read_key(parmfits, TFLOAT, "DX_MAX", &dx_max, NULL, &status);
    FITS_read_key(parmfits, TFLOAT, "DY_MAX", &dy_max, NULL, &status);
    FITS_close_file(parmfits, &status);

    /* Reject times when tracking is good, but target is out of aperture. */
    for (j=0; j< njitr; j++) {
	if (trkflg[j] >= trkflg_min  && 
		(dx_jit[j] < -dx_max || dx_jit[j] > dx_max ||
		 dy_jit[j] < -dy_max || dy_jit[j] > dy_max)) {
	    trkflg[j] = BAD_TRKFLG;
	}
    }

    /*
     * Map times in the jitter file to times in the timeline table.
     * If trkflg[j] = BAD_TRKFLG, set the jitter flag in the timeline
     * status array.
     */
    for (j=k=0; k < nsec; k++) {
	while ((float) time_jit[j+1] < ttime[k] && j+1 < njitr) j++;
	if (trkflg[j] == BAD_TRKFLG) {
	    status_flag[k] |= TEMPORAL_JITR;
	}
    }

    /* Deallocate storage space */
    free(time_jit);
    free(dy_jit);
    free(trkflg);
   
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");

    return (status);
}
