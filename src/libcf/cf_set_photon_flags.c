/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_set_photon_flags(fitsfile *infits, long nevents,
 *                                  float *photon_time, 
 *				    float *photon_weight,
 *				    unsigned char *photon_status,
 *				    unsigned char *photon_locflag,
 *                                  long nseconds, float *timeline_time,
 *                                  unsigned char *timeline_status);
 *
 * Description: Copy status flags from the timeline table to the photon list.
 *		Update header keywords.
 *
 *		For histogram data, we set only the DAY flag.
 *		Other situations are treated as follows:
 *
 *		DAY      All photons flagged if exposure >10% day
 *
 *		LIMB     True for all if true for one.
 *		SAA      Do not flag photons or modify EXPTIME.
 *		BRST     Instead, set header keyword EXP_STAT.
 *			 Ignore all violations during the 
 *			 first and last 20 seconds of an exposure.
 *
 *		HV	 Modify EXPTIME, but do not flag photons.
 *		JITR     Ditto.
 *
 *		OPUS     Ignore for now.
 *		USER	 Ignore for now.
 *
 * Arguments:   fitsfile *infits        Input FITS file pointer
 *              long     nevents        Number of points in the photon list
 *              float    *photon_time	Time array of photon events
 *              float    *photon_weight Scale factor for each photon
 *              unsigned char  *photon_status	Status flags for photon events
 *              unsigned char  *photon_locflag  Location flags for photon events
 *              long     nseconds       Number of points in the timeline table
 *              float    *timeline_time Time array of timeline list
 *              unsigned char  *timeline_status Status flags for timeline list
 *
 * Calls:
 *
 * Returns:	0 on success
 *
 * History:	10/31/02   1.1   jch    Initial coding
 * 		12/20/02   1.3   jch    Make sure the flags are "unsigned" char
 *		05/09/03   1.4   wvd	Change EXP* keywords to type long.
 *					Consider DAYNIGHT keyword in evaluation
 *					of NBADEVNT and EXP_BAD.
 *					Use frame_tolerance in time comparison.
 *					Use photon_locflag in counting NBADEVNT
 *               05/20/03  1.5   rdr    Add call to cf_proc_check
 *               05/29/03  1.6   wvd    Modify to handle HIST data:
 *					Use single flag for all photon events.
 *					Ignore LIMB, SAA, and BRST flags when
 *					counting EXP_BAD.  Scale NBADEVNT by
 *					photon_weight[i]/TOT_DEAD.
 *               06/02/03  1.7   wvd    Properly deal with GTI's.
 *               07/14/03  1.8   wvd    Read DAYNIGHT from SCRN_CAL and 
 *					populate IDF header keyword.
 *					Use FRAME_TOLERANCE from calfuse.h.
 *               09/15/03  1.9   wvd    Copy flags only to photons whose times
 *					are within 1s of timeline table time.
 *               10/02/03  1.10  wvd    Change timeline table: no breaks in
 *					time; add TEMPORAL_OPUS flag.
 *               10/13/03  1.11  wvd    Use same scheme for calculating
 *					EXP_BAD and exp_screen.
 *               04/20/04  1.12  bjg    Cosmetic change to prevent 
 *                                      compiler warning with gcc -Wall 
 *               05/20/04  1.13  wvd	Don't flag HIST photons as bad
 *					when HV or JITR are bad.
 *					We'll reduce EXPTIME instead.
 *               06/01/04  1.14  wvd	For HIST data, set only DAY flag.
 *					Write other flags to EXP_STAT keyword.
 *               06/28/04  1.15  wvd	If > 90% of exposure is lost to a
 *					burst, jitter, etc., write cause to
 *					trailer file.
 *               03/24/05  1.16  wvd	If EXP_STAT != 0, don't change it.
 *               01/01/07  1.17  wvd	Remove mention of FIFO flag.
 *               07/18/08  1.18  wvd	If EXP_STAT = 2, target is bright
 *					earth or airglow. Set limb-angle 
 *					flags to zero.
 *               08/22/08  1.19  wvd	Compare EXP_STAT to TEMPORAL_LIMB, 
 *					rather than a particular value.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include "calfuse.h"

#define HIST_PAD 20	/* Ignore first and last 20 s of each exposure. */

/*
 *  Compute integration time lost to a particular screening step.
 *  Ignore times outside of OPUS-defined good-time intervals.
 */
static long
time_loss(long nseconds, unsigned char *timeline_status,
	unsigned char criterion)
{
    long exp_screen = 0, j;

    for (j = 0; j < nseconds; j++) {
	if (timeline_status[j] & TEMPORAL_OPUS) continue;
	if (timeline_status[j] & criterion) exp_screen++;
    }

    return(exp_screen);
}


/*
 *  Compute integration time lost to a particular screening step,
 *  but ignore first and last HIST_PAD seconds of exposure.
 */
static long
time_loss_pad(long nseconds, unsigned char *timeline_status,
	unsigned char criterion)
{
    long exp_screen = 0, j;

    for (j = HIST_PAD; j < nseconds - HIST_PAD; j++) {
	if (timeline_status[j] & TEMPORAL_OPUS) continue;
	if (timeline_status[j] & criterion) exp_screen++;
    }
    return (exp_screen);
}


/*
 *  In HIST mode, generate a status flag good for all photons.
 *  Except for D/N flag, we ignore first and last 20 seconds of exposure.
 *  D/N flag is written to TEMPORAL_MASK.
 *  LIMB/SAA/BRST flags are written to EXP_STAT.
 */
static void
generate_temporal_mask(float rawtime, long nseconds, float *timeline_time,
	unsigned char *timeline_status, unsigned char *TEMPORAL_MASK,
	unsigned char *EXP_STAT, char *comment)
{
    long exp_day;

    *EXP_STAT = *TEMPORAL_MASK = 0;

   /*
    *  Day if exp_day > 10% of total exposure time.
    */
    exp_day=time_loss(nseconds, timeline_status, TEMPORAL_DAY);
    if (exp_day * 10 > rawtime) {
	*TEMPORAL_MASK = TEMPORAL_DAY;
	sprintf(comment, "Exposure more than 10 percent day");
	cf_verbose(1, comment);
    }
   /*
    * If exposure is less than 40 seconds long, we're done.
    */
    if (rawtime < 2. * HIST_PAD) return;

   /*
    * Are there limb-angle violations?
    */
    if (time_loss_pad(nseconds, timeline_status, TEMPORAL_LIMB)) {
    	*EXP_STAT |= TEMPORAL_LIMB;
	sprintf(comment, "Data suspect: Limb-angle violation");
	cf_verbose(1, comment);
    }
   /*
    * Do we ever enter the SAA?
    */
    if (time_loss_pad(nseconds, timeline_status, TEMPORAL_SAA)) {
	*EXP_STAT |= TEMPORAL_SAA;
	sprintf(comment, "Data suspect: SAA violation");
	cf_verbose(1, comment);
    }
   /*
    * Are there any bursts?
    */
    if (time_loss_pad(nseconds, timeline_status, TEMPORAL_BRST)) {
	*EXP_STAT |= TEMPORAL_BRST;
	sprintf(comment, "Data suspect: Burst detected");
	cf_verbose(1, comment);
    }

   /*
    * Is the detector voltage ever too low?
    * If so, just issue a warning; we'll modify the exposure time later.
    */
    if (time_loss_pad(nseconds, timeline_status, TEMPORAL_HV))
        cf_verbose(1, "EXPTIME modified: Detector voltage low.");
   /*
    * Was any time rejected by the jitter routine?
    * If so, just issue a warning; we'll modify the exposure time later.
    */
    if (time_loss_pad(nseconds, timeline_status, TEMPORAL_JITR))
        cf_verbose(1, "EXPTIME modified: Target out of aperture.");

    return;
}


int
cf_set_photon_flags(fitsfile *infits, long nevents, float *photon_time,
	    float *photon_weight, unsigned char *photon_status,
	    unsigned char *photon_locflag, long nseconds,
	    float *timeline_time, unsigned char *timeline_status)
{
    char CF_PRGM_ID[] = "cf_set_photon_flags";
    char CF_VER_NUM[] = "1.19";

    char  comment[FLEN_COMMENT], daynight[FLEN_KEYWORD]; 
    char  instmode[FLEN_KEYWORD], scrnfile[FLEN_KEYWORD];
    long  i, j, k;

    long nbadevnt=0;	/* number of events deleted in screening */
    long exp_bad=0;	/* integration time lost to screening */
    long exp_brst=0;	/* integration time lost to event bursts */
    long exp_hv=0;	/* integration time lost to low detector voltage */
    long exp_jitr=0;	/* integration time lost to jitter */
    long exp_limb=0;	/* integration time with low limb angle */
    long exp_saa=0;	/* integration time while in SAA */
    long exp_day=0;	/* integration time during day after screening */
    long expnight=0;	/* integration time during night after screening */
    int  day = FALSE, night = FALSE;
    int  errflg = 0;	/* value returned by cf_proc_check */
    int	 status = 0;	/* CFITSIO status */
    int  expstat;	/* initial value of EXP_STAT keyword */
    float deadtime, rawtime;
    fitsfile *scrnfits;
    unsigned char EXP_STAT, TEMPORAL_MASK;

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID) )) return errflg;

    /* Read and interpret DAYNIGHT keyword. Copy to IDF header. */

    FITS_read_key(infits, TSTRING, "SCRN_CAL", scrnfile, NULL, &status);
    FITS_open_file(&scrnfits, cf_parm_file(scrnfile), READONLY, &status);
    FITS_read_key(scrnfits, TSTRING, "DAYNIGHT", daynight, NULL, &status);
    FITS_close_file(scrnfits, &status);

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

    /*
     *  Read EXP_STAT, RAWTIME, INSTMODE and TOT_DEAD keywords.
     */
    FITS_read_key(infits, TINT, "EXP_STAT", &expstat, NULL, &status);
    FITS_read_key(infits, TFLOAT, "RAWTIME", &rawtime, NULL, &status);
    FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(infits, TFLOAT, "TOT_DEAD", &deadtime, NULL, &status);

    /*
     *  If EXP_STAT = TEMPORAL_LIMB, target is bright earth or airglow.  
     *  Compute EXP_LIMB, then set limb-angle flags to 0.
     */
    exp_limb=time_loss(nseconds, timeline_status, TEMPORAL_LIMB);
    if (expstat == (int) TEMPORAL_LIMB) 
	for (i = 0; i < nseconds-1; i++) 
	    timeline_status[i] &= ~TEMPORAL_LIMB;

    /*
     *  If in HIST mode:
     *  Generate TEMPORAL_MASK
     *  Copy status flags from TEMPORAL_MASK to photon list.
     *  Count bad photons.
     */
    if (!strncmp(instmode, "HIST", 4)) {
	generate_temporal_mask(rawtime, nseconds, timeline_time,
		timeline_status, &TEMPORAL_MASK, &EXP_STAT, comment);
	for (i = 0; i < nevents; i++) {
            photon_status[i] = TEMPORAL_MASK;
           /*
            *  A bad photon is one for which
            *  - something other than AIRGLOW is wrong; or
            *  - you don't want day, but this one is; or
            *  - you don't want night, but this one is.
            */
            if ((photon_locflag[i] & ~LOCATION_AIR) ||
                (!day   && (photon_status[i] & TEMPORAL_DAY)) ||
                (!night && (photon_status[i] ^ TEMPORAL_DAY)))
                nbadevnt += (long)(photon_weight[i]/deadtime + 0.5);
        }
	if (EXP_STAT && !expstat)
	    FITS_update_key(infits, TBYTE, "EXP_STAT", &EXP_STAT,
	    comment, &status);
    }

    /* 
     *  If in TTAG mode:
     *  Copy status flags from timeline table to photon list.
     *  Confirm that photon time is included in timeline table.
     *  Count bad photons.
     */
     else {
	for (i = k = 0; i < nevents; i++) {
            while (photon_time[i] > timeline_time[k+1]-FRAME_TOLERANCE &&
		k < nseconds-1) { k++; }

            if (photon_time[i] - timeline_time[k] < 2.0)
            	photon_status[i] = timeline_status[k];
	    else
		cf_if_warning("Time %g not included in timeline table",
		    photon_time[i]);

           /*
            *  A bad photon is one for which 
	    *  - something other than DAYNIGHT is wrong; or
	    *  - something other than AIRGLOW is wrong; or
	    *  - you don't want day, but this one is; or
	    *  - you don't want night, but this one is.
            */
            if ((photon_status[i] & ~TEMPORAL_DAY) ||
	        (photon_locflag[i] & ~LOCATION_AIR) ||
                (!day   && (photon_status[i] & TEMPORAL_DAY)) ||
                (!night && (photon_status[i] ^ TEMPORAL_DAY)))
	        nbadevnt++;
	}
     }

    /*
     *  Now count day, night, and bad seconds.
     *  If in HIST mode, mask out LIMB, SAA, and BRST flags.
     *  Exclude OPUS-defined bad times.
     */
    TEMPORAL_MASK = TEMPORAL_DAY;

    if (!strncmp(instmode, "HIST", 4)) {
        TEMPORAL_MASK |= TEMPORAL_LIMB;
        TEMPORAL_MASK |= TEMPORAL_SAA;
        TEMPORAL_MASK |= TEMPORAL_BRST;
    }

    for (j = 0; j < nseconds; j++) {
	if (timeline_status[j] & TEMPORAL_OPUS)
		continue;
	if (timeline_status[j] & ~TEMPORAL_MASK)
		exp_bad++;
	else if (timeline_status[j] & TEMPORAL_DAY)
		exp_day++;
	else
		expnight++;
    }
    /* 
     *  If user rejects either day or night data, add to EXP_BAD.
     *  If user rejects night data, set EXPNIGHT = 0.
     */
     if (!day) {
	exp_bad += exp_day;
     }
     if (!night) {
	exp_bad += expnight;
	expnight = 0;
     }

    /* 
     *  Compute EXP_BRST 
     */
    exp_brst=time_loss(nseconds, timeline_status, TEMPORAL_BRST);
    /* 
     *  Compute EXP_HV
     */
    exp_hv=time_loss(nseconds, timeline_status, TEMPORAL_HV);
    /* 
     *  Compute EXP_JITR 
     */
    exp_jitr=time_loss(nseconds, timeline_status, TEMPORAL_JITR);
    /* 
     *  Compute EXP_LIMB	(already done)
	exp_limb=time_loss(nseconds, timeline_status, TEMPORAL_LIMB);
     */
    /*
     *  Compute EXP_SAA 
     */
    exp_saa=time_loss(nseconds, timeline_status, TEMPORAL_SAA);

    /*
     *  If data are in time-tag mode and more than 90% of time was
     *  was rejected by screening, say why. 
     */
    if (!strncmp(instmode, "TTAG", 4)) {
	if (exp_brst > 0.9 * rawtime)
	    cf_verbose(1, "%d sec lost to bursts.", exp_brst);
	if (exp_hv > 0.9 * rawtime)
	    cf_verbose(1, "%d sec lost to low voltage.", exp_hv);
	if (exp_jitr > 0.9 * rawtime)
	    cf_verbose(1, "%d sec lost to jitter.", exp_jitr);
	if (exp_limb > 0.9 * rawtime)
	    cf_verbose(1, "%d sec lost to low limb angle.", exp_limb);
	if (exp_saa > 0.9 * rawtime)
	    cf_verbose(1, "%d sec lost to SAA.", exp_saa);
    }

    /* 
     *  Update the header keywords 
     */
    FITS_update_key(infits, TLONG, "NBADEVNT", &nbadevnt, NULL, &status);
    FITS_update_key(infits, TLONG, "EXP_BAD",  &exp_bad,  NULL, &status);
    FITS_update_key(infits, TLONG, "EXP_BRST", &exp_brst, NULL, &status);
    FITS_update_key(infits, TLONG, "EXP_HV",   &exp_hv,   NULL, &status);
    FITS_update_key(infits, TLONG, "EXP_JITR", &exp_jitr, NULL, &status);
    FITS_update_key(infits, TLONG, "EXP_LIM",  &exp_limb, NULL, &status);
    FITS_update_key(infits, TLONG, "EXP_SAA",  &exp_saa,  NULL, &status);
    FITS_update_key(infits, TLONG, "EXPNIGHT", &expnight, NULL, &status);
    FITS_update_key(infits, TSTRING, "DAYNIGHT", &daynight, NULL, &status);

    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
