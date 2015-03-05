/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Description: Routines used by cf_ttag_init and cf_hist_init.
 *              
 * History:     07/15/03  v1.1	wvd   Move routines from cf_ttag_init
 *		07/24/03   1.2  wvd   If tsunrise and tsunset get huge, 
 *					store them as floats, not integers.
 *		07/30/03   1.5  wvd   Omit TZERO and TSCALE for these 
 *					timeline table entries:
 *					HIGH_VOLTAGE
 *					LIF_CNT_RATE
 *					SIC_CNT_RATE
 *					FEC_CNT_RATE
 *					AIC_CNT_RATE
 *					BKGD_CNT_RATE
 *					This will effectively convert them
 *					to arrays of shorts.
 *					Add 0.5 to all count rates in 
 *					anticipation of truncation.
 *		08/22/03   1.7  wvd   Populate MIN_LIMB keyword.
 *		09/03/03   1.8  wvd   Implement cf_verbose throughout.
 *				      Add EXP_JITR and EXPNIGHT to IDF header.
 *		09/22/03   1.9  wvd   Correct comment field for YCENT1
 *		10/02/03   1.10 wvd   Modify timeline table:
 *				      - no time gaps
 *				      - set TEMPORAL_OPUS for times outside
 *					of OPUS good-time intervals.
 *				      Delete cf_get_times.
 *				      Use cf_read_col in cf_get_gti
 *					and cf_get_geocorona.
 *				      Don't modify gtis[0] for HIST data.
 *		10/13/03   1.11 wvd   Add 1s to end of timeline table for
 *					TTAG data.
 *		10/15/03   1.12 wvd   Swap voltages for detectors 2A and 2B
 *					if housekeeping file was written
 *					before April of 2003.
 *		10/22/03   1.13 wvd   Compute EXPNIGHT for raw data file
 *					and write to output file header.
 *					Add YQUALLIF & YQUALSIC keywords.
 *		11/25/03   1.14 wvd   Raise upper limits on LIF_CNT_RATE and
 *					SIC_CNT_RATE to 32000.  Raise limits
 *					on FEC_CNT_RATE and AIC_CNT_RATE
 *					to 64000.  Use TZERO = 32768 for
 *					FEC_CNT_RATE and AIC_CNT_RATE. 
 *					They must be read as floats.
 *              01/15/04   1.15 bjg   Now performs cnt0 -= 16777216 when
 *                                      cnt0 < cnt before computing 
 *                                      cnt_rate for the first cnt and 
 *                                      cnt0 in function 
 *                                      fill_cntrate_array
 *                                    cf_timeline 1.6 -> 1.7
 *                                      free the arrays only at the end
 *                                      version 1.6 tries to read tflag
 *                                      array after it has been freed
 *		02/06/04   1.16 wvd   Modify fill_cntrate_array() so
 *					that the output array no longer runs
 *					16 seconds behind the HSKP data.
 *		02/18/04   1.17 wvd   Voltages for detectors 2A and 2B
 *					are still swapped.  Set HKERR_ENDS
 *					to 20100000.
 *              03/04/04   1.18 bjg   Read begin and end counters as float
 *                                    Set countrate to zero when cntb is more
 *                                    than 16777216.
 *              03/19/04        bjg   FITS data types (TFLOAT, TSTRING, TBYTE)
 *                                    now match type of C containers (float, 
 *                                    char *, char) when reading and writing 
 *                                    in FITS header.
 *              04/06/04   1.19 bjg   Include ctype.h
 *                                    Functions now return EXIT_SUCCESS
 *                                    Remove unused variables
 *                                    Remove static keyword in struct key 
 *                                    definition
 *                                    If cnt_rate>32000 set cnt_rate=0
 *		05/07/04   1.20 wvd   Delete HKERR_ENDS.  If HSKPVERS keyword
 *				      is not present in HSKP file, test to
 *				      determine whether high-voltage arrays
 *				      for 2A and 2B are swapped.
 *		05/27/04   1.21 bjg   Test counter keywords. If bad, use
 *				      cnt_rate=NEVENTS/exptime or 0
 *		06/14/04   1.22 wvd   Add BRIT_OBJ keyword to header.
 *				      Rewrite fill_hv_array() to make it more
 *					robust.
 *				      If HV values in header are good, use
 *					them to determine whether HV 
 *                                    	arrays for 2A and 2B are swapped.
 *              10/17/04        bjg   Replaced lots of
 *                                    "while (val<array[i] && i<nmax)"
 *                                    with
 *                                    "while (i<nmax && val<array[i])"
 *              11/09/04   1.23 bjg   Fixed Timeline AIC countrate. 
 *                                    Sum of AIC countrates 
 *                                    for all 4 detectors
 *		02/09/05   1.24 wvd   Use count rates, rather than counts,
 *					when testing header keywords.
 *		02/15/05   1.25 wvd   If the count-rate arrays get too big,
 *					truncate and issue a warning.
 *		02/28/05   1.26 wvd   Fix bug in AIC count-rate test.
 *				      Lower AIC and FEC thresholds to
 *					0.8 * NEVENTS/EXPTIME.
 *				      If GTI values span more than 25 ks,
 *					make timeline table entries only for
 *					times in the GTI's.
 *				      Set OPUS flag for timeline-table 
 *					entries whose values exceed 25 ks.
 *				      Rewrite fill_hv_array() to make it more
 *					robust.
 *				      Set overflow count-rate values to zero.
 *				      Treat HV as array of shorts throughout.
 *		03/03/05   1.27 wvd   Raise exposure-length limit to 55 ksec.
 *		03/14/05   1.28 wvd   Modify gtis[0] only if TTPERIOD = 1.0
 *		03/17/05   1.29 wvd   Issue a warning if bigtime = TRUE.
 *		03/21/05   1.30 wvd   If there's a break in the timeline
 *					array, recalculate times of sunrise
 *					and sunset.  Delete bigtime and the
 *					mechanism for storing tsunrise and
 *					tsunset as floats.  Define lastsunrise
 *					and lastsunset relative to ptime[0],
 *					rather than mjd_start.  Write time
 *					between sunsets to header in ORBPERID.
 *		08/30/05   1.31 wvd   Allow for the possibility that some of
 *					the photons have absurd arrival times
 *					AND others are not included in the
 *					OPUS good-time intervals.
 *		08/31/05   1.32 wvd   In cf_set_background_limits, use LWRS
 *					limits for RFPT observations.
 *		11/03/05   1.33 wvd   When housekeeping file is missing and
 *					engineering snapshot times are
 *					corrupted, use EXPSTART and EXPEND
 *					to estimate eng_time. If EXPEND
 *					is corrupted, use EXPTIME.
 *		04/09/06   1.34 wvd   Flag as bad times after final OPUS GTI.
 *		08/21/06   1.35 wvd   Treat housekeeping time array as double.
 *					Don't add 0.5 to arrays in subroutine
 *					fill_cntrate_array.
 *		12/29/06   1.36 wvd   Correct bug in construction of AIC
 *					count-rate array.
 *		02/13/07   1.37 wvd   If TTPERIOD = 0, set to 1.
 *		03/07/07   1.38 wvd   Remove pos from call to space_vel.
 *		03/23/07   1.39 wvd   Give more detailed error message if
 *					housekeeping file is not usable.
 *		04/07/07   1.40 wvd   Clean up compiler warnings.
 *		05/18/07   1.41 wvd   If HV array in HSKP file consists of a 
 *					single 0 followed by all -1's, use 
 *					header info to populate timeline table.
 *					Rename fill_hv_array to hv_from_hskp
 *					and add hv_from_header.
 *
 ****************************************************************************/
 
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "calfuse.h" 
#include "sgp4.h"


struct key {
        char keyword[FLEN_KEYWORD];
        float value;
        };


/***********************************************************************
 *
 *  procedure to get the GTI values from the input file
 *
 **********************************************************************/

int cf_get_gti(fitsfile *infits, double **gtis, double **gtie) {

    char	instmode[FLEN_VALUE];
    int		ngti, status=0;
    float	ttperiod;

    cf_verbose(3, "Entering cf_get_gti.");

   /* Read header keywords. */
    FITS_movabs_hdu(infits, 1, NULL, &status);
    FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(infits, TFLOAT, "TTPERIOD", &ttperiod, NULL, &status);
    if (ttperiod < 1E-6) ttperiod = 1.0;

   /* Move to the second extension and read the GTI's */
    FITS_movabs_hdu(infits, 3, NULL, &status);
    ngti = cf_read_col(infits, TDOUBLE, "START", (void **) gtis);
    ngti = cf_read_col(infits, TDOUBLE, "STOP",  (void **) gtie);

   /* 
    *  If gtis[0] = 0.000, then it is probably wrong.
    *  Set it equal to the fractional part of gtie[0].
    *  In HIST mode, gtis[0] is always 0., so don't change it.
    *  This trick only works if TTPERIOD = 1.0.
    */
    if (**gtis < FRAME_TOLERANCE && !strncmp(instmode, "TTAG", 4)
	&& fabs(ttperiod - 1.0) < FRAME_TOLERANCE)
	**gtis = fmod(*gtie[0], 1.);

    cf_verbose(3, "Exiting cf_get_gti.");
    return ngti;
}

/***********************************************************************
 *
 * procedure to get a list of geocoronal line locations from the AIRG 
 *   calibration file
 *
 ***********************************************************************/

int
cf_get_geocorona(fitsfile *outfits, short **gxmin, short **gxmax,
	short **gymin, short **gymax)
{
   char airgfile[FLEN_VALUE]; 
   int nrow, status=0; 
   fitsfile *airgfits;

   cf_verbose(3, "Entering cf_get_geocorona.");

  /* Read the name of the AIRG calibration file to use and open it */
   FITS_movabs_hdu(outfits, 1, NULL, &status);
   FITS_read_key(outfits, TSTRING, "AIRG_CAL", airgfile, NULL, &status);
   FITS_open_file(&airgfits, cf_cal_file(airgfile), READONLY, &status); 
   cf_verbose(3, "AIRG_CAL = %s", airgfile);

  /* 
   *  Move to the first extension and read the positions of the
   *  geocoronal lines.
   */  
   FITS_movabs_hdu(airgfits, 2, NULL, &status);
   nrow = cf_read_col(airgfits, TSHORT, "XMIN", (void **) gxmin); 
   nrow = cf_read_col(airgfits, TSHORT, "XMAX", (void **) gxmax);
   nrow = cf_read_col(airgfits, TSHORT, "YMIN", (void **) gymin);
   nrow = cf_read_col(airgfits, TSHORT, "YMAX", (void **) gymax);

   FITS_close_file(airgfits, &status);
   cf_verbose(3, "Exiting cf_get_geocorona.");
   return nrow;
}


/***********************************************************************
 *
 * Procedures to compute times of most recent sunrise and sunset
 *
 ***********************************************************************/

SGP4   set_orbit_parms(fitsfile *);

static double
sunrise(double mjd_start, SGP4 sgp4)
{
    double ang_sep, mjd, ptime, pos[3], vel[3];
    int    i, dn_flag, dn_flag_last;

    dn_flag_last = 1; /* Assume that we begin in night. */

    for (i=1; i>-8000; i--) {
        ptime = (double) i;
        mjd = mjd_start + ptime/86400.0;

        /* Get the state vector at time mjd */
        SGP4_getStateVector(sgp4, mjd, pos, vel);
        SGP4_precess(pos, mjd, MJD2000);

        dn_flag = eclipse(pos, mjd, &ang_sep); /* 0=day, 1=night */

	/* We're going backward in time, so sunrise is day into night. */
        if ((dn_flag == 1) && (dn_flag_last == 0)) break;

        dn_flag_last=dn_flag;
    }
    /* Historically, we've defined sunrise as the time increment after
	the change, so add 1 to the number we've just calculated. */
    return (ptime + 1.);
}


static double
sunset(double mjd_start, SGP4 sgp4)
{
    double ang_sep, mjd, ptime, pos[3], vel[3];
    int    i, dn_flag, dn_flag_last;

    dn_flag_last = 0; /* Assume that we begin in day. */

    for (i=1; i>-8000; i--) {
	ptime = (double) i;
	mjd = mjd_start + ptime/86400.0;

	/* Get the state vector at time mjd */
	SGP4_getStateVector(sgp4, mjd, pos, vel);
	SGP4_precess(pos, mjd, MJD2000);

	dn_flag = eclipse(pos, mjd, &ang_sep);  /* 0=day, 1=night */

	/* We're going backward in time, so sunset is night into day. */
	if ((dn_flag == 0) && (dn_flag_last == 1)) break;

	dn_flag_last = dn_flag;
    }
    /* Historically, we've defined sunset as the time increment after
	the change, so add 1 to the number we've just calculated. */
    return (ptime + 1.);
}

/****************************************************************************/

static void
fill_cntrate_array(long nsam_hk, long *hk_colval, double *time_hk,
		   long ntime, double *ttime, float *tarray)
{
    /*
     *  Procedure to fill an array with count rate housekeeping (HK)
     *  information, which is tabulated on an approx 16s time interval.
     *
     *  nsam_hk   = number of time samples in the housekeeping file
     *  hk_colval = array containing the housekeeping data
     *  time_hk   = time from the exposure start for each second in the
     *              houskeeping file
     *  ntime     = number of tabulated rows in the timeline array
     *  ttime     = tabulated times in the timeline array
     *  tarray    = timeline array to be filled, one sample per second
     *              starting at the start of the exposure
     */

    double hk_time, hk_time0;
    long hk_cnt, hk_cnt0;
    long i;	/* index through timeline table */
    long k;	/* index through housekeeping array */
    float cnt_rate = 0.;

    i = k = 0;	/* Initialize indices */

    /*
     *  Find the first non-negative value in the HK array
     */
    while(k < nsam_hk && hk_colval[k] < 0) k++;
    hk_cnt0  = hk_colval[k];
    hk_time0 = time_hk[k];

    /*
     *  Begin big loop through the timeline table
     */
    while (i < ntime) {

	/* Find the next non-negative value in the HK array. */
	k++;
	while(k < nsam_hk && hk_colval[k] < 0) k++;

	/* If we've run out of HK entries,
		fill in the rest of the timeline and quit. */
	if (k == nsam_hk) {
	    for ( ; i < ntime; i++)
		tarray[i] = cnt_rate;
	    break;
	}

	/* Otherwise, calculate the count rate between
		times hk_time0 and hk_time. */
	hk_cnt = hk_colval[k];
	hk_time = time_hk[k];
	if (hk_cnt < hk_cnt0)
	    hk_cnt0 -= 16777216;
	cnt_rate = (hk_cnt - hk_cnt0) / (hk_time - hk_time0);

	/* Populate the timeline table for times less than hk_time */
	while (i < ntime && cf_nlong(ttime[i]) < hk_time) 
	    tarray[i++] = cnt_rate;

	/* Save latest HK time and count values. */
	hk_cnt0  = hk_cnt;
	hk_time0 = hk_time;
    }
}


/****************************************************************************/

static void
hv_from_header(fitsfile *outfits, char *det, char *side, long ntime, short *hv)
{

	char  hv_key[FLEN_CARD];
	int   det_hv, hv_flag, status=0;
	long  i;

        fits_read_key(outfits, TINT, "HV_FLAG", &hv_flag, NULL, &status);
	if (hv_flag == -1) {
           det_hv = -999;
           cf_verbose(1, "Bad HV keywords: setting HV to -999 in timeline table.");
	}
	else {
           sprintf(hv_key, "DET%.1sHV%.1sH", det, side);
           FITS_read_key(outfits, TINT, hv_key, &det_hv, NULL, &status);
	}
	if (hv_flag == 0)
	   cf_verbose(1, "Applying max voltage value to whole exposure.");

        cf_verbose(2, "High voltage = %d", det_hv);
	for (i=0; i<ntime; i++)
	    hv[i] = det_hv;
}


/****************************************************************************/

int
hv_from_hskp(long nsam_hk, long *hk_colval, double *time_hk,
	      long ntime, double *ttime, short *tarray)
{
    /*
     *  Procedure to fill an array with high voltage housekeeping (HK)
     *  information, which is tabulated on an approx 16s time interval.
     *
     *  nsam_hk   = number of time samples in the housekeeping file
     *  hk_colval = array containing the housekeeping data
     *  time_hk   = time from the exposure start for each second in the
     *              housekeeping file
     *  ntime     = number of tabulated rows in the timeline array
     *  ttime     = tabulated times in the timeline array
     *  tarray    = timeline array to be filled, one sample per second
     *              starting at the start of the exposure
     */

    long i;     /* index through timeline table */
    long k;     /* index through housekeeping array */
    short hv = 0; 
    
    i = k = 0;  /* Initialize indices */

    /* 
     *  Scan the input HK array.  If it consists of a zero followed by
     *  all -1's, exit with an error.
     */

    if (hk_colval[0] == 0) {
	for (i = 1; i < nsam_hk; i++) if (hk_colval[i] != -1) break;
	if (i == nsam_hk) return (-1);
    }

    /* 
     *  Find the first non-negative value in the HK array
     */

    while(k < nsam_hk && hk_colval[k] < 0) k++;
    hv = hk_colval[k];

    /*
     *  Begin big loop through the timeline table
     */

    while (i < ntime) {

	/* Find the next non-negative value in the HK array */
	k++;
	while(k < nsam_hk && hk_colval[k] < 0) k++;

	/* If we've run out of HK entries, fill in the rest of the
	   timeline and quit. */
	if (k == nsam_hk) {
	    for ( ; i < ntime; i++)
		tarray[i] = hv;
	    break;
	}

	/* Otherwise, populate the timeline table for times less than hk_time */
	while (i < ntime && ttime[i] < time_hk[k])
	    tarray[i++] = hv;

	/* Save latest HV value from HK array */
	hv = hk_colval[k];
    }

    return (0);
}


/****************************************************************************/

int cf_timeline(fitsfile *outfits)
{
  /*************************************************************
   *
   *   Procedure to fill in the third extension of the IDF,
   *   which contains information about orbital parameters,
   *   count rates, high voltage info, and Y centroids as
   *   a function of time during the observation. The data
   *   are tabulated once each second.
   *
   **************************************************************/

    char  hkexists[FLEN_VALUE], det[FLEN_VALUE], side[2];
    char  fmt_byte[FLEN_CARD], fmt_float[FLEN_CARD], fmt_short[FLEN_CARD];
    char  instmode[FLEN_VALUE], keyword[FLEN_KEYWORD], card[FLEN_CARD];
    char  volt_cal[FLEN_VALUE];
    unsigned char *tflag;
    int   hv_flag, ngti, n_bad_times, TIMES_TOO_BIG = FALSE;
    int   status=0, anynull, intnull=0, dn_flag_last=0;
    int   biglif = FALSE, bigsic = FALSE, bigfec = FALSE, bigaic = FALSE; 
    long  i, j, k, expnight, nevents, ntime, vmax;
    short *hv;
    float *lifcr, *siccr, *feccr, *aiccr, *bkgdcr;
    float first_time, last_time, bad_times[1024], *time=NULL;
    float *ycentl, *ycents;
    float max_lif, max_sic, max_fec, max_aic;
    double *gtis, *gtie, min_limb=999., period;
    double mjd_start, mjd_end, mjd, *ptime, *lon, *lat, *limbang;
    double *vorb, *tsunrise, lastsunrise=0, *tsunset, lastsunset=0;
    orbital orb;
    SGP4  sgp4;
    struct key tscal[16], tzero[16];
    long nevts;

    char * dets[]={"1A","1B","2A","2B"};

    char extname[]="TIMELINE";
    int tfields=16;  /* output table will have 16 columns */
    char *ttype[]={"TIME", "STATUS_FLAGS", "TIME_SUNRISE", "TIME_SUNSET",
		   "LIMB_ANGLE", "LONGITUDE", "LATITUDE", "ORBITAL_VEL",
		   "HIGH_VOLTAGE", "LIF_CNT_RATE", "SIC_CNT_RATE",
		   "FEC_CNT_RATE", "AIC_CNT_RATE", "BKGD_CNT_RATE",
                   "YCENT_LIF","YCENT_SIC"};

    char *tform[16]; /* we will define tform later, when the number
			of photons is known */

    char *tunit[]={"seconds", "unitless", "seconds", "seconds", "degrees",
                   "degrees", "degrees",  "km/s", "unitless", "counts/sec",
                   "counts/sec", "counts/sec", "counts/sec", "counts/sec",
                   "pixels","pixels" };

    char CF_PRGM_ID[] = "cf_timeline";
    char CF_VER_NUM[] = "1.41";

    /* Initialize error checking. */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");

    /*  Set the orbital parameters keywords. */
    FITS_movabs_hdu(outfits, 1, NULL, &status);
    FITS_read_key(outfits, TDOUBLE, "EXPSTART", &mjd_start, NULL, &status);
    FITS_read_key(outfits, TDOUBLE, "EXPEND",   &mjd_end, NULL, &status);
    FITS_read_key(outfits, TDOUBLE, "RA_TARG",  &orb.ra_ap, NULL, &status);
    FITS_read_key(outfits, TDOUBLE, "DEC_TARG", &orb.dec_ap, NULL, &status);
    FITS_read_key(outfits, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(outfits, TSTRING, "VOLT_CAL", volt_cal, NULL, &status);

    cf_velang(outfits, (mjd_end+mjd_start)/2.0);

    sgp4 = set_orbit_parms(outfits);

    /* Read the good time intervals */
    ngti = cf_get_gti(outfits, &gtis, &gtie);
    for (i=0; i<ngti; i++)
	cf_verbose(3, "GTI start = %0.1f, end= %0.1f", gtis[i], gtie[i]);

    /* 
     * An exposure should never span more than 55 ks.  If it appears to
     * do so, include only the good-time intervals in the timeline table.
     * Otherwise, include all times between EXPSTART and EXPEND, whether
     * or not there are any associated photons.
     */
    if (gtie[ngti-1] - gtis[0] > MAX_EXPTIME) {
	TIMES_TOO_BIG = TRUE;
	cf_if_warning("Exposure appears to span %0.0f ks.  Truncating "
		"timeline table.", (gtie[ngti-1] - gtis[0])/1000.);
    }

    /* For HIST data, use the good-time intervals
	to set the times in the timeline table. */ 
    if (!strncmp(instmode, "HIST", 4)) {
        if (TIMES_TOO_BIG) {
	    ntime = 0;
	    for (i = 0; i < ngti; i++) 
		ntime += (long) (gtie[i] - gtis[i] + 1.5);
	    ptime = (double *) cf_malloc(sizeof(double) * ntime);
	    k = 0;
	    for (i = 0; i < ngti; i++) {
		ntime = (long) (gtie[i] - gtis[i] + 1.5);
		for (j = 0; j < ntime; j++)
		    ptime[k++] = gtis[i] + j;
	    }
	    ntime = k;
	}
	else {
	    ntime = (long) (gtie[ngti-1] - gtis[0] + 1.5);
	    ptime = (double *) cf_malloc(sizeof(double) * ntime);
	    for (i = 0; i < ntime; i++)
		ptime[i] = gtis[0] + i;
	}
    }

    /* For TTAG data, use the photon times themselves. */
    else {
	FITS_movabs_hdu(outfits, 2, NULL, &status);
	nevents = cf_read_col(outfits, TFLOAT, "TIME", (void **) &time);
	FITS_movabs_hdu(outfits, 1, NULL, &status);
	n_bad_times = 0;
        if (TIMES_TOO_BIG) {
	    last_time = 0;
	    i = nevents - 1;
	    while (time[i] > MAX_EXPTIME && i > 0) {
		if (fabs(time[i] - last_time) > FRAME_TOLERANCE && n_bad_times < 1024) {
		    bad_times[n_bad_times++] = time[i];
		    last_time = time[i];
		}
		i--;
	    }
	    last_time = time[i] + 1.;
	}
	else
	    last_time = time[nevents-1] + 1.;
	first_time = time[0];
	if (first_time < fmod(last_time, 1.))
		first_time = fmod(last_time, 1.) - 1.;
	ntime = (long) (last_time - first_time + 1.5) + n_bad_times;
	ptime = (double *) cf_malloc(sizeof(double) * ntime);
	for (i = 0; i < ntime - n_bad_times; i++)
	    ptime[i] = first_time + i;
	for (j = n_bad_times-1; i < ntime; i++, j--)
	    ptime[i] = bad_times[j];
	free(time);
    }
    cf_verbose(2, "Timeline table will contain %ld entries", ntime);

    /* Now generate all of the other arrays in the timeline table. */

    lon   = (double *) cf_calloc(ntime, sizeof(double));
    lat   = (double *) cf_calloc(ntime, sizeof(double));
    limbang  = (double *) cf_calloc(ntime, sizeof(double));
    vorb  = (double *) cf_calloc(ntime, sizeof(double));
    tsunrise = (double *) cf_calloc(ntime, sizeof(double));
    tsunset = (double *) cf_calloc(ntime, sizeof(double));
    hv    = (short *) cf_calloc(ntime, sizeof(short));
    lifcr = (float *) cf_calloc(ntime, sizeof(float));
    siccr = (float *) cf_calloc(ntime, sizeof(float));
    feccr = (float *) cf_calloc(ntime, sizeof(float));
    aiccr = (float *) cf_calloc(ntime, sizeof(float));  
    bkgdcr = (float *) cf_calloc(ntime, sizeof(float));
    tflag = (unsigned char *) cf_calloc(ntime, sizeof(unsigned char));
    ycentl = (float *)cf_calloc(ntime, sizeof(float));
    ycents = (float *)cf_calloc(ntime, sizeof(float));

    for (i=0; i<ntime; i++) {
	int   dn_flag, day_limb;
	float dx;
	double geo_lon, geo_lat, gmst, pos[3], vel[3], zdist, ang_sep;

	mjd = mjd_start + ptime[i]/86400.0;

	/* get the state vector at time mjd */
	SGP4_getStateVector(sgp4, mjd, pos, vel);
	SGP4_precess(pos, mjd, MJD2000);
	SGP4_precess(vel, mjd, MJD2000);

	dx = space_vel(vel, orb.ra_ap, orb.dec_ap);

	state_geod(pos, mjd, &geo_lon, &geo_lat, &gmst);
	limbang[i] = (double) state_limb(pos, mjd, orb.ra_ap, orb.dec_ap, 
					 &zdist, &day_limb);
	if (min_limb > limbang[i]) min_limb = limbang[i];

	dn_flag = eclipse(pos, mjd, &ang_sep);  /* 0=day, 1=night */

	/* The first time through, and after any breaks in the timeline,
	 * compute times of most recent sunrise and sunset. */
	if ((i == 0) || (ptime[i] - ptime[i-1] - 1. > FRAME_TOLERANCE)) {
	    lastsunrise = ptime[i] + sunrise(mjd, sgp4);
	    lastsunset  = ptime[i] + sunset(mjd, sgp4);
	    if (lastsunset > lastsunrise) dn_flag_last = 1 ;
	    else dn_flag_last = 0 ;
            /*  The first time through, compute the orbital period and
	     *	write it to file header.  */
	    if (i == 0) {
		mjd += (lastsunset - ptime[i] + 7000.)/86400.0;
		period = 7000. + sunset(mjd, sgp4);
		FITS_movabs_hdu(outfits, 1, NULL, &status);
		FITS_update_key(outfits, TDOUBLE, "ORBPERID", &period,
		"[s] Estimate of orbital period", &status);
	    }
	}
	/* If we've moved from day into night (or night into day), update
	   last sunrise/sunset times. */
	else if ((dn_flag == 0) && (dn_flag_last == 1))
	    lastsunrise = ptime[i];
	else if ((dn_flag == 1) && (dn_flag_last == 0)) 
      	    lastsunset = ptime[i]; 

	/* If this is daytime, set the day bit in TFLAG array. */
	if (dn_flag == 0)
	    tflag[i] |= TEMPORAL_DAY;

	lon[i] = geo_lon;
	lat[i] = geo_lat;
	vorb[i] = dx;
	tsunrise[i] = ptime[i] - lastsunrise;
	tsunset[i]  = ptime[i] - lastsunset;
	dn_flag_last = dn_flag;
    }


    /* 
     * Set the OPUS flag for times outside of the OPUS good-time
     * intervals.  Since there are no photons associated with the
     * last second of a GTI, we make sure that it is flagged as
     * bad.  The first second of a good-time interval is good.
     */
    j = 0;
    for (i = 0; i < ntime; i++) {
	while (j < ngti-1 && ptime[i] > gtie[j] - 0.5)
	    j++;
	if (ptime[i] < gtis[j] - FRAME_TOLERANCE ||
	    ptime[i] > gtie[j] + FRAME_TOLERANCE)
	    tflag[i] |= TEMPORAL_OPUS;
    }
    tflag[ntime-1] |= TEMPORAL_OPUS;

    /*
     * No exposure should be longer than 55 ks.  For HIST data, we'll just
     * put up with it, but for TTAG data, we'll flag as bad any timeline
     * table entries with absurd arrival times.  We use the OPUS flag for
     * this.  What we're really saying is that we don't believe that the 
     * times associated with these photons are correct.
     */
    for (i = 0; i < ntime; i++)
	if (ptime[i] > MAX_EXPTIME)
	    tflag[i] |= TEMPORAL_OPUS;

    /* Fill housekeeping columns (count rates and high voltage). */
    FITS_movabs_hdu(outfits, 1, NULL, &status);
    FITS_read_key(outfits, TSTRING, "HKEXISTS", hkexists, NULL, &status);
    FITS_read_key(outfits, TSTRING, "DETECTOR", det, NULL, &status);
    strncpy(side, det+1, 1);

    /* Populate MIN_LIMB keyword. */
    FITS_update_key(outfits, TDOUBLE, "MIN_LIMB", &min_limb, NULL, &status);

    /* If no housekeeping file exists, use header info to fill timeline. */
 loop:
    if(!strncasecmp(hkexists, "N", 1)) {
	char  cntb_key[FLEN_CARD], cnte_key[FLEN_CARD];
	float  cntb, cnte;
	float cnt_rate, cnt_rate2, eng_time, exp_time;
	double ctimeb, ctimee;

        cf_verbose(1, "No housekeeping file: filling timeline with header info.");

	/* move to the top level header to read the relevant keywords */ 
        FITS_movabs_hdu(outfits, 1, NULL, &status);
	FITS_read_key(outfits, TFLOAT, "EXPTIME", &exp_time, NULL, &status);
	FITS_read_key(outfits, TLONG, "NEVENTS", &nevts, NULL, &status);

	/* calculate engineering count rates */ 
	FITS_read_key(outfits, TDOUBLE, "CTIME_B", &ctimeb, NULL, &status);
	FITS_read_key(outfits, TDOUBLE, "CTIME_E", &ctimee, NULL, &status);
	eng_time = (ctimee-ctimeb) * 86400.;

        /* If eng_time is unreasonable, use the EXPTIME keyword */
	if (eng_time < 1) {
	  cf_if_warning("Engineering snapshot time less than or equal to zero.");
	  if ((eng_time = (mjd_end-mjd_start) * 86400.) < MAX_EXPTIME)
	     cf_if_warning("Estimating LiF, SiC, FEC and AIC count rates from EXPSTART and EXPEND.");
	  else {
	     eng_time = exp_time;
	     cf_if_warning("Estimating LiF, SiC, FEC and AIC count rates from EXPTIME.");
	  }
	}

	/* SiC count rate timeline */
        sprintf(cntb_key, "C%.2sSIC_B", det);
        FITS_read_key(outfits, TFLOAT, cntb_key, &cntb, NULL, &status);
        sprintf(cnte_key, "C%.2sSIC_E", det);
        FITS_read_key(outfits, TFLOAT, cnte_key, &cnte, NULL, &status); 
				
	if ((cnte<0) || (cntb<0) || (cntb>cnte) || eng_time < 1. ||
	    (cnte==0xEB90EB90) || (cnte==0xDEADDEAD) ||
	    (cntb==0xEB90EB90) || (cntb==0xDEADDEAD) ||
	    ((cnt_rate = (cnte - cntb) / eng_time) > 32000)) {
	    cf_if_warning("Bad SiC counter.  SiC count rate will be set to zero.");
	    cnt_rate=0;
	}
	
	cf_verbose(2, "SiC count rate = %d", (int) cnt_rate);
	for (i=0; i<ntime; i++)
	    siccr[i] = cnt_rate;

	/* LiF count rate timeline */
        sprintf(cntb_key, "C%.2sLIF_B", det);
        FITS_read_key(outfits, TFLOAT, cntb_key, &cntb, NULL, &status);
        sprintf(cnte_key, "C%.2sLIF_E", det);
        FITS_read_key(outfits, TFLOAT, cnte_key, &cnte, NULL, &status);
      
        if ((cnte<0) || (cntb<0) || (cntb>cnte) || eng_time < 1. ||
            (cnte==0xEB90EB90) || (cnte==0xDEADDEAD) ||
            (cntb==0xEB90EB90) || (cntb==0xDEADDEAD) ||
            ((cnt_rate = (cnte - cntb) / eng_time) > 32000)) {
            cf_if_warning("Bad LiF counter.  LiF count rate will be set to zero.");
            cnt_rate=0;
        }

	cf_verbose(2, "LiF count rate = %d", (int) cnt_rate);
	for (i=0; i<ntime; i++)
	    lifcr[i] = cnt_rate;

	/* FEC count rate timeline */
        sprintf(cntb_key, "C%.2sFE_B", det);
        FITS_read_key(outfits, TFLOAT, cntb_key, &cntb, NULL, &status);
        sprintf(cnte_key, "C%.2sFE_E", det);
        FITS_read_key(outfits, TFLOAT, cnte_key, &cnte, NULL, &status);
      
	if ((cnte<0) || (cntb<0) || (cntb>cnte) || eng_time<1 ||
	    (cnte==0xEB90EB90) || (cnte==0xDEADDEAD) ||
	    (cntb==0xEB90EB90) || (cntb==0xDEADDEAD) ||
	    ((cnt_rate = (cnte - cntb) / eng_time) > 64000) ||
	    (cnt_rate < 0.8*nevts/exp_time)) {
	    cf_if_warning("Bad FEC counter"); 
	    cf_if_warning("-- FEC count rate will be computed from NEVENTS and EXPTIME.");
	    cf_if_warning("-- Electronic deadtime correction will be underestimated.");
	    cf_if_warning("-- Y stretch will be underestimated.");
	    cnt_rate=nevts/exp_time;
	}			
				
	cf_verbose(2, "FEC count rate = %d", (int) cnt_rate);
	for (i=0; i<ntime; i++)
	    feccr[i] = cnt_rate;

	/* AIC count rate timeline */
	cnt_rate=0;
	for (i=0;i<4;i++){
	  sprintf(cntb_key, "C%.2sAI_B", dets[i]);
	  FITS_read_key(outfits, TFLOAT, cntb_key, &cntb, NULL, &status);
	  sprintf(cnte_key, "C%.2sAI_E", dets[i]);
	  FITS_read_key(outfits, TFLOAT, cnte_key, &cnte, NULL, &status);			

	  if (((cnte<0) || (cntb<0) || (cntb>cnte) || (eng_time<1) ||
	    (cnte==0xEB90EB90) || (cnte==0xDEADDEAD) ||
	    (cntb==0xEB90EB90) || (cntb==0xDEADDEAD) ||
	    ((cnt_rate2 = (cnte - cntb)/ eng_time) > 64000)) ||
	    ((!(strncasecmp(dets[i],det,2))) && (cnt_rate2<0.8*nevts/exp_time))) {
	      cf_if_warning("Bad AIC counter %s",dets[i]); 
	      cf_if_warning("-- AIC count rate will be computed from NEVENTS and EXPTIME.");
	      cf_if_warning("-- IDS deadtime correction will be underestimated.");
	      cnt_rate=nevts/exp_time;
	      break;
	  }	
	  cnt_rate+=cnt_rate2;
	}
       
	cf_verbose(2, "AIC count rate = %d", (int) cnt_rate);
	for (i=0; i<ntime; i++)
	    aiccr[i] = cnt_rate;

	/* High voltage timeline */
	hv_from_header(outfits, det, side, ntime, hv);
    }
    else {
	/*
	 *  Housekeeping file exists.  Use it to calculate timeline values.
	 */
	char  cntb_key[FLEN_VALUE], hk_file_name[FLEN_VALUE];
	int   hk_colnum, hskpvers, MUST_TEST_HV=FALSE, swap_HV=FALSE;
	long  nsam_hk, *hk_colval;
	double *hk_mjd, *time_hk;
	fitsfile *hskpfits;
	float * aiccr2;

	FITS_read_key(outfits, TSTRING, "HSKP_CAL", hk_file_name, NULL, &status);
	/* If you can't open the housekeeping file, try a lower-case filename.*/
	if (fits_open_file(&hskpfits, hk_file_name, READONLY, &status)) {
            status = 0;
            cf_verbose(3, "Can't find housekeeping file %s", hk_file_name);
            hk_file_name[0] = (char) tolower(hk_file_name[0]);
            cf_verbose(3, "Will look for file named %s", hk_file_name);

	    /* If you still can't open the housekeeping file, use header keywords. */
            if (fits_open_file(&hskpfits, hk_file_name, READONLY, &status)) {
		status = 0;
                cf_verbose(3, "Can't find housekeeping file %s", hk_file_name);
		strncpy(hkexists, "N", 1);
		goto loop;
	    }
	    /* If lower-case filename works, write it to the file header. */
            FITS_update_key(outfits, TSTRING, "HSKP_CAL", hk_file_name, NULL, &status);
        }
	cf_verbose(1, "Housekeeping file exists: using it to fill timeline.");

	/* If we're on side 2 and the HSKPVERS keyword does not exist, then we'll
		need to check whether the HV arrays are swapped. */
        if (*det == '2' &&
	    fits_read_key(hskpfits, TINT, "HSKPVERS", &hskpvers, NULL, &status)) {
	    status = 0;
	    MUST_TEST_HV = TRUE;
	}

        FITS_movabs_hdu(hskpfits, 2, NULL, &status);
        FITS_read_key(hskpfits, TLONG, "NAXIS2", &nsam_hk, NULL, &status);
	cf_verbose(3, "Number of samples in the HK file = %ld", nsam_hk);

	/*  Allocate space to hold the column contents */
	hk_colval = (long *) cf_calloc(nsam_hk, sizeof(long));
        hk_mjd    = (double *) cf_calloc(nsam_hk, sizeof(double));
	time_hk   = (double *) cf_calloc(nsam_hk, sizeof(double));

	/* Read in the various columns */
        
	/* First the MJD values - convert to time from start in seconds  */
	FITS_get_colnum(hskpfits, TRUE, "MJD", &hk_colnum, &status);
        FITS_read_col(hskpfits, TDOUBLE, hk_colnum, 1L, 1L, nsam_hk, &intnull,
		      hk_mjd, &anynull, &status);
        for (i=0; i<nsam_hk; i++)
	    time_hk[i] = (hk_mjd[i]-mjd_start) * 86400.;

	/* SiC count rate timeline */
        sprintf(cntb_key, "I_DET%.1sCSIC%.1s", det, side);

	/* Check to see whether the column exists. If not, then 
           reset the housekeeping flag and do the analysis using 
           keyword values in the input file header. If so, then check
           whether it contains real values. If not, then use the 
           keywords */
	if(!fits_get_colnum(hskpfits, TRUE, cntb_key, &hk_colnum, &status)) {
           FITS_read_col(hskpfits, TLONG, hk_colnum, 1L, 1L, nsam_hk, &intnull,
		      hk_colval, &anynull, &status);
           vmax = 0;
           for (j=0; j< nsam_hk; j++)
	       if (hk_colval[j] > vmax)
		   vmax = hk_colval[j];
	   if (vmax <= 0) {
              cf_if_warning("Data missing from housekeeping column.  "
		"Will treat file as missing.");
              strncpy(hkexists, "N", 1);
              status=0;
              goto loop; }
            else fill_cntrate_array(nsam_hk, hk_colval, time_hk, ntime,
				    ptime, siccr);}
	 else {
	      cf_if_warning("Data column missing from housekeeping file.  "
		"Will treat file as missing.");
              strncpy(hkexists, "N", 1);
              status=0;
              goto loop;
	    }
	cf_verbose(3,"SiC count-rate info read from housekeeping file.");

	/* LiF count rate timeline */
        sprintf(cntb_key, "I_DET%.1sCLIF%.1s", det, side);
	FITS_get_colnum(hskpfits, TRUE, cntb_key, &hk_colnum, &status);
        FITS_read_col(hskpfits, TLONG, hk_colnum, 1L, 1L, nsam_hk, &intnull,
		      hk_colval, &anynull, &status);
        fill_cntrate_array(nsam_hk, hk_colval, time_hk, ntime, ptime, lifcr);
	cf_verbose(3,"LiF count-rate info read from housekeeping file.");

	/* FEC count rate timeline */
	sprintf(cntb_key, "I_DET%.1sCFE%.1s", det, side);
	FITS_get_colnum(hskpfits, TRUE, cntb_key, &hk_colnum, &status);
	FITS_read_col(hskpfits, TLONG, hk_colnum, 1L, 1L, nsam_hk, &intnull,
		      hk_colval, &anynull, &status);
        fill_cntrate_array(nsam_hk, hk_colval, time_hk, ntime, ptime, feccr);
	cf_verbose(3,"FEC count-rate info read from housekeeping file.");

 	/* AIC count rate timeline */ 
	aiccr2 = (float *) cf_calloc(ntime, sizeof(float));
	for (j=0;j<ntime;j++) aiccr[j]=0;
	for (i=0;i<4;i++){
	  sprintf(cntb_key, "I_DET%.1sCAI%.1s", dets[i], dets[i]+1);
	  FITS_get_colnum(hskpfits, TRUE, cntb_key, &hk_colnum, &status);
	  FITS_read_col(hskpfits, TLONG, hk_colnum, 1L, 1L, nsam_hk, &intnull,
			hk_colval, &anynull, &status); 
	  fill_cntrate_array(nsam_hk, hk_colval, time_hk, ntime, ptime, aiccr2); 
	  for (j=0;j<ntime;j++) aiccr[j]+=aiccr2[j];
	}
	cf_verbose(3,"AIC count-rate info read from housekeeping file.");

	/* If MUST_TEST_HV = TRUE, determine whether HV for 2A and 2B are swapped. */
	if (MUST_TEST_HV) {
	    char  mjd_str[FLEN_VALUE], full_str[FLEN_VALUE], saa_str[FLEN_VALUE];
	    double mjd_volt;
	    int full, hdr_max2a, max2a = 0, saa;
	    long n, nhk, *hv2a;
	    fitsfile *voltfits;

	    cf_verbose(3, "Testing whether HV arrays are swapped.");

	    /* Read HV array for 2A from HSKP file; find max value. */
	    nhk = cf_read_col(hskpfits, TLONG, "I_DET2HVBIASAST", (void **) &hv2a);
	    for (n = 0L; n < nhk; n++) if (max2a < hv2a[n]) max2a = hv2a[n];
	    free (hv2a);
	    
	    /*
	     *  If HV values in file header are good, read 2A max.
	     *	Compare with value from HSKP file.
	     */
	    fits_read_key(outfits, TINT, "HV_FLAG", &hv_flag, NULL, &status);
	    if (hv_flag > -1) {
		FITS_read_key(outfits, TINT, "DET2HVAH", &hdr_max2a, NULL, &status);
		cf_verbose(3, "   Max for 2A: header says %d, housekeeping file says %d.",
		    hdr_max2a, max2a);
		if (abs(max2a - hdr_max2a) < 5)
		    cf_verbose(3, "   HV arrays are OK.  No need to swap.");
		else
		    swap_HV = TRUE;
	    }

	    /*
	     *  If HV values in file header are bad, compare with expected
	     *  values from VOLT_CAL file.
	     */
	    else {
		FITS_open_file(&voltfits, cf_cal_file(volt_cal), READONLY, &status);
	        n = 0L;
	        do {
	            n++;
	            sprintf(mjd_str, "MJD%ld", n);
	            FITS_read_key(voltfits, TDOUBLE, mjd_str, &mjd_volt, NULL, &status);
	        } while (mjd_start > mjd_volt);
	        n--;

	        sprintf(full_str, "FULL%ld", n);
	        sprintf(saa_str, "SAA%ld", n);
	        FITS_read_key(voltfits, TINT, full_str, &full, NULL, &status);
	        FITS_read_key(voltfits, TINT, saa_str, &saa, NULL, &status);
	        FITS_close_file(voltfits, &status);

	        /* If max for 2A matches expected full or SAA voltage -- and
		    we're on side 2A -- then we're OK.  Otherwise, must swap. */
	        cf_verbose(3, "   2A max = %d; Expected values for %s: full = %d, SAA = %d",
		    max2a, det, full, saa);
	        if ((abs(max2a - full) < 10 || abs(max2a - saa) < 5) &&
		    (*side == 'A'))
		    cf_verbose(3, "   HV arrays are OK.  No need to swap.");
	        else swap_HV = TRUE;
	    }
	}

 	/* High voltage timeline */
	if (swap_HV) {
	    if (*side == 'A') sprintf(cntb_key, "I_DET2HVBIASBST");
	    else              sprintf(cntb_key, "I_DET2HVBIASAST");
	    cf_verbose(3, "   Swapping HV values for detectors 2A and 2B");
	}
        else sprintf(cntb_key, "I_DET%.1sHVBIAS%.1sST", det, side);
	FITS_get_colnum(hskpfits, TRUE, cntb_key, &hk_colnum, &status);
        FITS_read_col(hskpfits, TLONG, hk_colnum, 1L, 1L, nsam_hk, &intnull,
		      hk_colval, &anynull, &status);
        if (!hv_from_hskp(nsam_hk, hk_colval, time_hk, ntime, ptime, hv))
	    cf_verbose(3,"HV info read from housekeeping file.");
        else {
	    cf_verbose(1,"Housekeeping file corrupted: reading HV info from file header.");
	    hv_from_header(outfits, det, side, ntime, hv);
	}

	free(time_hk);
	free(hk_mjd);
	free(hk_colval);
	FITS_close_file(hskpfits, &status);
    }

    /*
     * Set TFORM, TSCALE, and TZERO values for output table.
     */
    sprintf(fmt_byte,  "%ldB", ntime);
    sprintf(fmt_float, "%ldE", ntime);
    sprintf(fmt_short, "%ldI", ntime);

    /* First set default values. */
    tform[0] = fmt_float;
    tform[1] = fmt_byte;
    for (i=2; i<tfields; i++)
	tform[i] = fmt_short;

    /* For most arrays, we keep one decimal place. */
    for (i=2; i<tfields; i++) {
	sprintf(tscal[i].keyword, "TSCAL%ld", i+1);
	tscal[i].value = 0.1;
	sprintf(tzero[i].keyword, "TZERO%ld", i+1);
	tzero[i].value = 0.;
    }

    /* Now we start tinkering. */

    /* The orbital velocity is small, so we can keep two decimal places. */
    tscal[7].value = 0.01;

    /* FEC_CNT_RATE and AIC_CNT_RATE can get big, so we
        use TZERO AND TSCALE to compress them. */
    tscal[11].value = 1;
    tscal[12].value = 2;
    tzero[11].value = 32768;
    tzero[12].value = 65536;

    /* If a count rate is high, set to zero. */
    max_lif = max_sic = 32767.;
    max_fec = 65535.;
    max_aic = 131070.;

    for (i = 0; i < ntime; i++) {
	if (lifcr[i] > max_lif) {
		lifcr[i] = 0.;
		biglif = TRUE;
	}
	if (siccr[i] > max_sic) {
		siccr[i] = 0.;
		bigsic = TRUE;
	}
	if (feccr[i] > max_fec) {
		feccr[i] = 0.;
		bigfec = TRUE;
	}
	if (aiccr[i] > max_aic) {
		aiccr[i] = 0.;
		bigaic = TRUE;
	}
    }

	if (biglif) cf_if_warning("LIF_CNT_RATE out of bounds.  Setting bad values to zero.");
	if (bigsic) cf_if_warning("SIC_CNT_RATE out of bounds.  Setting bad values to zero.");
	if (bigfec) cf_if_warning("FEC_CNT_RATE out of bounds.  Setting bad values to zero.");
	if (bigaic) cf_if_warning("AIC_CNT_RATE out of bounds.  Setting bad values to zero.");

    /*
     * Write the timeline table to the output file.
     */
    FITS_movabs_hdu(outfits, 3, NULL, &status);
    FITS_create_hdu(outfits, &status);

    FITS_create_tbl(outfits, BINARY_TBL, 1L, tfields, ttype, 
                    tform, tunit, extname, &status);

    /* Write TSCALE and TZERO entries to header */
    for (i=2; i<tfields; i++) {

        /* Omit TSCALE and TZERO for these six arrays. */
        if (i ==  2) continue;  /* TIME_SUNRISE */
        if (i ==  3) continue;  /* TIME_SUNSET  */
        if (i ==  8) continue;  /* HIGH_VOLTAGE */
        if (i ==  9) continue;  /* LIF_CNT_RATE */
        if (i == 10) continue;  /* SIC_CNT_RATE */
        if (i == 13) continue;  /* BKGD_CNT_RATE */

	sprintf(keyword, "TUNIT%ld", i+1);
	if (fits_read_keyword(outfits, keyword, card, NULL, &status))
		cf_if_fits_error(status);
	FITS_insert_key_flt(outfits, tscal[i].keyword, tscal[i].value,
		-1, NULL, &status);
	FITS_insert_key_flt(outfits, tzero[i].keyword, tzero[i].value,
		-5, NULL, &status);
    }

    FITS_write_col(outfits, TDOUBLE, 1, 1L, 1L, ntime, ptime, &status);
	cf_verbose(3, "ptime array written to timeline table");
    FITS_write_col(outfits, TBYTE,   2, 1L, 1L, ntime, tflag, &status);
	cf_verbose(3, "tflag array written to timeline table");
    FITS_write_col(outfits, TDOUBLE, 3, 1L, 1L, ntime, tsunrise, &status);
	cf_verbose(3, "tsunrise array written to timeline table");
    FITS_write_col(outfits, TDOUBLE, 4, 1L, 1L, ntime, tsunset, &status);
	cf_verbose(3, "tsunset array written to timeline table");
    FITS_write_col(outfits, TDOUBLE, 5, 1L, 1L, ntime, limbang, &status);
	cf_verbose(3, "limbang array written to timeline table");
    FITS_write_col(outfits, TDOUBLE, 6, 1L, 1L, ntime, lon, &status);
	cf_verbose(3, "lon array written to timeline table");
    FITS_write_col(outfits, TDOUBLE, 7, 1L, 1L, ntime, lat, &status);
	cf_verbose(3, "lat array written to timeline table");
    FITS_write_col(outfits, TDOUBLE, 8, 1L, 1L, ntime, vorb, &status);
	cf_verbose(3, "vorb array written to timeline table");
    FITS_write_col(outfits, TSHORT,  9, 1L, 1L, ntime, hv, &status);
	cf_verbose(3, "hv array written to timeline table");
    FITS_write_col(outfits, TFLOAT, 10, 1L, 1L, ntime, lifcr, &status);
	cf_verbose(3, "lifcr array written to timeline table");
    FITS_write_col(outfits, TFLOAT, 11, 1L, 1L, ntime, siccr, &status);
	cf_verbose(3, "siccr array written to timeline table");
    FITS_write_col(outfits, TFLOAT, 12, 1L, 1L, ntime, feccr, &status);
	cf_verbose(3, "feccr array written to timeline table");
    FITS_write_col(outfits, TFLOAT, 13, 1L, 1L, ntime, aiccr, &status);
	cf_verbose(3, "aiccr array written to timeline table");
    FITS_write_col(outfits, TFLOAT, 14, 1L, 1L, ntime, bkgdcr, &status);
	cf_verbose(3, "bkgdcr array written to timeline table");
    FITS_write_col(outfits, TFLOAT, 15, 1L, 1L, ntime, ycentl, &status);
	cf_verbose(3, "ycentl array written to timeline table");
    FITS_write_col(outfits, TFLOAT, 16, 1L, 1L, ntime, ycents, &status);
	cf_verbose(3, "ycents array written to timeline table");

    /* Compute EXPNIGHT for raw data file and write to output header. */
    expnight = 0;
    for (i=0; i<ntime; i++)
	if (!tflag[i]) expnight++;
    FITS_movabs_hdu(outfits, 1, NULL, &status);
    FITS_update_key(outfits, TLONG, "EXPNIGHT", &expnight, NULL, &status);
    cf_verbose(3, "For raw data, EXPNIGHT = %d", expnight);

    free(ptime);
    free(lon);
    free(lat);
    free(limbang);
    free(vorb);
    free(tsunrise);
    free(tsunset);
    free(hv);
    free(lifcr);
    free(siccr);
    free(feccr);
    free(aiccr);
    free(bkgdcr);
    free(tflag);
    free(ycentl);
    free(ycents);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
    return EXIT_SUCCESS; 
}


/* Add new header keywords - needed only for old versions of OPUS */

int cf_add_header_keywords(fitsfile *outfits)
{
    char  hkexists[FLEN_CARD];
    int status=0;
    float dum=0., rawtime;
    short sdum=0;
    long  ldum=0;
    short bkgd_num=0;

    cf_verbose(3, "Entering cf_add_header_keywords");

    /* If keyword RAWTIME does not exist, create with value of EXPTIME */
    fits_read_key(outfits, TFLOAT, "RAWTIME", &rawtime, NULL, &status);
    if (status) {
	status = 0;
        FITS_read_key(outfits, TFLOAT, "EXPTIME", &rawtime, NULL, &status);
        FITS_update_key(outfits, TFLOAT, "RAWTIME", &rawtime,
		"[s] Exposure duration of raw data file", &status);
    }

    /* If keyword HKEXISTS is missing, create with value of "NO" */
    fits_read_key(outfits, TSTRING, "HKEXISTS", hkexists, NULL, &status);
    if (status != 0) {
	status=0;
	cf_verbose(1, "HKEXISTS keyword missing: assume no housekeeping file.");
	FITS_update_key(outfits, TSTRING, "HKEXISTS", "NO",
	    "Housekeeping data file exists", &status);
    }

    FITS_update_key(outfits, TSTRING, "BRIT_OBJ", "NO", "Not an over-bright observation", &status);
    FITS_update_key(outfits, TSTRING, "HSKP_CAL", "        ", "Housekeeping data file", &status);
    FITS_update_key(outfits, TSTRING, "DAYNIGHT", "BOTH", "Use only DAY, NIGHT or BOTH", &status);
    FITS_update_key(outfits, TLONG,  "EXP_HV", &ldum, "[s] Integration time lost to low voltage", &status);
    FITS_update_key(outfits, TLONG,  "EXP_JITR", &ldum, "[s] Integration time lost to jitter", &status);
    FITS_update_key(outfits, TLONG,  "EXPNIGHT", &ldum, "[s] Integration time during night after screening", &status);
    FITS_update_key(outfits, TFLOAT, "FPADXLIF", &dum, "[pixels] Correction for FPA position",&status);
    FITS_update_key(outfits, TFLOAT, "FPADXSIC", &dum, "[pixels] Correction for FPA position",&status);
    FITS_update_key(outfits, TFLOAT, "YCENT1", &dum, "[pixels] Y centroid of LIF HIRS aperture",&status);
    FITS_update_key(outfits, TFLOAT, "YCENT2", &dum, "[pixels] Y centroid of LIF MDRS aperture",&status);
    FITS_update_key(outfits, TFLOAT, "YCENT3", &dum, "[pixels] Y centroid of LIF LWRS aperture",&status);
    FITS_update_key(outfits, TFLOAT, "YCENT5", &dum, "[pixels] Y centroid of SIC HIRS aperture",&status);
    FITS_update_key(outfits, TFLOAT, "YCENT6", &dum, "[pixels] Y centroid of SIC MDRS aperture",&status);
    FITS_update_key(outfits, TFLOAT, "YCENT7", &dum, "[pixels] Y centroid of SIC LWRS aperture",&status);
    FITS_update_key(outfits, TSTRING, "YQUAL1", "        ", "Quality of Y centroid value (LIF HIRS)",&status);
    FITS_update_key(outfits, TSTRING, "YQUAL2", "        ", "Quality of Y centroid value (LIF MDRS)",&status);
    FITS_update_key(outfits, TSTRING, "YQUAL3", "        ", "Quality of Y centroid value (LIF LWRS)",&status);
    FITS_update_key(outfits, TSTRING, "YQUAL5", "        ", "Quality of Y centroid value (SIC HIRS)",&status);
    FITS_update_key(outfits, TSTRING, "YQUAL6", "        ", "Quality of Y centroid value (SIC MDRS)",&status);
    FITS_update_key(outfits, TSTRING, "YQUAL7", "        ", "Quality of Y centroid value (SIC LWRS)",&status);
    FITS_update_key(outfits, TFLOAT, "YGEO_LIF", &dum, "Y centroid of geocoronal lines in target spectrum", &status);
    FITS_update_key(outfits, TFLOAT, "YGEO_SIC", &dum, "Y centroid of geocoronal lines in target spectrum", &status);
    FITS_update_key(outfits, TSHORT, "BKGD_NUM", &bkgd_num, "Number of background regions defined", &status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN0", &sdum, "[pixels] bkgd sample region 0 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN1", &sdum, "[pixels] bkgd sample region 1 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN2", &sdum, "[pixels] bkgd sample region 2 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN3", &sdum, "[pixels] bkgd sample region 3 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN4", &sdum, "[pixels] bkgd sample region 4 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN5", &sdum, "[pixels] bkgd sample region 5 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN6", &sdum, "[pixels] bkgd sample region 6 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN7", &sdum, "[pixels] bkgd sample region 7 - lower limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX0", &sdum, "[pixels] bkgd sample region 0 - upper limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX1", &sdum, "[pixels] bkgd sample region 1 - upper limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX2", &sdum, "[pixels] bkgd sample region 2 - upper limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX3", &sdum, "[pixels] bkgd sample region 3 - upper limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX4", &sdum, "[pixels] bkgd sample region 4 - upper limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX5", &sdum, "[pixels] bkgd sample region 5 - upper limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX6", &sdum, "[pixels] bkgd sample region 6 - upper limit",&status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX7", &sdum, "[pixels] bkgd sample region 7 - upper limit",&status);

    cf_verbose(3, "Exiting cf_add_header_keywords");
    return EXIT_SUCCESS;
}
 
/* Set limits to the background sample regions */

int cf_set_background_limits(fitsfile *outfits) {

    char aper[FLEN_VALUE];
    char bchrfile[FLEN_FILENAME];
    int status=0;
    long npts;
    short bkgd_num=3, *ylim;
    fitsfile *bchrfits;

    cf_verbose(3, "Entering cf_set_background_limits");

/*  Read the target aperture from the input file header. */
    FITS_read_key(outfits,TSTRING,"APERTURE",aper, NULL, &status);
    cf_verbose(3, "aperture = %s", aper); 

/*  If APERTURE = RFPT, use same background region as for LWRS. */
    if (!strncmp (aper, "RFPT", 4)) {
	strncpy(aper, "LWRS", 4);
	cf_verbose(1, "APERTURE = RFPT.  Will treat as LWRS observation.");
    }

/*  Read the name of the BCHR parameter file and open it. */
    FITS_read_key(outfits, TSTRING, "BCHR_CAL",bchrfile, NULL, &status);
    cf_verbose(3,"BCHR parameter file = %s ",bchrfile);
    FITS_open_file(&bchrfits,cf_parm_file(bchrfile), READONLY, &status); 

/*  From the first extension of the BCHR_CAL file, read the limits
    to the background regions.  Close the file */
    FITS_movabs_hdu(bchrfits, 2, NULL, &status);
    npts = cf_read_col(bchrfits, TSHORT, aper, (void **) &ylim);
    FITS_close_file(bchrfits, &status);

    cf_verbose(3, "Number of background regions: %d", bkgd_num);
    cf_verbose(3, "Limits of background regions: %d, %d, %d, %d, %d, %d",
      ylim[0],ylim[1],ylim[2],ylim[3],ylim[4],ylim[5]);

/*  Update the header keywords */
    FITS_update_key(outfits, TSHORT, "BKGD_NUM", &bkgd_num, NULL, &status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN0", &ylim[0], NULL, &status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN1", &ylim[2], NULL, &status);
    FITS_update_key(outfits, TSHORT, "BKG_MIN2", &ylim[4], NULL, &status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX0", &ylim[1], NULL, &status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX1", &ylim[3], NULL, &status);
    FITS_update_key(outfits, TSHORT, "BKG_MAX2", &ylim[5], NULL, &status);

    free (ylim);
    cf_verbose(3, "Exiting cf_set_background_limits");
    return EXIT_SUCCESS ;
}
