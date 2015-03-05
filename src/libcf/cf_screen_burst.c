/****************************************************************************
 *              Johns Hopkins University
 *              Center for Astrophysical Sciences
 *              FUSE
 ****************************************************************************
 *
 * Synopsis:    cf_screen_burst(infits, nevents, ptime, x, y, locflags,
 *                              gti, nsec, ttime, sflag, aic_rate, bkgd) 
 *
 * Description: Determines the times for bursts and sets the burst
 *              flag in the timeline status array so that the busts
 *              can be removed later.
 *
 * Arguments:   fitsfile  infile        Pointer to Intermediate Data File 
 *              long      nevents       Number of photons
 *              float     ptime         Time of event
 *              float     *x, *y        X and Y positions of the photon
 *              byte      locflags       flags on photon location 
 *                                        (LOC_FLGS column in IDF)
 *              GTI       gti           Structure with Good Time Intervals
 *              long      nsec          Length of timeline array
 *              float     ttime         Timeline array
 *              byte      sflag         Status flags in timeline table
 *		float	  aic_rate	AIC count rate
 *              short     bkgd          Count rate for the background 
 *                                      sample region 
 *
 * Algorithm used:
 *
 *     1. Generate an array (flgarr) with 1 second intervals that 
 *        contains flags indicating whether that time interval is good 
 *        or not. All elements of this array are originally set to negative
 *        numbers, indicating bad times.
 *     2. Read in the good time intervals from the initial screening
 *        analysis and set all good time elements in flgarr to positive 
 *        values
 *     3. Determine the median count per binned element in the entire time 
 *        sequence (ignoring times with no data) and set all elements with
 *        a count > 5 times this median value to negative values (indicating
 *        a burst at that time).
 *     4. Do a median filter of the modified time sequence (ignoring the neg
 *        numbers) and identify all times where the count per bin is greater
 *        than a specified value or is more than a specified number of 
 *        standard deviations from the median. These points are flagged with
 *        negative numbers to indicate bursts. The screening parameters used
 *        are specified in the screening parameter file.
 *     5. Iterate time-sequence analysis until no more points are removed.
 *     6. Transfer the burst times from the cnt array to flgarr
 *     7. Go through all detected photons. If the time of arrival for that 
 *        photon has been flagged as bad in FLGARR then ignore it,  
 *        otherwise copy the information to the output file.   
 *
 *
 * Calibration files:
 *	BCHR_CAL: location of the background sample regions
 *	SCRN_CAL: burst screening params
 *
 * Parameter file keywords used to control the burst screening process
 *
 *    MNCNT : minimum enhancement of the background (in counts) needed to 
 *            flag a burst (currently = 5)
 *    STDRE : minimum number of standard deviations from the local mean 
 *            for a burst to be detected (currently=5)
 *    NBIN  : binning factor (in seconds) for the tiome sequence before 
 *            analysis. (currently=15 seconds)
 *    NSMED : Interval (in seconds) used in determining the median filter 
 *            (currently=600)
 *    SRCFRAC :  Fraction of the source intensity required for a 
 *               burst to be removed. SCRFRAC=0.01 requires that a burst 
 *               exceed 1% of the integrated source intensity before being 
 *               removed (default 0.001)
 *
 * Returns:     0 on success
 *
 * History:    10/29/02   v1.1   RDR    Adapted from cf_ttag_burst
 *             11/14/02   v1.2   peb    Tidied code and made it compatible 
 *                                      with cf_screen_burst
 *             12/09/02   v1.3   RDR    - Adjusted procedure for filling
 *                                      burst flags in the timeline
 *                                      - Corrected error in the procedure
 *                                      that tabulates bursts properties.
 *                                      - Used LOCATION flags to exclude
 *                                      geocoronal emission and photons
 *                                      off the active detector.
 *             12/20/02   v1.4   rdr    changed status flag arrays to 
 *                                      unsigned char
 *             01/24/03   v1.5   rdr    corrected problem with high 
 *                                      count rates
 *             03/01/03   v1.6   wvd    Correct use of pointer in FITS_read_key
 *             05/20/03   v1.7   rdr    Add cf_proc_check
 *             05/22/03   v1.8   wvd    Exit if EXPTIME !> 0.
 *					Change cf_error_init to stderr.
 *					Implement cf_verbose throughout.
 *             04/06/03   v1.9   rdr    Exclude times which have previously 
 *                                      had a screening bit set
 *             09/03/03   v1.10  wvd    Delete extraneous print statements
 *             09/10/03   v1.11  wvd    Change background array to type short
 *					Set burst flag only when flgarr[i]
 *					< -1, not 0, to prevent bursts 
 *					from echoing LIMB/SAA/HV flags.
 *             09/15/03   v1.12  wvd    Use ttime[i] to populate bkgd array.
 *					Old scheme didn't work.
 *             10/06/03   v1.13  wvd    Change locflgs to locflags throughout.
 *					Use bit-wise logic to test sflag array.
 *             06/04/04   v1.14  wvd    Clean up i/o.
 *             02/02/05   v1.15  wvd	Read AIC count rate from timeline table
 *					rather than file header keywords.
 *             02/02/05   v1.16  wvd	Generate bkgd array for the entire
 *					exposure, not just good times.
 *             08/30/05   v1.17  wvd	Ignore photons with arrival times
 *					greater than MAX_EXPTIME.
 *             04/09/06   v1.18  wvd	In subroutine median_filter, require
 *					that nwin >= 1.
 *             06/12/06   v1.19  wvd	Call cf_verbose rather than 
 *					cf_if_warning.
 *
 ***************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "calfuse.h"

/*
 *  Function to determine the median value of an array.  
 *  Arguments are the array,
 *  the number of points in the array (npt) and the maximum count to
 *  use in the histogram analysis (nmax). The median value is returned.
 *
 *  The analysis uses a histogram approach to determine the median
 *  value.  Negative values of the array have been flagged as bursts
 *  and should be ignored in determining the median value.  If all of
 *  the points in the array have been flagged, returns zero.
 */

static int
median(int *array, int npts, int nmax)
{
    int i, ngood=0, *hist, medval;

    /*
     * Set up an array to contain the histogram of the data.
     */
    hist = (int *) cf_calloc((nmax+1), sizeof(int));
    /*
     *  Fill the histogram with the data.  
     *  Negative points are not counted.
     */
    for (i=0; i< npts; i++) {
	int ndx = array[i];
	if (ndx > nmax-1)
	    ndx = nmax-1;
	if (ndx >= 0 ) {
	    ngood ++;
	    hist[ndx]++;
	}
    }

    /*
     *  Determine the median value.  Set the data to 0 if there
     *  are fewer than 2 points.
     */
    if (ngood > 2 ) {
	int hnum= (int) ngood/2, atot=0;
	i=-1;
	while ((atot <= hnum) && (i < nmax) )
	    atot += hist[++i];
	medval = (int) i;
    }
    else
	medval = 0;

    free(hist) ;

    return medval;
}

/* 
 *  Function to construct a median filter for an array:
 *  array  = array for which the median filter is to be constructed
 *  npts   = number of points in the array
 *  nflt   = number of points over which to determine the median value
 *  medflt = array containing the median filter.
 */

static void
median_filter(int *array, int npts, int nflt, int *medflt)
{
    int cntmax=0, i, hwidth, nwin, *win, medval, nsam ;
    /*
     *  Determine the maximum value of the array.
     */
    for (i = 0; i < npts; i++)
	if (array[i] > cntmax)
	    cntmax=array[i];

    cf_verbose(3,"Constructing median filter - max cnt = %d",cntmax) ;
    
    /*
     *  Determine the properties of the window over which to calculate
     *  the median values. Make sure that the window is not larger
     *  than the width of the data.
     */
    if (nflt > npts)
	hwidth = (int) ( npts/2. ) - 1;
    else
	hwidth = (int) (0.5 + nflt/2. );
    nwin = (int) 2*hwidth+1;
    if (nwin < 1) nwin = 1;
    win = (int *) cf_calloc(nwin, sizeof(int));
    /*
     *  Construct the median filter.
     */
    for (i = 0; i < npts; i++) {
	int j, nmin, nmax;
	/*
	 *  Determine the limits to the sample.
	 */
	nmin = i - hwidth;
	if (nmin < 0)
	    nmin=0;
	nmax = nmin+nwin;
	if (nmax > npts-1 )
	    nmin = npts-nwin-1;

	/* Fill the sample array */
        nsam=0 ;
	for (j=0; j<=2*hwidth; j++) {  
	    int ndx = nmin+j;
	    if (ndx > npts-1)
		ndx = npts-1;
	    win[j] = array[ndx];
            if (array[ndx] > 0) nsam++ ;
	}

	/* Get the median value for the sample array */
	medval = median(win, nwin, cntmax);

	/*
        cf_verbose(6,"at point %d, nmin=%d, nmax=%d, nsam=%d, medval=%d",
		   i, nmin, nmax, nsam, medval) ;
	*/

	/* Update the median filter */
	medflt[i] =  medval;
    }

    free(win) ;

}

/*
 *  Function to take an array marked with the times and intensities of
 *  the bursts and identify individual bursts events and output these
 *  to an ASCII file for later analysis.
 */

static void
tab_burst(int *bursts, int npts, int nbin, int mjd, long tst, char *outfile)
{
    FILE *output;
    int i, bflg=-1, bdur=0, bcnts=0, bst=0;
    /*
     *  Go through the bursts array and determine the times and
     *  properties of the bursts.
     */
    output = fopen(outfile, "w");

    cf_verbose(2, "BURST SUMMARY") ;
    for (i=0; i<npts; i++) {
	if (bursts[i] > 0 && bflg < 0) {
	    bflg=1;
	    bst= i*nbin;
	}
	if (bursts[i] > 0 && bflg > 0) {
	    bdur++;
	    bcnts = bcnts + bursts[i];
	}
	if (bursts[i] < 0 && bflg > 0) {
	  cf_verbose(2, "start=%d,  duration=%ld seconds,  counts=%d",
		 bst, bdur*nbin, bcnts) ; 
 	   fprintf(output, " %d   %ld    %d   %d \n",
	       mjd, bst+tst, bdur*nbin, bcnts);
  	    bflg=-1;
	    bdur=0;
	    bcnts=0;
	}
    }

    /* Check to see whether a burst is in progress at the end of the obs */
    if (bflg > 0) {
       	  cf_verbose(2, "start=%d,  duration=%ld seconds,  counts=%d",
		 bst, bdur*nbin, bcnts) ; 
	  fprintf(output, " %d   %ld    %d   %d\n",
		    mjd, bst+tst, bdur*nbin, bcnts);
    }

}


int
cf_screen_burst(fitsfile *infits, long nevents, float *ptime, float *x, 
		float *y, unsigned char *locflags, GTI *gti, long nsec, 
                float *ttime, unsigned char *sflag, float *aic_rate,
		short *bkgd)
{
    char CF_PRGM_ID[] = "cf_screen_burst";
    char CF_VER_NUM[] = "1.19";

    short ylim[8], k;
    int   anynull, errflg=0, intnull=0, cflg=-1, nrowbkg=0, ngood;
    int   npts, nmax, nflt, iflg, numflg, ndxs, ndxe, ndx, ndxb;
    int   *bcnts, *cnts, *gcnts, *medflt, *flgarr, *bursts, lim_col;
    int   nflg, niter, mjd , status=0, medval;
    long  frow=1, felem=1, tst , i, j;
    long  stdrej, mncnt, nsmed, mnburst, nbin, nflgarr;
    float exptime, max_aic_rate, bkgsf, bkgsft, srcfrac, mxtime;
    double expstrt;
    char  burst_tab[FLEN_CARD];
    char  expid[FLEN_CARD], progid[FLEN_CARD], targid[FLEN_CARD];
    char  scobsid[FLEN_CARD], scrnfile[FLEN_CARD], aperf[FLEN_CARD];
    char  det[FLEN_CARD], bchrfile[FLEN_CARD], aper[5];
    fitsfile *scrnfits, *bchrfits;

    /*
     *  Initialize error checking
     */ 
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    /*
     *  Read header keywords.
     */
    FITS_read_key(infits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    FITS_read_key(infits, TDOUBLE, "EXPSTART", &expstrt, NULL, &status);
    FITS_read_key(infits, TSTRING, "PRGRM_ID", progid, NULL, &status);
    FITS_read_key(infits, TSTRING, "TARG_ID", targid, NULL, &status);
    FITS_read_key(infits, TSTRING, "SCOBS_ID", scobsid, NULL, &status);
    FITS_read_key(infits, TSTRING, "EXP_ID", expid, NULL, &status);
    FITS_read_key(infits, TSTRING, "DETECTOR", det, NULL, &status);

    /*
     * Exit if EXPTIME is 0.
     */ 
    if (exptime < 1) {
	cf_verbose(1, "EXPTIME = %f.  Exiting.", exptime);
	return status;
    }

    /*
     *  Determine the Julian date and start time from expstrt.
     */
    mjd = (int) expstrt;
    tst = (int) ((expstrt - mjd)*86400.);
    det[1] = tolower(det[1]);

    /*  Ignore events with arrival times greater than MAX_EXPTIME. */
    i = nevents - 1;
    while (ptime[i] > MAX_EXPTIME && i > 0) i--;
    nevents = i + 1;
    i = nsec - 1;
    while (ttime[i] > MAX_EXPTIME && i > 0) i--;
    nsec = i + 1;
  
    /*
     *  Determine the maximum count rate for the exposure.
     *  If the count rate is too high, there is the possibility of
     *  losing counts and of dropouts - Check the count rate and set
     *  a flag if it is too high.
     */
    max_aic_rate = -999.;
    for (i = 0; i < nsec; i++)
	if (aic_rate[i] > max_aic_rate) max_aic_rate = aic_rate[i];
    cf_verbose(2, "Max AIC count rate = %d", cf_nint(max_aic_rate));
    if(max_aic_rate > 4000) cflg=1;

    /*
     *  Read the screening parameter file.
     */
    FITS_read_key(infits, TSTRING, "SCRN_CAL", scrnfile, NULL, &status);
    cf_verbose(3, "Screening parameters from %s\n", scrnfile);
    FITS_open_file(&scrnfits, cf_parm_file(scrnfile), READONLY, &status);
    FITS_read_key(scrnfits, TLONG, "MNCNT", &mncnt, NULL, &status);
    FITS_read_key(scrnfits, TLONG, "STDREJ", &stdrej, NULL, &status);
    FITS_read_key(scrnfits, TLONG, "NBIN", &nbin, NULL, &status);
    /*
     *  Set the binning factor to 1 for high count rate data.
     */
    if (cflg > 0) {
        nbin=1;
        cf_verbose(1, "High count-rate observation: using 1 sec time bins") ;
    }
    FITS_read_key(scrnfits,TLONG,"NSMED",&nsmed, NULL, &status);
    /*
     *  srcfrac=required fraction of the source intensity for a
     *  burst to be counted, i.e. the burst must exceed this fraction
     *  before it is deemed to be significant.
     */
    FITS_read_key(scrnfits, TFLOAT, "SRCFRAC", &srcfrac, NULL, &status);
    FITS_read_key(infits, TSTRING, "APERTURE", aperf, NULL, &status);
    FITS_close_file(scrnfits, &status);

    cf_verbose(3, "Screening parameters used:");
    cf_verbose(3, "Minimum count rate = %ld", mncnt);
    cf_verbose(3, "Detection criterion (num sigma) = %ld", stdrej);
    cf_verbose(3, "Number of seconds to bin the time sequence = %ld", nbin);
    cf_verbose(3, "Interval (seconds) used in median filter = %ld", nsmed);
    cf_verbose(3, "Fraction of source intensity for burst = %f\n", srcfrac);
    /*
     *  Select the background extraction regions - RFPT implies
     *  that there is no source, so the entire detector can be used.
     */
    strncpy(aper, aperf, 5);
    cf_verbose(3, "Aperture designation = %s", aper);
    if (strcmp ( aper , "RFPT") == 0 ) {
	cf_verbose(3, "APERTURE = RFPT.  No source in aperture.");
	ylim[0] = 0;
	ylim[1] = 1;
	ylim[2] = 2;
	ylim[3] = 1020;
	ylim[4] = 1021;
	ylim[5] = 1022;
	ylim[6] = 1022;
	ylim[7] = 1022;
    }
    else {
	/*
	 *  open the bchr parameter file
	 */
	FITS_read_key(infits, TSTRING, "BCHR_CAL", bchrfile, NULL, &status);
	cf_verbose(3, "BCHR parameter file: %s", bchrfile);
	FITS_open_file(&bchrfits, cf_parm_file(bchrfile), READONLY, &status);

	/* Move to the first extension of the BCHR file, which contains
	   the limits to the background regions */
	FITS_movabs_hdu(bchrfits, 2, NULL, &status);

	/* Get the column number for the data to be extracted */
	FITS_get_colnum(bchrfits, CASEINSEN, aper, &lim_col, &status);
	
	/* Read the data from the table */
	FITS_read_col(bchrfits, TSHORT, lim_col, frow, felem, 8, &intnull, 
		      ylim, &anynull, &status);

	/* Close the bchr file */
	FITS_close_file(bchrfits, &status);
    }
    cf_verbose(3, "Y limits of background region: %d-%d, %d-%d, %d-%d\n",
	   ylim[0], ylim[1], ylim[2], ylim[3], ylim[4], ylim[5]);
    /*
     *  Now get the time sequence and a good time interval flag array.
     */
    cf_verbose(3, "nevents = %ld", nevents);
    mxtime=ptime[nevents-1];
    nflgarr = cf_nint(mxtime);
    npts= (int) mxtime/nbin;
    cf_verbose(3, "Number of points in the (binned) time sequence = %d", npts);
    cf_verbose(3, "Number of points in the good time flag array = %ld", nflgarr);
    /*
     *  Generate arrays to contain the time sequence in the
     *  background sample regions (cnts) as well as the entire
     *  detector (gcnts). Also generate an array to contain only
     *  the bursts (bursts).
     *  Array bcnts will contain the background count rate.
     */
    bcnts  = (int *) cf_calloc(npts, sizeof(int));
    cnts   = (int *) cf_calloc(npts, sizeof(int));
    gcnts  = (int *) cf_calloc(npts, sizeof(int));
    bursts = (int *) cf_calloc(npts, sizeof(int));
    /* Generate an array to contain the good time flags - this is an array 
       with 1s time intervals. Each interval is set to a negative number (-1)
       if it contains a bad time interval or a positive number (1) if it 
       contains a good time interval. The contents of this array are used to
       screen the photons and to determine the good time intervals */
    flgarr = (int*) cf_calloc(nflgarr, sizeof(int) );
    /* Generate an array with 1s time intervals to contain the day/night
       information - this will be 1 for night and 0 for day */

    /* Initially flag all of the points with a negative number.
       The good time intervals will be determined later */
    for (i=0; i<npts; i++) {
      cnts[i] = -10 ;
      gcnts[i] = 0 ;
      bursts[i] = -10 ;
    }
    for (i=0; i<nflgarr; i++) flgarr[i] = -1 ;
    /*
     *  Set the counts array to zero (i.e. unflag the point) for all
     *  good time points.
     */
    ngood = 0 ;
    for (i=0; i < gti->ntimes;  i++) {
	 ndxs= (int) gti->start[i]/nbin;
	   if (ndxs < 0 ) ndxs = 0 ;
	 ndxe = (int) (0.5 + gti->stop[i]/nbin);
           if (ndxe > npts-1) ndxe = npts-1;
	for (j=ndxs; j<= ndxe; j++) {
	    ngood++;
	    cnts[j]=0;
	}
    }
    cf_verbose(2, "Number of good time points = %d\n", ngood);

    /* Set the good time flags array to a positive value for all good 
       time points  */

    for (i=0; i < gti->ntimes;  i++) {
	ndxs= (int) gti->start[i];
	    if (ndxs < 0) ndxs = 0 ;
	ndxe = (int) gti->stop[i];
	   if (ndxe > nflgarr-1) ndxe=nflgarr-1;
	for (j=ndxs; j<= ndxe; j++)
	    flgarr[j]=1;
    }

    /* Check the timeline-table screening flag and mark all times that have
       a screening bit set (except for day/night) */
     for (i=0; i<nsec; i++) {
       ndx = (int) (0.5+ttime[i]) ;
	if( ndx > nflgarr-1) ndx = nflgarr-1;
        if( ndx < 0) ndx = 0 ;
       ndxb = (int) (0.5 + (float) ndx/ (float) nbin) ;
        if(ndxb > npts-1) ndxb = npts-1 ;
        if (ndxb < 0) ndxb=0 ;
       if (sflag[i] & ~TEMPORAL_DAY) {
	  flgarr[ndx] = -1 ;
          cnts[ndxb] = -10 ;
	}
    }

    /* Generate the time sequence for the total counts from the detector. 
       Use only the active area and ignore geocoronal emission. */

    for (i=0; i < nevents; i+=1)
    if( (locflags[i] & LOCATION_SHLD) == 0 && 
        (locflags[i] & LOCATION_AIR) == 0 ) {
	    int ndx = (int) ptime[i]/nbin;
	    if (ndx > npts-1)
		ndx=npts-1;
	    if (ndx < 0)
		ndx=0;
	    gcnts[ndx] += 1;
	}
    
    /* Determine the number of rows in the background sample region */
    for (k=0; k < 3; k++)
	nrowbkg += ylim[2*k+1] - ylim[2*k];
    /*  
     *  Determine the scaling factor needed to convert the measured
     *  counts in the background regions to an estimate of the
     *  background counts over the active aperture.
     */
    bkgsf= (float) ylim[6]/nrowbkg;
    bkgsft = (float) ylim[7]/nrowbkg;

    cf_verbose(2, "Rows in target aperture = %d", ylim[6]);
    cf_verbose(2, "Rows on the detector = %d", ylim[7]);
    cf_verbose(2, "Rows in background region = %d", nrowbkg);
    cf_verbose(2, "Scale factor to aperture = %f ", bkgsf);
    cf_verbose(2, "Scale factor to entire detector = %f", bkgsft);
        
    /* Generate the time sequence using only the counts from the three
       selected background regions - ignore geocoronal lines and times 
       that have been flagged as bad (with negative count values) */ 
    for (k=0; k<3; k++) {  
	short ylow=ylim[2*k];
	short yhigh=ylim[2*k+1];
	for (i=0; i < nevents; i+=1) {
	  if( (locflags[i] & LOCATION_SHLD) == 0 && 
              (locflags[i] & LOCATION_AIR) == 0 && 
              (y[i] > ylow) && (y[i] < yhigh) ) {
		int ndx = (int) ptime[i]/nbin;
		if (ndx > npts-1)
		    ndx=npts-1;
		if (ndx < 0)
		    ndx=0;
		if (cnts[ndx] >= 0) cnts[ndx] += 1;
		bcnts[ndx] ++;
	    }
	}
    }
    
    /* Fill background count rate array with cnts array, and divide by nbin
       to give count rate */
/*	THIS SCHEME WORKS ONLY IF THERE ARE NO TIME BREAKS IN THE TIMELINE
 *	TABLE.  SINCE TIME BREAKS ARE COMMON, MUST USE TTIME ARRAY TO MAP
 *	CNTS ARRAY TO BKGD ARRAY.	WVD 09/15/03
 *   for (i=0; i<npts; i++) {
 *	int ndxs= i*nbin;
 *	int ndxe= ndxs + nbin;
 *	if (ndxe > nsec-1)
 *	    ndxe=nsec-1;
 *	for (j=ndxs;j<ndxe; j++)
 *	    bkgd[j] = (short) ((float) cnts[i] / nbin + 0.5);
 *   }
 */
    for (i=0; i<nsec; i++) {
       ndx = cf_nint(ttime[i]/nbin) ;
	if( ndx > nflgarr-1)
	    ndx = nflgarr-1;
/*	if (cnts[ndx] > 0.)
  	    bkgd[i] = (short) ((float) cnts[ndx] / nbin + 0.5);
	else
	    bkgd[i] = 0;
*/
	bkgd[i] = bcnts[ndx] / nbin;
    }
    free (bcnts);

    /* Remove the background counts (scaled to the entire detector) from the 
       measured gross count rate  */
    for (i=0; i < npts; i++) if (cnts[i] > 0) {
	gcnts[i] = gcnts[i] - bkgsft*cnts[i];
	if (gcnts[i] < 0)
	    gcnts[i] = 0;
    }

    /* Check for zeros in the count array for high count rate data - these
       indicate dropouts. If one is found, flag it with a negative number so
       that the program will ignore it in future analysis */
    numflg=0;
     if (cflg > 0) {
	for (i=0; i < npts; i++) {
	    if (cnts[i] <= 0.1 ) {
		numflg += 1;
		cnts[i]=-10;
	    }
	}
    }
    cf_verbose(2, "Number of dropouts found = %d", numflg);

    /* Determine the median value of the array */
    nmax=2000;
    medval = median(cnts, npts, nmax);
    cf_verbose(2, "Median value of the entire original array = %d ", medval);

    /* Screen the time series for very large count rates (> 3 x median) */
    cf_verbose(2, "Searching for large bursts ") ;
    for ( i=0; i<npts; i++)
	if (cnts[i] > 4*medval) {
	  cf_verbose(2,"point rej: time=%d, cnts=%d, 4*medval=%d ",
		     i*nbin, cnts[i], 4*medval) ;
	    bursts[i] = cnts[i]-medval;
	    cnts[i] = -100;
	}

    /* Do a median filter on the pre-filtered array */

    /* Generate the array to contain the median filter */
    medflt = (int *) cf_calloc(npts, sizeof(int) );

    /*
     *  Now specify the parameters of the filter - 
     *  nsmed = number of seconds over which to determine median values
     *  nflt = number of points in the cnts array over which to
     *  determine the median value
     */
    nflt = nsmed / nbin;
    niter=0;
 
    cf_verbose(3,"Searching for small bursts - nflt = %d", nflt) ;

    do {
	median_filter(cnts, npts, nflt, medflt);

	/* Throw out all points which are more than (stdrej) standard
	   deviations above the median. Also require that the counts exceed
	   a certain minimum value. */
	nflg=0;
	mnburst=mncnt*nbin;
	for (i=0; i<npts; i++) {
	    int dcnt=cnts[i]-medflt[i];
	    if ( (dcnt > mnburst) && 
		 (dcnt > stdrej*sqrt((double) cnts[i])) && 
		 (dcnt*bkgsf > srcfrac*gcnts[i]) ) {
		nflg++;
		cf_verbose(2, "point rej: time=%ld, cnts=%d, medflt=%d, "
		       "dif=%d, gcnts= %d",
		       i*nbin, cnts[i], medflt[i], dcnt, gcnts[i]);
		cnts[i] = -100;
		bursts[i] = dcnt;
	    }
	}
	niter++;
	cf_verbose(3,"After iteration %d, we found %d times affected by bursts",
                   niter, nflg) ;
	/*  Iterate the median filter if any additional points have
	 *  been modified*/ 
    } while(nflg > 0 && niter < 10);

    /* Print out the results of the screening 
       cf_verbose(5, " "); 
       for (i=0; i<npts; i++) 
           cf_verbose(5, "At time %4d, cnts= %4d, medflt= %4d, dif= %4d, gcnts= %4d",
		i*nbin, cnts[i],medflt[i], cnts[i]-medflt[i], gcnts[i]);
       cf_verbose(3, " "); 
    */

    /* Use the regions of the cnts array equal to -100 to specify
       the burst times in the flgarr array */
    for (i=0; i<npts; i++ ) {
	if (cnts[i] < -20) {
	    int ndxs = i*nbin, ndxe;
	    if (ndxs > nflgarr-nbin-2) 
		ndxs=nflgarr-nbin-2;
	    ndxe = ndxs + nbin - 1;
	    for (j=ndxs; j<=ndxe; j++)
		flgarr[j] = -10;
	}
    }

    /* tabulate the burst times */
    sprintf(burst_tab, "%4s%2s%2s%3s%2s_bursts.dat", progid, targid, scobsid, expid, det);
    tab_burst(bursts, npts, nbin, mjd, tst, burst_tab);
    cf_verbose(3, "Burst properties written to file: %s", burst_tab);

    /* Check the last 2*nbin elements of flgarr. If any of these are set as
       bad, then set all the end points as bad. The binned count time sequence 
       will not sample the end of the ovarall time sequence if there is not an 
       integral number of bins present. We also do not need the last bin if the
       second-to-last is bad. */

    i=nflgarr-2*nbin;
    iflg=10;
    do{ 
	if (flgarr[i] < -1) {		/* 09/12/03 - wvd - Change from 0 */
	    for (j=i; j<nflgarr; j++)
		flgarr[j]=-10;
	    iflg = -10;
	};
	i++;
    } while (iflg > 0 && i < nflgarr);

    /* set the timeline status flag (sflag - TEMPORAL_BRST) for all 
       times containing bursts */
    for (i=0; i<nsec; i++) {
       ndx = (int) (0.5+ttime[i]) ;
	if( ndx > nflgarr-1)
	    ndx = nflgarr-1;
	if (flgarr[ndx] < -1)		/* 09/12/03 - wvd - Change from 0 */
	    sflag[i] = sflag[i] | TEMPORAL_BRST;
    }
    /* deallocate storage space */
    free(cnts) ;
    free(gcnts) ;
    free(flgarr) ;
    free(medflt) ;
    free(bursts) ;

    /* Enter a time stamp into the log */
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing");

    return (status);
} 
