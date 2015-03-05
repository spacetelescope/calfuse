/**************************************************************************
 *           Johns Hopkins University
 *           Center for Astrophysical Sciences
 *           FUSE
 **************************************************************************
 *
 * Synopsis:    cf_calculate_y_centroid(infits, nevents, weight, x, y,
 *			channel, timeflags, locflags)
 *
 * Description: Determines the y centroid of the target and airglow spectra
 *              in each aperture.  The value written to the header depends
 *		on the quality flag for each channel in the IDF header.
 *		If QUALITY = HIGH, the target centroid is used; if MEDIUM,
 *		the airglow centroid is used; if LOW, the default centroid
 *		(from the CHID_CAL file) is used.
 *
 * Arguments:   fitsfile infits   Pointer to the input Intermediate Data File
 *              long     nevents  Number of photon events in the file
 *		float	 weight   Scale factor for each photon
 *              float    *x, *y   X and Y positions of the photon
 *       unsigned char  *channel  Aperture associated with each photon
 *       unsigned char *timeflags Time flags - used to identify photons
 *                                  arriving during bursts, etc.
 *       unsigned char *locflags  Location flags - used to identify photons
 *                                  in geocoronal line regions
 *
 * Calibration files required: 
 *
 * Returns: 0 on success
 *                                
 * HISTORY:     10/16/02    v1.1   RDR  started work
 *              12/02/02    v1.2   RDR  cleaned up variables
 *		12/12/02    v1.3   wvd	added include files
 *              12/26/02    v1.4   RDR  corrected error in calculating centroid
 *              02/24/03    v1.5   peb  changed include file to calfitsio.h
 *		03/11/03    v1.6   wvd  Changed screen to unsigned char
 *              03/25/03    v1.7   rdr  ignore pinhole aperture
 *		04/03/03    v1.9   wvd  Simplify calculation of centroid.
 *					Include airglow in centroid
 *					calculation of non-target apertures.
 *					Replace printf() with cf_verbose().
 *		04/17/03    v1.10  wvd  Fix bug in calculation of non-target
 *					centroids.
 *		04/18/03    v1.11  wvd  Call cf_proc_update only if final_call
 *					is TRUE.  Scale YCENT by photon weight
 *					for use with HIST data.  Changed
 *					screen to locflags.
 *              05/20/03    v1.12  rdr  Added call to cf_proc_check
 *              05/22/03    v1.14  wvd	cf_error_init to stderr
 *              08/21/03    v1.16  wvd	Change channel to unsigned char
 *              09/02/03    v1.17  bjg  Now read user-specified centroid 
 *                                      information from PARM_CAL and
 *                                      expected centroids from CHID_CAL.
 *              09/02/03    v1.18  wvd	Tidy up code.
 *              09/17/03    v1.19  wvd	Return -1 if a target centroid
 *					is too far from expected value.
 *              10/28/03    v1.20  wvd	Rewrite program with new logic.
 *					Delete use of final_call.
 *              11/12/03    v1.21  wvd	If no photons in channel, don't
 *					crash, just use default centroid.
 *					If quality = LOW, ygeo = ycent.
 *              04/09/04    v1.22  bjg  Include string.h and stdio.h
 *              06/04/04    v1.23  wvd	Use photon time flags to exclude
 *					bursts, etc. from calculation.
 *              08/18/04    v1.24  wvd	Test all spectra to see whether
 *					centroid is out of bounds.
 *		03/11/05    v1.25  wvd  Tabulate centroids as doubles,
 *					then convert to floats.
 *					Write target centroid values
 *					to trailer file if verbose >= 2.
 *		04/15/05    v1.26  wvd  Clean up i/o.
 *		04/26/05    v1.27  wvd  Don't populate YGEO keywords, as
 *					they are wrong for point sources.
 *					Consider only events with good
 *					LOCATION flags.
 *					Use the largest aperture with a valid
 *					YGEO to determine offset from tabulated
 *					position.  Use this shift to calculate 
 *					ygeo for the smaller apertures. 
 *		05/18/05    v1.28  wvd  Test both SPEX_LIF and SPEX_SIC.
 *					Don't override user-specified centroid,
 *					even if it's far from expected value.
 *					Fix bug in tabulation of expected
 *					Y centroids.
 *		02/01/06    v1.29  wvd  If ycent is too far from default value,
 *					use default value, not airglow centroid.
 *					Allow computation of yshift from
 *					airglow photons in target aperture.
 *		11/25/08    1.30 wvd	In HIST mode, there are no spectra in
 *					other apertures to cause confusion,
 *					so we need not require that the target
 *					centroid lie within X pixels of either
 *					the airglow or tabulated centroid.
 *
 **************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "calfuse.h"

int cf_calculate_y_centroid(fitsfile *infits, long nevents, float *weight,
     float *x, float *y, unsigned char *channel, unsigned char *timeflags,
     unsigned char *locflags)
{

     char CF_PRGM_ID[] = "cf_calculate_y_centroid";
     char CF_VER_NUM[] = "1.30";
 
     char aperture[FLEN_VALUE], instmode[FLEN_VALUE];
     char ycentname[FLEN_KEYWORD], file_name[FLEN_VALUE];
     char key[FLEN_KEYWORD], quality[8][FLEN_VALUE];
     double dncent, dngeo, dycent, dygeo;
     float ycent, ygeo, yshift, ycent_tab[8];
     int chan_num, errflg=0, status=0, active_ap[2];
     int got_shift, i, hdu, target_ap, user;
     int spex_sic, spex_lif, emax_sic, emax_lif, emax, extended;
     long n;
     fitsfile *parm_cal, *chid_cal;

  /* Initialize error checking */ 
     cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

  /* Enter a time stamp into the log */
     cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

  /* Check whether routine is appropriate for this data file. */
     if ( (errflg = cf_proc_check(infits, CF_PRGM_ID)) ) return errflg;

  /* Read the source aperture from the file header. */
     FITS_read_key(infits, TSTRING, "APERTURE", aperture, NULL, &status);
     extended=cf_source_aper(infits, active_ap) ;

  /* Were data taken in HIST or TTAG mode? */
     FITS_read_key(infits, TSTRING, "INSTMODE", instmode, NULL, &status);
     
  /* Read centroid quality keywords from file header. */
     for (i = 1; i < 8; i++) {
	if (i == 4) continue;
	sprintf(key, "YQUAL%1d", i);
	FITS_read_key(infits, TSTRING, key, &quality[i][0], NULL, &status);
     }
     
  /* Read extraction parameters from PARM_CAL file. */
     FITS_read_key(infits, TSTRING, "PARM_CAL", file_name, NULL, &status);
     FITS_open_file(&parm_cal, cf_parm_file(file_name), READONLY, &status);
     FITS_read_key(parm_cal,TINT,"SPEX_SIC",&spex_sic,NULL,&status);
     FITS_read_key(parm_cal,TINT,"SPEX_LIF",&spex_lif,NULL,&status);
     FITS_read_key(parm_cal,TINT,"EMAX_SIC",&emax_sic,NULL,&status);
     FITS_read_key(parm_cal,TINT,"EMAX_LIF",&emax_lif,NULL,&status);
     FITS_close_file(parm_cal,&status);
     
  /* Read expected centroids from CHID_CAL file.  Use extended-aperture
     values except for apertures containing point-source targets.  */
     FITS_read_key(infits, TSTRING, "CHID_CAL", file_name, NULL, &status);
     FITS_open_file(&chid_cal, cf_cal_file(file_name), READONLY, &status);
     for (i = 1; i < 8; i++) {
     	if (i==4) continue;
        if ((i == active_ap[0] || i == active_ap[1]) && !extended) hdu=i+1;
	else hdu=i+9;
     	FITS_movabs_hdu(chid_cal, hdu, NULL, &status);
     	FITS_read_key(chid_cal,TFLOAT,"CENTROID",ycent_tab+i,NULL,&status);
     }
     FITS_close_file(chid_cal,&status);

  /* Determine the centroids for each channel */
     yshift = 0;
     got_shift = FALSE;
     for (chan_num = 7; chan_num > 0; chan_num--) {
	if (chan_num == 4) {
	    yshift = 0;
	    got_shift = FALSE;
	    continue;
	}
        dngeo = dncent = 0.;
	dygeo = dycent = 0.;
	sprintf (ycentname, "YCENT%d", chan_num);

        /* Is this a target aperture? */
        if (chan_num == active_ap[0] || chan_num == active_ap[1])
            target_ap = TRUE;
        else
            target_ap = FALSE;

	/* If not, and we're in histogram mode, skip to the next channel. */
	if (!target_ap && !strncmp(instmode, "HIST", 4)) continue;

	/* Calculate separate target and geocoronal centroids. */
        for (n = 0; n < nevents; n++) {
	    if (channel[n] == chan_num && !(timeflags[n] & ~TEMPORAL_DAY) &&
		!(locflags[n] & ~LOCATION_AIR)) { 

		if(locflags[n] & LOCATION_AIR) {  /* Airglow photon */
		    dygeo += y[n] * weight[n];
		    dngeo += weight[n];
		}
		else {				  /* Target photon */
		    dycent += y[n] * weight[n];
		    dncent += weight[n];
		}
	    }
        }

	/* If we can't find any photons but expect to, set quality to LOW. */
	if ((quality[chan_num][0] == 'H' && dncent < 1) ||
	    (quality[chan_num][0] == 'M' && dngeo  < 1)) {
	    quality[chan_num][0] = 'L';
	    sprintf(key, "YQUAL%1d", chan_num);
	    FITS_update_key(infits, TSTRING, key, "LOW", NULL, &status) ;
	    cf_verbose(1, "No photon events in channel %d", chan_num);
	}

        /* Compute centroids of target and airglow spectra */
	ycent = ygeo = 0; 
	if (dncent > 0) ycent = dycent / dncent; 
	if (dngeo  > 0) ygeo  = dygeo / dngeo; 

	/* If yshift is already known, use it to calculate ygeo.
	 * If not, and ygeo is valid for this aperture, compute yshift.
	 */
	if (got_shift)
	    ygeo = yshift + ycent_tab[chan_num];
	else if (dngeo > 0) {
	    yshift = ygeo - ycent_tab[chan_num];
	    got_shift = TRUE;
	}

	/* Now decide which values to write to the header. */
	  
	/* If we are in the target aperture and the user has specified
         * YCENT for this channel, use it.
	 */
        if (target_ap && ((chan_num<5 && spex_lif>=0 && spex_lif<1024) ||
            (chan_num>4 && spex_sic>=0 && spex_sic<1024))) {
            if (chan_num<5) ycent=spex_lif;
            else ycent=spex_sic;
	    cf_verbose(1, "Channel %d: User-specified Y centroid = %g",
		chan_num, ycent);
	    user = TRUE;
	}
	else
	    user = FALSE;

	/* If the centroid quality is HIGH, use the calculated target
	 * value, regardless of whether we are in a target aperture.
	 * (This happens by default, so we don't have to do anything.)
	 */

	/* If the centroid quality is MED, use the airglow centroid. */
	if (quality[chan_num][0] == 'M')
	    ycent = ygeo;

	/* If the quality is LOW, use the default value. */
	else if (quality[chan_num][0] == 'L')
	    ycent = ycent_tab[chan_num];

        /* 
	 *  Test whether YCENT is too far from the expected value.
	 *  If so, use the tabulated value and set the quality to LOW.
	 *  Don't test HIST data. (wvd, 11/25/2008)
	 */
	if (chan_num<5) emax=emax_lif; else emax=emax_sic;
	if (fabs(ycent-ycent_tab[chan_num]) > emax && !user && !strncmp(instmode, "TTAG", 4)) {
	    if(target_ap) {
        	cf_verbose(1, "Channel %d: computed centroid (%.1f) is too far "
		    "from expected value.", chan_num, ycent);
		cf_verbose(1, "Channel %d: using default centroid of %.1f",
		    chan_num, ycent_tab[chan_num]);
	    }
	    else {
		cf_verbose(2, "Channel %d: using default centroid of %.1f",
		    chan_num, ycent_tab[chan_num]);
	    }
	    ycent = ycent_tab[chan_num];
	    sprintf(key, "YQUAL%1d", chan_num);
	    FITS_update_key(infits, TSTRING, key, "LOW", NULL, &status) ;
	}
	    
        /* Write centroid of target spectrum to the header */
        FITS_update_key(infits, TFLOAT, ycentname, &ycent, NULL, &status) ;
        cf_verbose(2, "For channel %d, Y centroid = %.1f", chan_num, ycent) ;

     }			/* End of loop over channels. */

     cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
     cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing"); 

     return status;
}
