/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE 
 *****************************************************************************
 *
 * Synopsis:    cf_find_spectra(fitsfile *header, long nevents,
 *		float *weight, float, *xfarf, float *yfarf, char *chan,
 *		unsigned char *timeflags, unsigned char *locflags,
 *		int airglow_centroid)
 *
 * Description: Determines the spectral centroid for each aperture.
 *		First, we compute the airglow centroid for the LWRS
 *		aperture and determine its offset relative to the
 *		tabulated value.  For the MDRS and HIRS apertures, we
 *		apply this offset to the tabulated centroid to get the
 *		airglow centroid.  For each aperture, we search for a
 *		target spectrum in a narrow range around the airglow
 *		centroid.  If we find it, we use it; otherwise, we use
 *		the airglow centroid.  If airglow_centroid = TRUE, use
 *		airglow centroid for the target spectrum.
 *
 *
 * Arguments:   fitsfile  *header       Pointer to the location of the FITS
 *                                      header of the Intermediate Data File
 *              long      nevents       Number of photons in the file
 *              float     *weight       scale factor for each photon event
 *              float     *xfarf        X position of each photon event
 *              float     *yfarf        Y position of each photon event
 *          unsigned char *chan         Pointer to the array containing
 *                                      channel information for each event
 *          unsigned char *timeflags    Pointer to the array containing
 *					time flags for each event
 *          unsigned char *locflags     Pointer to the array containing
 *					location flags for each event
 *		int  airglow_centroid	If TRUE, use airglow lines to compute
 *					centroid of target spectrum.
 *
 * Calls:       
 *
 * Return:      0 on success
 *
 * History:     04/18/03    1.1  wvd    Initial work
 *              05/20/03    1.2  rdr    Added call to cf_proc_check 
 *              07/29/03    1.3  rdr    Deal with case where no source is
 *                                      present (i.e. Aperture = RFPT)
 *              08/21/03    1.6  wvd    Change channel to unsigned char.
 *              09/08/03    1.7  wvd    Fix typo in verbose statement.
 *              09/18/03    1.8  wvd    Stop iteration if
 *					cf_calculate_y_centroid returns -1.
 *              10/22/03    1.9  wvd    Rewrite program.
 *              11/12/03    1.10 wvd    In HIST mode, use default centroid
 *					for non-target channels.
 *              04/09/04    1.11 bjg    Include stdlib.h and stdio.h
 *                                      Change formats in printf to 
 *                                      match arg types.
 *              04/19/04    1.12 wvd    Sum background only over detector
 *					active area.  Scale intrinsic count
 *					rate by EXPTIME and detector width.
 *              05/03/04    1.13 wvd    Move call to cf_identify_channel
 *					from this routine to cf_remove_motions.
 *              06/02/04    1.14 wvd    Don't need PHA array; can use locflags
 *					instead.  Use time flags to exclude
 *					bursts, etc.
 *              04/21/05    1.15 wvd    Rewrite program using airglow lines
 *					(if available) to constrain search
 *					for target centroid.  Don't bother
 *					reading model background files.
 *              06/03/05    1.16 wvd    Fix bug in tabulation of exected Y
 *					centroids. Remove unused variables.
 *              02/01/06    1.17 wvd    Tighten region over which centroids
 *					are computed.
 *					If target centroid is too far from
 *					airglow centroid, reject it.
 *					If ycent is too far from default value,
 *					use default, not airglow value.
 *					Change ysum and norm to doubles.
 *              08/30/06    1.18 wvd    Added airglow_centroid argument.
 *		03/14/08    1.19 wvd	In HIST mode, there are no spectra in
 *					other apertures to cause confusion,
 *					so we need not require that the 
 *					target centroid be within 30 pixels
 *					of the airglow value.
 *		11/25/08    1.20 wvd	Don't compare HIST centroids to
 *					either airglow or tabulated centroid.
 *
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "calfuse.h"

#define YLENGTH		NYMAX/16

int
cf_find_spectra(fitsfile *header, long nevents,
	float *weight, float *xfarf, float *yfarf, unsigned char *channel,
	unsigned char *timeflags, unsigned char *locflags, int airglow_centroid)
{
    char aper[FLEN_VALUE], instmode[FLEN_VALUE], quality[FLEN_VALUE];
    char key[FLEN_KEYWORD], filename[FLEN_VALUE], detector[FLEN_VALUE];
    double norm, ysum;
    float bkgd, peak, ycent, yshift=0;
    float yairg[NYMAX], ydist[NYMAX], ysigma[NYMAX];
    float ybin[YLENGTH], sbin[YLENGTH], ycent_tab[8];
    float dy, sigma, xmin, xmax, ycent_air;
    int  active_ap[2], airglow=FALSE, target_ap, status=0;
    int  bkgd_num, bkgd_min[8], bkgd_max[8];
    int  chan_num, errflg = 0, lwrs;
    int  chan_order[] = { 3, 2, 1, 7, 6, 5 };
    int  center, high, low, npeak, hdu;
    int  spex_sic, spex_lif, emax, emax_sic, emax_lif, extended;
    long i, j, k;
    fitsfile  *chidfits, *parmfits;

    char CF_PRGM_ID[] = "cf_find_spectra";
    char CF_VER_NUM[] = "1.20";

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Check whether routine is appropriate for this data file. */
    if ( (errflg = cf_proc_check(header, CF_PRGM_ID)) ) return errflg;

    /* Initialize Y-distribution arrays. */
    for (i = 0; i < NYMAX; i++)
	yairg[i] = ydist[i] = ysigma[i] = 0.;
    for (i = 0; i < YLENGTH; i++)
	ybin[i] = sbin[i] = 0.;

    /* Read header keywords from IDF. */
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(header, TSTRING, "DETECTOR", detector, NULL, &status);
    FITS_read_key(header, TSTRING, "APERTURE", aper, NULL, &status);
    extended = cf_source_aper(header, active_ap);

    FITS_read_key(header, TINT, "BKGD_NUM", &bkgd_num, NULL, &status);
    for (i = 0; i < bkgd_num; i++) {
	sprintf(key, "BKG_MIN%ld", i);
        FITS_read_key(header, TINT, key, bkgd_min+i, NULL, &status);
	sprintf(key, "BKG_MAX%ld", i);
        FITS_read_key(header, TINT, key, bkgd_max+i, NULL, &status);
    }

    /* Read expected channel centroids from CHID_CAL file. */
    FITS_read_key(header, TSTRING, "CHID_CAL", filename, NULL, &status);
    FITS_open_file(&chidfits, cf_cal_file(filename), READONLY, &status);
    for (i = 1; i < 8; i++) {
	if (i==4) continue;
        if ((i == active_ap[0] || i == active_ap[1]) && !extended) hdu=i+1;
	else hdu=i+9;
	FITS_movabs_hdu(chidfits, hdu, NULL, &status);
	FITS_read_key(chidfits, TFLOAT, "CENTROID", ycent_tab+i, NULL, &status);
    }
    FITS_read_key(chidfits, TFLOAT, "XMIN", &xmin, NULL, &status);
    FITS_read_key(chidfits, TFLOAT, "XMAX", &xmax, NULL, &status);
    FITS_close_file(chidfits, &status);

    /* 
     * Generate photon and airglow arrays as a function of YFARF.
     * Ignore events near the detector edge.
     */
    xmin += 250;
    xmax -= 250;
    for (i = 0; i < nevents; i++) {
	if (!(timeflags[i] & ~TEMPORAL_DAY) &&
	    (xfarf[i] > xmin) && (xfarf[i] < xmax)) {
	    if (!locflags[i])
	        ydist[cf_nint(yfarf[i])] += weight[i];
	    else if (locflags[i] == LOCATION_AIR)
	        yairg[cf_nint(yfarf[i])] += weight[i];
	}
    }
    /* Set error array equal to photon counts (= variance). */
    for (i = 0; i < NYMAX; i++)
	ysigma[i] = ydist[i];

    /* If exposure was taken in TTAG mode, estimate BKGD level */
    bkgd = 0.;
    if (!strcmp(instmode, "TTAG")) {
	k = 0;
	for (i = 0; i < bkgd_num; i++) {
	    for (j = bkgd_min[i]; j <= bkgd_max[i]; j++) {
		bkgd += ydist[j];
		k++;
	    }
	}
	bkgd /= k;
	cf_verbose(3, "Mean background level: %.1f", bkgd);
    }

    /* Subtract mean background and bin data by 16 pixels. */
    for (i = 0; i < NYMAX; i++) {
	j = i / 16;
	ybin[j] += ydist[i] - bkgd;
	sbin[j] += ysigma[i];
    }
    for (i = 0; i < YLENGTH; i++) sbin[i] = sqrt(sbin[i]);

    /* Read extraction parameters from PARM_CAL file. */
    FITS_read_key(header, TSTRING, "PARM_CAL", filename, NULL, &status);
    FITS_open_file(&parmfits, cf_parm_file(filename), READONLY, &status);
    FITS_read_key(parmfits,TINT,"SPEX_SIC",&spex_sic,NULL,&status);
    FITS_read_key(parmfits,TINT,"SPEX_LIF",&spex_lif,NULL,&status);
    FITS_read_key(parmfits,TINT,"EMAX_SIC",&emax_sic,NULL,&status);
    FITS_read_key(parmfits,TINT,"EMAX_LIF",&emax_lif,NULL,&status);
    FITS_close_file(parmfits,&status);

    /* Determine centroid for each channel. */
    for (k = 0; k < 6; k++) {
	chan_num = chan_order[k];

	/* Is this the LWRS aperture? */
	if (chan_num == 3 || chan_num == 7)
	    lwrs = TRUE;
	else 
	    lwrs = FALSE;

	/* Is it the target aperture? */
	if (chan_num == active_ap[0] || chan_num == active_ap[1])
	    target_ap = TRUE;
	else
	    target_ap = FALSE;

        /* If we're in histogram mode, and this is not the target aperture,
	   move to next aperture. */
        if (!strncmp(instmode, "HIST", 4) && !target_ap) {
	    if (lwrs) yshift = 0;
	    continue;
	}

	/* First, we compute a centroid from the airglow lines.
	   For the LWRS aperture, search within 70 pixels of
	   the expected centroid.  Compute the offset between the
	   measured and tabulated centroids.  For the MDRS and HIRS
	   apertures, apply the offset computed from the LWRS aperture
	   to the tabulated centroid.
	   Disregard regions near the bottom and top of the detector. */

	if (lwrs) {
	    center = cf_nint(ycent_tab[chan_num]);
	    dy = 70;
            if ((low  = center-dy) < 16) low = 16;
            if ((high = center+dy) > NYMAX-50) high = NYMAX-50;
	    norm = ysum = 0.;
	    for (i = low; i <= high; i++) {
		ysum += i * yairg[i];
		norm += yairg[i];
	    }

	    /* If the airglow lines contain at least 33 events, use their
	       centroid.  If not, use the tabulated centroid. */
	    if (norm > 33.) {
		ycent = ysum / norm;
		yshift = ycent - ycent_tab[chan_num];
		airglow = TRUE;
	    }
	    else {
		ycent = ycent_tab[chan_num];
		yshift = 0;
		airglow = FALSE;
	    }
	}
 		/* For MDRS and HIRS apertures */
	else
	    ycent = yshift + ycent_tab[chan_num];

	/* Remember airglow centroid */
	ycent_air = ycent;

	/* Now we find the centroid for the target spectrum.  If we are
	   in a target aperture and user has specified YCENT, use it. */

	if (target_ap && ((chan_num<5 && spex_lif>=0 && spex_lif<1024) ||
	    (chan_num>4 && spex_sic>=0 && spex_sic<1024))) {
            if (chan_num<5) ycent=spex_lif;
            else ycent=spex_sic;
            sprintf(key, "YCENT%1d", chan_num);
            FITS_update_key(header, TFLOAT, key, &ycent, NULL, &status);
	    sprintf(quality, "HIGH");
            sprintf(key, "YQUAL%1d", chan_num);
            FITS_update_key(header, TSTRING, key, quality, NULL, &status);
            cf_verbose(1, "Channel %d: User-specified Y centroid = %.2f",
                chan_num, ycent);
	    continue;
        }

	/* If we must find the target spectrum ourselves, use the binned
	   data array to improve the S/N.  Search within 2 binned pixels
	   of ycent, however it was computed. */

	center = cf_nint(ycent / 16.);
	dy = 2;
	if ((low = center - dy) < 3) low = 3;
	if ((high = center + dy) > YLENGTH-4) high = YLENGTH-4;

	npeak = 0;
	peak = -1.0E5;
	for (i = low; i <= high; i++) {
	    if (ybin[i] > peak) {
		peak = ybin[i];
		npeak = i;
	    }
	}

	/* If peak is significant, compute Y centroid for events falling
	   within 40 pixels of peak.  If this number is within 30 pixels
	   of airglow centroid, use it.  If not, use airglow centroid.
	   Note that we require a higher significance for the SiC LWRS
	   channel on side 1, because it sits on an elevated background. */

	if (!strncmp(detector, "1", 1) && chan_num == 7) sigma = 9.;
	else sigma = 5.;

	sprintf(quality, "UNKNOWN");

	if (peak > sigma * sbin[npeak] && !airglow_centroid) {
	    center = npeak * 16 + 8;
	    dy = 40;
	    if ((low = center - dy) < 32) low = 32;
	    if ((high = center + dy) > NYMAX-50) high = NYMAX-50;
	    norm = ysum = 0.;
	    for (i = low; i <= high; i++) {
		ysum += i * ydist[i];
		norm += ydist[i];
	    }
	    ycent = ysum / norm;

	    /* Don't test offset for HIST data. (wvd, 03/14/2008) */
	    if (fabs(ycent - ycent_air) < 30 || !strncmp(instmode, "HIST", 4)) {
	    	sprintf(quality, "HIGH");
	    	cf_verbose(3, "Channel %d: default centroid: %.2f",
			chan_num, ycent_tab[chan_num]);
	    	if (airglow) cf_verbose(3, "Channel %d: airglow centroid: %.2f",
			chan_num, ycent_air);
	    	cf_verbose(3,"Channel %d: counts in peak = %.0f, limit = %.0f",
			chan_num, peak, sigma * sbin[npeak]);
	    	cf_verbose(2, "Channel %d: using target centroid = %.2lf",
			chan_num, ycent);
	    }
	}
	/* If peak is not significant or ycent differs too much from
	   the airglow centroid, use the airglow centroid. */
	if (!strncmp(quality, "UNKNOWN", 1)) {
	    if (airglow) {
		sprintf(quality, "MEDIUM");
	        cf_verbose(3, "Channel %d: default centroid: %.2f",
		    chan_num, ycent_tab[chan_num]);
		cf_verbose(3,"Channel %d: peak = %.0f, limit = %.0f",
		    chan_num, peak, sigma * sbin[npeak]);
		cf_verbose(2, "Channel %d: using airglow centroid = %.2lf",
		    chan_num, ycent);
	    }
	    else {
	        sprintf(quality, "LOW");
		cf_verbose(3,"Channel %d: peak = %f, limit = %f",
		    chan_num, peak, sigma * sbin[npeak]);
		cf_verbose(2, "Channel %d: using default centroid = %.2lf",
		    chan_num, ycent);
	    }
	}

	/* If centroid lies too far from the default value, use default. */
        if (chan_num<5) emax=emax_lif; else emax=emax_sic;
	    /* Don't test offset for HIST data. (wvd, 11/25/2008) */
        if (fabs(ycent - ycent_tab[chan_num]) > emax && !strncmp(instmode, "TTAG", 4)) {
	    if (target_ap) {
        	cf_verbose(1, "Computed centroid (%.1lf) is too far from "
		    "expected value.", ycent);
	        cf_verbose(1, "Channel %d: using default centroid of %.1f",
		    chan_num, ycent_tab[chan_num]);
	    }
	    else {
	        cf_verbose(2, "Channel %d: using default centroid of %.1f",
		    chan_num, ycent_tab[chan_num]);
	    }
            ycent = ycent_tab[chan_num];
	    sprintf(quality, "LOW");
        }

	/* Write results to file header */
        sprintf(key, "YCENT%1d", chan_num);
        FITS_update_key(header, TFLOAT, key, &ycent, NULL, &status);
        sprintf(key, "YQUAL%1d", chan_num);
        FITS_update_key(header, TSTRING, key, quality, NULL, &status);
    }

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");

    return status;
}
