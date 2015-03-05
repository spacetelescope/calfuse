/*****************************************************************************
 *	        Johns Hopkins University
 *	        Center For Astrophysical Sciences
 *	        FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_mirror_motion(header, nevents, time, x, y, channel,
 *			ntimes, ttime, tsunset)
 *
 * Description: Calculates the shift of the spectrum in both X and Y caused
 *              by thermal changes in the mirror position and corrects
 *              the X and Y coordinates of each photon event. 
 *
 *              fitsfile  *header       Pointer to the location of the FITS
 *                                      header of the Intermediate data File
 *              long      nevents       Number of photons in the file
 *              int       *time         Time stamp (in seconds) since the
 *                                      start of the exposure
 *              float     *x            Position of the photon (in pixels)
 *              float     *y            Position of the photon (in pixels)
 *     unsigned char      channel       Channel id of each photon
 *              long      ntimes        Number of entries in timeline table
 *              float     ttime         Time array of timeline table
 *              short     tsunset       Time since last sunset
 *
 * Returns:     0 on success
 *
 * History:	09/06/02   1.1  RDR     Begin work, adapted from cf_make_shift
 *		03/01/03   1.3  wvd	Correct use of pointer in FITS_read_key()
 *              04/21/03   1.4  wvd	Use tsunset array from timeline table.
 *					Do not assume that ttime is continuous.
 *					Interpolate between tabulated shifts.
 *              05/20/03   1.5  rdr     Added call to cf_proc_check
 *              08/21/03   1.6  wvd     Change channel to unsigned char.
 *              06/22/04   1.7  wvd     Estimate time between sunsets using
 *					orbit period in file header.
 *              07/14/04   1.8  rdr     Correct line which selects calibration
 *              03/22/05   1.9  wvd     Change tsunset from float to short.
 *					Read orbital period from file header.
 *              04/15/05   1.10 wvd     If ORBPERID keyword is not present,
 *					assume that it is 6000 s.
 *              12/20/05   1.11 wvd     Omit correction for an extended
 *					source.
 *              04/07/07   1.12 wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

int
cf_mirror_motion(fitsfile *header, long nevents, float *time, float *x,
		 float *y, unsigned char *channel, long ntimes, float *ttime,
                 short *tsunset) {

    char CF_PRGM_ID[] = "cf_mirror_motion";
    char CF_VER_NUM[] = "1.12";

    fitsfile *mmfits;
    int    errflg=0, status=0, anynull=0, hdutype, nullval=0;
    int    active_ap[2], extended;
    long   j, k, ndx;
    char   det[FLEN_VALUE], mmcal[FLEN_VALUE];
    float  *tdxlif, *tdxsic, *dxlif, *dylif, *dxsic, *dysic;
    float  frac, period, w0, w1, w99;

    /* Initialize error checking. */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Check whether routine is appropriate for this data file. */
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Read detector keyword from header. */
    FITS_read_key(header, TSTRING, "DETECTOR", det, NULL, &status);

    /* Determine the channel numbers of the active apertures */
    extended = cf_source_aper(header, active_ap);
    
    /* If source is extended, exit now. */
    if (extended) {
	cf_verbose(1, "Extended source.  Omitting photon shift.");
        cf_proc_update(header, CF_PRGM_ID, "SKIPPED");
	return (status);
    }

    /* 
     *  Read orbital period from file header.
     */
    fits_read_key(header, TFLOAT, "ORBPERID", &period, NULL, &status);
    if (status) {
	status = 0;
	period = 6000;
        cf_verbose(1, "Keyword ORBPERID not found; assuming 6000 seconds");
    }
    else
	cf_verbose(2, "Estimated orbital period is %.2f seconds", period);

    /* Allocate space for shift arrays. */
    dxlif = (float *) cf_calloc(ntimes, sizeof(float));
    dylif = (float *) cf_calloc(ntimes, sizeof(float));
    dxsic = (float *) cf_calloc(ntimes, sizeof(float));
    dysic = (float *) cf_calloc(ntimes, sizeof(float));

    /* Read the name of the mirror-motion calibration file */
    FITS_read_key(header, TSTRING, "MIRR_CAL", mmcal, NULL, &status); 
    cf_verbose(3, "Mirror motion calibration file = %s", mmcal);

    /*  Open the calibration file and read in the information relating 
     *  x shifts as a function of time (in minutes) from sunset for one orbit
     *  (100 minutes)
     */
    FITS_open_file(&mmfits, cf_cal_file(mmcal), READONLY, &status);
    tdxsic = (float *) cf_calloc(100, sizeof(float)); 
    tdxlif = (float *) cf_calloc(100, sizeof(float)); 
    cf_verbose(3,"detector= %s",det) ;
    if (!strncmp(det, "1", 1) ) {
      cf_verbose(3,"Correcting detector 1") ;
	/*
	 *  NOTE: LiF 1 is used as reference and so, by definition, has no
	 *  mirror motions. Thus, on side 1 only the SiC spectrum will move.
	 */
	FITS_movabs_hdu(mmfits,2, &hdutype, &status);
	FITS_read_img(mmfits, TFLOAT, 1,100, &nullval, tdxsic, &anynull,
		      &status);
    }
    else {
      cf_verbose(3,"Correcting detector 2") ;
	FITS_movabs_hdu(mmfits,3, &hdutype, &status);
	FITS_read_img(mmfits, TFLOAT, 1,100, &nullval, tdxlif, &anynull,
		      &status);
	FITS_movabs_hdu(mmfits,4, &hdutype, &status);
	FITS_read_img(mmfits, TFLOAT, 1,100, &nullval, tdxsic, &anynull,
		      &status);
    }
    FITS_close_file(mmfits, &status);
    /* Determine the shifts for each second of the observation. */
    for (k=0; k<ntimes; k++) {
	frac = tsunset[k] / period * 100.;
	ndx = (long) frac;
	if (ndx >= 99) {
	    w99 = 100. - frac;
	    w0 = frac - 99.;
	    dxlif[k] = w0 * tdxlif[0] + w99 * tdxlif[99];
	    dxsic[k] = w0 * tdxsic[0] + w99 * tdxsic[99];
	}
	else {
	    w0 = (float) (ndx + 1) - frac;
	    w1 = frac - (float) (ndx);
	    dxlif[k] = w0 * tdxlif[ndx] + w1 * tdxlif[ndx + 1];
	    dxsic[k] = w0 * tdxsic[ndx] + w1 * tdxsic[ndx + 1];
	}
    }
    
    /* Apply the calculated shifts to the data. */
    for (j=k=0; j<nevents; j++) {
        while(ttime[k+1] - FRAME_TOLERANCE < time[j] && k+1 < ntimes) k++;
	/*
	 *  Shift only photons in the active aperture.
	 *  active_ap[0] = lif, active_ap[1] = sic 
	 */
	if (channel[j] == active_ap[0]) {
            x[j] += dxlif[k];
            y[j] += dylif[k];

        }
	if (channel[j] == active_ap[1]) {
            x[j] += dxsic[k];
            y[j] += dysic[k];
        }
    }

    /* Release array storage */
    free(tdxlif);
    free(tdxsic);
    free(dxlif);
    free(dylif);
    free(dxsic);
    free(dysic);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");

    return status;
}
