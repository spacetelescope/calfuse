/**************************************************************************
 *           Johns Hopkins University
 *           Center for Astrophysical Sciences
 *           FUSE
 *************************************************************************
 *
 * Synopsis:    cf_target_count_rate(header, nevents, ptime, weight,
 *              channel, locflags, ntimes, ttime, rate_lif, rate_sic)
 *
 * Description: Determines the count rate in the target aperture for both
 *		LiF and SiC channels.  Excludes airglow lines.
 *
 * Arguments:   fitsfile  *header :   pointer to IDF header
 *              long      nevents :   number of photon events in the file
 *              float     *ptime  :   detection time for each photon
 *              float     *weight :   weight associated with each photon
 *     unsigned char      *channel:   aperture associated with each photon
 *     unsigned char     *locflags:   location flag array - flags the 
 *                                    geocoronal photons
 *              long      ntimes    :   number of seconds tabulated in timeline
 *              float     *ttime  :   tabulated times in the timeline
 *              float	*rate_lif, *rate_sic :
 *                                    count rates through the Lif and SiC
 *				      apertures, excluding geocoronal emission
 *
 * Calibration files required:     None
 *
 * Returns:	0 on success
 *
 *                                
 * HISTORY:     03/03/03    v1.1   RDR    started work
 *		03/05/03    v1.2   wvd    change name of subroutine & install
 *              03/10/03    v1.21  rdr    changed specification of locflags
 *                                        from char to unsigned char
 *              05/20/03    v1.3   rdr    Added call to cf_proc_check
 *              05/22/03    v1.5   wvd    Direct cf_error_init to stderr
 *              06/02/03    v1.6   wvd    Implement cf_verbose throughout.
 *              06/04/03    v1.7   wvd    Revise scheme for summing counts.
 *              08/20/03    v1.8   wvd    Use FRAME_TOLERANCE from calfuse.h,
 *					  change channel to unsigned char.
 *              09/15/03    v1.9   wvd    Test that photon times match
 *					  timeline times to within 1 sec.
 *              10/06/03    v1.10  wvd    Changed screen to locflags throughout.
 *              12/29/06    v1.11  wvd    Sum weights rather than counts.
 *					  Convert output arrays to type float.
 *					  Scale HIST count rate by DET_DEAD.
 *              04/07/07    v1.12  wvd	  Clean up compiler warnings.
 *              04/07/07    v1.13  wvd	  Clean up compiler warnings.
 *
 **************************************************************************/

#include <string.h>
#include "calfuse.h"

int
cf_target_count_rate(fitsfile *header, long nevents, float *ptime, 
	float *weight, unsigned char *channel, unsigned char *locflags,
	long ntimes, float *ttime, float *rate_lif, float *rate_sic) {

    char CF_PRGM_ID[] = "cf_target_count_rate" ;
    char CF_VER_NUM[] = "1.13" ;
 
    char instmode[FLEN_VALUE];
    int errflg=0, status=0, active_ap[2] ;
    float delta_max = 1.0 + FRAME_TOLERANCE;
    float det_dead;
    long j, k ;
 
    /* Initialize error checking. */ 
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Confirm that this subroutine should be run. */
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Determine the source aperture from the header data. */
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    (void) cf_source_aper(header, active_ap) ;

    /*
     * For histogram data, scale the LiF and SiC count rates, which come
     * from the housekeeping data file, by DET_DEAD.
     */
    if (!strncmp(instmode, "HIST", 4)) {
	FITS_read_key(header, TFLOAT, "DET_DEAD", &det_dead, NULL, &status);
	for (k = 0; k < ntimes; k++) rate_lif[k] *= det_dead;
	for (k = 0; k < ntimes; k++) rate_sic[k] *= det_dead;
    }

    /* 
     * For time-tag data, count the photons arriving during each second
     * of the exposure.  Exclude events near airglow regions.
     */
    else {
	for (k = 0; k < ntimes; k++) rate_lif[k] = rate_sic[k] = 0.; 
	for (j=k=0; j<nevents; j++) {
            while(ttime[k+1]-FRAME_TOLERANCE < ptime[j] && k+1 < ntimes)
		k++;
            if (!(locflags[j] & LOCATION_AIR) &&
		(ptime[j] - ttime[k] < delta_max)) {
		if (channel[j] == active_ap[0]) rate_lif[k] += weight[j] ;
		if (channel[j] == active_ap[1]) rate_sic[k] += weight[j] ;
	    }
	}
    }

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing"); 
    return (status);
}
