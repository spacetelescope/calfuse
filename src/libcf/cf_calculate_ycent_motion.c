/**************************************************************************
 *           Johns Hopkins University
 *          Center for Astrophysical Sciences
 *           FUSE
 *************************************************************************
 *
 * synopsis:    cf_calculate_ycent_motion(header, nevents, ptime, y, channel,
 *              locflags, nsec, ttime, ycentl, ycents)
 *
 * Description: Determines the Y centroid of the emission as a function
 *              of time within the target apertures
 *
 * Arguments:   fitsfile  *header :   pointer to Intermediate Data File
 *              long      nevents :   number of photon events in the file
 *              float     *ptime  :   detection time for each photon
 *              float     *y      :   y position of each photon
 *       unsigned char    *channel:   aperture associated with each photon
 *       unsigned char   *locflags:   location flag array
 *              long      nsec    :   number of seconds tabulated in timeline
 *              float     *ttime  :   tabulated times in the timeline
 *                                    geocoronal photons
 *              float     *ycentl, *ycents : 
 *                                    y centroids of the Lif and SiC apertures,
 *                                    tabulated once per second throughout
 *                                    the observation
 *             
 *
 * Calibration files required:     None
 *
 * Returns:	0 on success
 *
 *                                
 * HISTORY:     11/05/02   v1.1   RDR   started work
 *              11/27/07   v1.2   RDR   corrected a bug in determining ycent
 *                                        cleaned up variables
 *              12/04/02   v1.3   RDR   locflags out geocoronal emission
 *                                        before determining centroids
 *              12/12/02   v1.4   wvd   change order of subroutine arguments
 *		01/17/03   v1.5   wvd   call cf_update_proc()
 *              01/20/03   v1.6   rdr   correct error in calculating ycent
 *                                        during the last time interval
 *              02/24/03   v1.7   peb   Changed include file to calfitsio.h
 *		03/11/03   v1.8   wvd   Changed locflags to unsigned char
 *              05/20/03   v1.9   rdr   Added call to cf_proc_check
 *		05/22/03   v1.10  wvd   cf_error_init to stderr
 *		06/04/03   v1.11  wvd   Implement cf_verbose throughout.
 *		08/21/03   v1.12  wvd   Change channel to unsigned char.
 *		10/06/03   v1.13  wvd   Change screen to locflags throughout.
 *              10/21/04   v1.14  bjg   Corrected several bugs
 *              04/07/07   v1.15  wvd	Clean up compiler warnings.
 *              04/07/07   v1.16  wvd	Clean up compiler warnings.
 *
 **************************************************************************/

#include "calfuse.h"


int
cf_calculate_ycent_motion(fitsfile *header, long nevents, float *ptime, 
			  float *y, unsigned char *channel, unsigned char *locflags,
			  long nsec, float *ttime, float *ycentl, float *ycents)
{
    char CF_PRGM_ID[] = "cf_calculate_ycent_motion";
    char CF_VER_NUM[] = "1.16";
 
    char ycentname[FLEN_CARD];
    int chan_num, status=0, j, src_type, active_ap[2], nave=500;
    int errflg=0;
    long i, ndx, ndxs, nsam;
    float ycent, ptst, dt, dt1;
    float ycentdef[2];
    float *ycentptr[2];

    ycentptr[0]=ycentl;
    ycentptr[1]=ycents;

    /* Initialize error checking */ 
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Check whether routine is appropriate for this data file. */
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* nave = number of photons to average in determining the centroid.
       Should eventually be read from a parameter file. */
    nave = 500;

    /* Determine the source aperture and type from the header data.
       src_type = 0 for a point source and 1 for an extended source. */

    src_type = cf_source_aper(header, active_ap);
    cf_verbose(3, "active_ap = %d and %d,   src_type = %d",
	   active_ap[0], active_ap[1], src_type);


    /* Determine the centroids.  Do LiF and SiC apertures separately. */

    for (j=0; j<2; j++) {

 
      /* Select the channel number of the appropriate active aperture */
      chan_num = active_ap[j];
      sprintf(ycentname, "YCENT%1d", chan_num);
      FITS_read_key(header, TFLOAT, ycentname, &ycentdef[j], NULL, &status);

      nsam = 0;
      ycent = 0.;
      ptst = ptime[0];
      ndxs = 0;

      /* Go through all of the photons, determine the y centroids,
	 and fill in the arrays.  Avoid airglow lines. */ 
      for (i=1; i<nevents; i++) {
	dt = ptime[i] - ptst;
	dt1 = ptime[i] - ptime[i-1];
	if (channel[i] == chan_num && !(locflags[i] & LOCATION_AIR)) {
	  ycent += y[i];
	  nsam += 1;
	}
	if ((nsam > nave && dt > 1.) || i == nevents-1 || dt1 > 10)  {
	  
	  if (nsam==0) ycent = ycentdef[j];
	  else ycent = ycent / nsam;
	  
	  ndx = ndxs;
	  
	  while ( ndx < nsec && (ndx==0 || ttime[ndx-1] < ptime[i])) {
	    
	    if (nsam > nave) ycentptr[j][ndx] = ycent;
	    else {
	      if (ndx==0) ycentptr[j][ndx] = ycentdef[j];
	      else ycentptr[j][ndx] = ycentptr[j][ndx-1]; 
	    }
	    
	    
	    ndx += 1; 
	  }
	  
	  ndxs=ndx;
	  if (ndxs > nsec-1) ndxs = nsec-1; 
	  nsam = 0;
	  ycent=0.;
	  ptst=ptime[i];
	}
      }
    }
    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished processing"); 
    return (status);
}
