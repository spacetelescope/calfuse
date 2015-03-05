/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_doppler_and_heliocentric(*infits, nevents, *photon_time, 
 *		  	*channel, *lambda, nseconds, *timeline_time, 
 *			*timeline_velocity);
 *
 * Description: Doppler-correct each photon for both heliocentric and orbital
 *		motions.  Resulting wavelength scale is heliocentric.
 *
 * Arguments:	*infits			Input FITS file pointer
 * 		nevents			number of points in the photon list
 * 		*photon_time		time array of photon events
 * 		*channel		channel number in the photon list
 * 		*lambda			wavelength in the photon list
 * 		nseconds		number of points in the timeline table
 * 		*timeline_time		time array of timeline list
 * 		*timeline_velocity	orbital velocity in km/s
 *
 * Returns:	O upon successful completion
 *
 * History:     12/09/2002   jch   1.1	Initial coding
 *	        12/11/2002   wvd   1.2	Move Doppler calculation to separate
 *					  loop.
 *              12/11/2002   wvd   1.3  Change nevents and nseconds to long.
 *              12/20/2002   wvd   1.4  Change channel to unsigned char.
 *              02/24/03     peb   1.5  Change include file to calfusettag.h
 *                                      and calfitsio.h
 *              03/11/2003   wvd   1.6  Change channel to char.
 *              05/20/2003   rdr   1.7  Added call to cf_proc_check
 *              09/17/2003   wvd   1.8  Return errflg from cf_prock_check
 *              10/21/2003   wvd   1.9  Change channel to unsigned char.
 *              04/21/2004   bjg   1.10 Cosmetic change to prevent warning
 *                                      with gcc -Wall
 *              11/17/2005   wvd   1.11 Change velocity and doppler to doubles.
 *
 ****************************************************************************/

#include <stdlib.h>
#include "calfuse.h"

char CF_PRGM_ID[] = "cf_doppler_and_heliocentric";
char CF_VER_NUM[] = "1.11";


int cf_doppler_and_heliocentric(fitsfile *infits, long nevents,
	float *photon_time, unsigned char *channel, float *lambda,
	long nseconds, float *timeline_time, float *timeline_velocity)
{
	int	errflg=0, status=0;
	long	i, k;
	float	v_helio;	/* heliocentric velocity in km/sec */
	double	velocity, *doppler;

	/* Enter a timestamp into the log. */
    	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
    	cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

        if ((errflg = cf_proc_check(infits, CF_PRGM_ID) )) return errflg;

	FITS_read_key(infits, TFLOAT, "V_HELIO", &v_helio, NULL, &status);

	/* Compute array of Doppler corrections. */
	doppler = (double *) cf_malloc(sizeof(double) * nseconds);
	for (k = 0; k < nseconds; k++) {
		velocity = timeline_velocity[k] + v_helio;
		doppler[k] = 1. + velocity / C;
	}

	/* Apply Doppler correction to each photon assigned to an aperture. */
	k = 0;
	for (i = 0; i < nevents; i++)
	    if (channel[i] > 0) {
	        while (photon_time[i] > timeline_time[k+1] - FRAME_TOLERANCE
		    && k < nseconds-1) k++;
	        lambda[i] *= doppler[k];
	    }
	free (doppler);

	/* Update processing flags. */
    	cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");
    	return (status);
}
