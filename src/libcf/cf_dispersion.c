/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Description: Assigns wavelength to each photon in data stream.  
 *                    
 * Returns:     0 upon successful completion.
 *
 * History:     05/15/06   1.1   wvd	Adapt from cf_astigmatism_and_dispersion.c
 *					Program incorporates cf_x2lambda.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_dispersion";
static char CF_VER_NUM[] = "1.1";


int
cf_dispersion(fitsfile *infits, long nevents, float *x,
		   unsigned char *channel, float *lambda)
{
    char  	wave_file[FLEN_VALUE];
    float	*wavelength=NULL;
    float	x_lower, x_upper, weight_l, weight_u, wave_l, wave_u;
    int	  	ap, errflg=0, nwave, xl_int, xu_int, status=0;
    long	k;
    fitsfile 	*wavefits;

    /* Enter a timestamp into the log. */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    if ((errflg = cf_proc_check(infits, CF_PRGM_ID))) return errflg;

    FITS_read_key(infits, TSTRING, "WAVE_CAL", wave_file, NULL, &status);
    FITS_open_file(&wavefits, cf_cal_file(wave_file), READONLY, &status);

    /* Cycle through channels, skipping #4 (pinhole) */
    for (ap = 1; ap < 8; ap++) {
	if (ap == 4) continue;

    	FITS_movabs_hdu(wavefits, ap+1, NULL, &status);
	nwave = cf_read_col(wavefits, TFLOAT, "WAVELENGTH", (void **) &wavelength);

	/* Compute wavelength for each event by interpolation. */
        for (k = 0; k < nevents; k++) {
	    if (channel[k] == ap) {
		x_lower = floor(x[k]);
		x_upper = ceil(x[k]);
		xl_int	= (int) (x_lower + 0.5);
		xu_int	= (int) (x_upper + 0.5);
		if (xl_int >= 0 && xu_int < nwave) {
		    if (xl_int == xu_int)
			lambda[k] = wavelength[xl_int];
		    else {
		        weight_l = x_upper - x[k];
		        weight_u = x[k] - x_lower;
		        wave_l = wavelength[xl_int];
		        wave_u = wavelength[xu_int];
		        lambda[k] = wave_l * weight_l + wave_u * weight_u;
		    }
		}
		else
		    channel[k] = 0;
	    }
	}
        /* Space for wavelength array is allocated in each loop. */
	free(wavelength);
    }
    FITS_close_file(wavefits, &status);

    /* Update processing flags. */
    cf_proc_update(infits, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");

    return status;
}
