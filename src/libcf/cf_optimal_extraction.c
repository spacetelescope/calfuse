/*****************************************************************************
 *              Johns Hopkins University 
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_optimal_extraction(fitsfile *header, int optimal, 
 *              int aperture, float *weight, float *y, unsigned char *channel,
 *              float *lambda, long ngood, long *good_index,
 *		float *barray, float *bpmask, int pny, float pycent,
 *		float *parray, long nout, float *wave_out, float **flux_out,
 *		float **sigma_out, long **counts_out, float **weights_out, 
 *		float **bkgd_out, short **bpix_out)
 *
 * Description: Perform standard or optimal extraction of target spectrum.
 *
 * Note:	We use a modified version of the optimal-extraction
 *		algorithm described by Keith Horne (PASP, 98, 609, 1986).
 *
 * Arguments:   fitsfile  *header       Pointer to IDF FITS file header
 *		int	  optimal	TRUE if optimal extraction requested.
 *		int	  aperture	Target aperture (values 1-8)
 *		float	  *weight       Scale factor for each photon
 *		float	  *y            Final y position of each photon
 *		unsigned char *channel	Channel ID of each photon
 *		float 	  *lambda       Wavelength for each photon
 *		long	  ngood		Length of array good_index 
 *		long	  *good_index	Indices of screened photon events
 *		float     *barray       2-D background array
 *                                       (flux-calibrated and binned)
 *              float     *bpmask       2-D array containing bad pixel locations
 *                                      Same size and location as the background array
 *	        int       pny           Size in Y of weights array
 *	        float     pycent        Y centroid of weights array
 *	        float     parray        2-D weights array (binned)
 *	        long      nout          Length of output arrays
 *	        float     *wave_out	Output wavelength array
 *	        float     **flux_out    Output flux array
 *	        float     **sigma_out   Output error array
 *	        long      **counts_out  Output counts array
 *	        float     **weights_out Output weights array
 *	        float     **bkgd_out 	Output background array
 *              short     **bpix_out    Output bad pixel spectrum
 *
 * Returns:     0 on success
 *
 * History:     01/28/03   1.1   wvd    Initial coding
 *              02/28/03   1.2   peb    Changed flux_out, var_out, counts_out
 *                                      and weight_out argument variables to
 *                                      pointer-to-pointers, so the arrays can
 *                                      be returned to the calling routine.
 *                                      Also tidied the code using lint.
 *		03/11/03   1.3   wvd    Update documentation regarding use of
 *					unsigned char
 *		03/12/03   1.4   wvd    Change weights_out to sum of weights
 *		03/19/03   1.5   wvd    Divide flux_out and var_out by WPC
 *              04/04/03   1.6   rdr    - correct bugs in optimum extraction
 *                                      - adjust criterion for redetermining
 *                                        y centroid
 *		05/07/03   1.7   wvd	Change ergcm2s to ergcm2 throughout. 
 *					Divide by EXPTIME.  Use cf_verbose.
 *		05/22/03   1.8   wvd	If EXPTIME = 0, set flux and error to 0.
 *              05/28/03   1.9   rdr    - For HIST data, spread counts
 *					  along the y dimension of
 *                                        the data and variance arrays
 *                                      - Incorporate the bad pixel mask into
 *                                        optimal extraction algorithm
 *                                      - Return bad-pixel spectrum
 *              09/30/03   1.10  wvd    Negative values of YCENT are 
 *					meaningless.  Use S/N of data to 
 *					determine significance of YCENT.
 *					Exclude airglow lines from S/N.
 *              10/08/03   1.11  wvd    Change counts_out array to type long.
 *					Add 0.5 to it before rounding.
 *					Limit Y coordinates of probability
 *					array before comparing with bkgd.
 *              10/31/03   1.12  wvd    Read centroid quality flags from IDF
 *					header.  Do not attempt to calculate
 *					spectral centroid.
 *              11/03/03   1.13  wvd    Modify calculation of variance estimate.
 *					Modify scheme for smoothing HIST data.
 *					Iterate only to 20% accuracy.
 *              12/05/03   1.14  wvd    Check value of SPECBINY before
 *					smoothing HIST data in the Y dimension.
 *              03/16/04   1.15  wvd    In HIST mode, assume counts are
 *					distributed across one X pixel.
 *					Delete wave_out -- not used.
 *					After 50 iterations, abandon optimal
 *					extraction and use boxcar extraction
 *					instead.
 *					Don't test size of barray,
 *					as it is now the same as parray.
 *					Change pycent from int to float.
 *              04/06/04   1.16  bjg    Include string.h
 *                                      Remove unused variables
 *              04/21/04   1.17  wvd	Set dispapix = fabs(dispapix).
 *					Iterate to 5% accuracy.
 *					In boxcar extraction, use parray to
 *					set extraction limits.
 *					For TTAG data, compute counts_out
 *					directly from photon list.
 *              04/26/04   1.18  wvd	Perform all calculations with
 *					weights array, rather than ergcm2.
 *					Return counts spectrum.
 *					Convert from variance to sigma here,
 *					rather than in cf_write_spectra.
 *					Iterate to 10% accuracy.
 *              05/17/04   1.19  wvd	If weights_out[j] = 0, set variance
 *					= bkgd.  Iterate optimal extraction
 *					to accuracy of 0.01 counts. In boxcar
 *					extraction, reduce parray limit to
 *					1E-4 to better match fluxes from 
 *					optimal extraction for bright stars.
 *              08/13/04   1.20  wvd	Use QUALITY array to set X limits 
 *					of extraction region.  Insure that
 *					sigma_out array is non-zero.
 *              11/22/05   1.21  wvd	Scale variance estimate by bad-pixel
 *					mask to better approximate observed
 *					counts.
 *              06/12/06   1.22  wvd	Call cf_verbose rather than
 *					cf_if_warning.
 *
 ****************************************************************************/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "calfuse.h"

static char 	CF_PRGM_ID[] = "cf_optimal_extraction";
static char 	CF_VER_NUM[] = "1.22";

int
cf_optimal_extraction(fitsfile *header, int optimal, int aperture, 
	float *weight, float *y, unsigned char *channel,
	float *lambda, long ngood, long *good_index,
	float *barray, float *bpmask, 
	int pny, float pycent, float *parray, long nout, float *wave_out,
	float **flux_out, float **sigma_out, long **counts_out,
	float **weights_out, float **bkgd_out, short **bpix_out)

{
    char	instmode[FLEN_VALUE], yquality_str[FLEN_KEYWORD], 
		ycent_str[FLEN_KEYWORD], yquality[FLEN_VALUE];
    int   	status=0;
    int		nloop;			/* Iterations of optimal
						extraction loop */
    int		ycent_uncertain = FALSE;/* TRUE if centroid quality
						less than HIGH */

    long	i, 			/* index through photon list */
		ii,			/* ii = good_index[i] */
		j, 			/* X index through 2-D arrays */
		jmin, jmax,		/* limits of output spectrum */
		k; 			/* Y index through 2-D arrays */
    
		 
    float	exptime,		/* exposure time */
		tdead,			/* mean dead-time correction */
		w0,			/* minimum value of output wave array */
		wpc,			/* wavelength increment per output bin */
		*flux_old,		/* flux values from previous iteration */
		*flux_box,		/* flux values from boxcar extraction */
		*var_out,		/* 1-D variance array */
		*var_box,		/* variance values from boxcar extraction */

		ycent;			/* Y centroid of target spectrum */

    double      numerator, denominator, /* Used in optimal-extraction calculations */
		var_num, weights_num,
		bkgd_num,
		*data,			/* 2-D data array */
		*variance,		/* 2-D variance array */
		*w;			/* weights array for variance calculation */

    int		l, lmin, lmax,		/* These variables are	*/
    		specbiny, specbiny2;	/* are used to smooth	*/
    float 	psum, pval;		/* HIST data in the	*/
    double	vval;			/* Y dimension.		*/

    long	jj, jjj, m, m_min, n;	/* These variables are used to smooth */
    double	edge[2], frac[2], x; 	/* HIST data in the X dimension.      */
    double	dispapix;		/* Mean dispersion per pixel */
    char	wave_file[FLEN_VALUE];  /* Name of WAVE_CAL file */
    int		hdunum;			/* HDU number of WAVE_CAL file */
    fitsfile	*wavefits;		/* Pointer to WAVE_CAL file */

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    cf_verbose(3, "Entering cf_optimal_extraction");

    /* 
     *  Check type of observation (HIST or TTAG) and binning factor.
     */
    FITS_movabs_hdu(header, 1, NULL, &status) ;
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(header, TINT, "SPECBINY", &specbiny, NULL, &status);
    specbiny2 = specbiny/2 ;
    if (specbiny2 < 1) specbiny2 = 1 ;
    cf_verbose(3, "Aperture = %d, instrument mode = %s, binning in Y = %d",
	aperture, instmode, specbiny) ;

    /*
     *  Read EXPTIME and TOT_DEAD from file header.
     */
    FITS_read_key(header, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    FITS_read_key(header, TFLOAT, "TOT_DEAD", &tdead, NULL, &status);
    if (tdead < 1.0) tdead = 1.0;
    cf_verbose (3, "EXPTIME = %.0f, TOT_DEAD=%g", exptime, tdead);

    /*
     *  Read wavelength parameters W0 and WPC from file header.
     */
    FITS_read_key(header, TFLOAT, "WPC", &wpc, NULL, &status);
    FITS_read_key(header, TFLOAT, "W0", &w0, NULL, &status);
    cf_verbose (3, "Wavelength info: w0=%g, dw=%g", w0, wpc) ;

    /*
     *  Read spectral centroid YCENT and YQUAL from file header.
     */
    sprintf(ycent_str, "YCENT%1d", aperture);
    FITS_read_key(header, TFLOAT, ycent_str, &ycent, NULL, &status);
    sprintf(yquality_str, "YQUAL%1d", aperture);
    FITS_read_key(header, TSTRING, yquality_str, yquality, NULL, &status);
    if (*yquality != 'H') ycent_uncertain = TRUE;
    cf_verbose (3, "%s = %g, quality = %s", ycent_str, ycent, yquality);

    /*
     * If HIST data, read DISPAPIX from header of WAVE_CAL file.
     */
    if (!strncmp(instmode,"H",1)) {
	hdunum = aperture + 1;
	FITS_read_key(header, TSTRING, "WAVE_CAL", wave_file, NULL, &status);
	FITS_open_file(&wavefits, cf_cal_file(wave_file), READONLY, &status);
	FITS_movabs_hdu(wavefits, hdunum, NULL, &status);
	FITS_read_key(wavefits, TDOUBLE, "DISPAPIX", &dispapix, NULL, &status);
	dispapix = fabs(dispapix);
	FITS_close_file(wavefits, &status);
    }

    /*
     *  Allocate space for all arrays.
     */
     data = (double *) cf_calloc(pny*nout, sizeof(double));
     variance = (double *) cf_calloc(pny*nout, sizeof(double));
     w = (double *) cf_calloc(nout, sizeof(double));

     flux_old = (float *) cf_calloc(nout, sizeof(float));
     flux_box = (float *) cf_calloc(nout, sizeof(float));
     var_box  = (float *) cf_calloc(nout, sizeof(float));

     *flux_out = (float *) cf_calloc(nout, sizeof(float));
     var_out = (float *) cf_calloc(nout, sizeof(float));
     *counts_out = (long *) cf_calloc(nout, sizeof(long));
     *weights_out = (float *) cf_calloc(nout, sizeof(float));
     *bkgd_out = (float *) cf_calloc(nout, sizeof(float));
     *bpix_out = (short *) cf_calloc(nout, sizeof(short)) ;

    /* 
     *  The quality array (bpix_out) is the product of the probability
     *  array and the bad-pixel mask.  We use it to set the limits of
     *  the extraction region.
     */
     cf_verbose(3, "Generating bad-pixel spectrum.") ;
     for (j=0; j<nout; j++) {
         float bpix=0. ;
         for (k=0; k<pny; k++)
	   bpix += (parray[k*nout + j] * bpmask[k*nout + j]);
         (*bpix_out)[j] = (short) (100. * bpix + 0.5) ;
     }
     jmin = 0L; jmax = nout;
     for (j=0; j<nout; j++) {
	if ((*bpix_out)[j] > 0) {
	    jmin = j;
	    break;
	}
     }
     for (j = nout-1; j >= 0; j--) {
	if ((*bpix_out)[j] > 0) {
	    jmax = j+1;
	    break;
	}
     }

    /*
     *  Construct data and variance arrays from photon lists.
     *  Arrays are pny pixels in Y (to match probability array) 
     *  and nout pixels in X (to match output wavelength array).  
     *  We must thus shift individual photons by PYCENT - YCENT
     *  pixels in Y.
     */

     cf_verbose(3, "Filling the data array") ;
     for (i = 0; i < ngood; i++) {
	ii = good_index[i];

	if (channel[ii] != aperture) continue;

	j = (long) ((lambda[ii] - w0)/wpc + 0.5);
	if (j < jmin || j >= jmax) continue;

	k = (long) (y[ii] - ycent + pycent + 0.5);
	if (k < 0 || k >= pny) continue;

        if (!strncmp(instmode,"T",1)) {
	   /* TTAG data */
	    (*counts_out)[j]++;
	    data[k*nout + j] += weight[ii];
	    variance[k*nout + j] += weight[ii] * weight[ii];
	}
	else { 
	   /* If HIST data are already smoothed in Y, it's easy: */
	   if (specbiny == 1) {
	       data[k*nout + j] += weight[ii];
	       variance[k*nout + j] += weight[ii] * tdead;
	   }
	   /* If HIST data are binned, smooth in both X and Y dimensions. */
	    else {
		m_min = 0;
		edge[0] = ((double) j - 0.5) * wpc + w0;
		edge[1] = ((double) j + 0.5) * wpc + w0;
		if ((x = (lambda[ii] - edge[0]) / dispapix) < 0.5) {
		   frac[0] = 0.5 - x;
		   frac[1] = 0.5 + x;
		   jj = j - 1;
		   n = 2;
		   if (jj == -1) m_min = 1;
		}
		else if ((x = (edge[1] - lambda[ii]) / dispapix) < 0.5) {
		   frac[0] = 0.5 + x;
		   frac[1] = 0.5 - x;
		   jj = j;
		   n = 2;
		   if (jj + 1 == nout) n = 1;
		}
		else {
		   frac[0] = 1.0;
		   jj = j;
		   n = 1;
		}

		psum = 0.;
		vval = weight[ii] * tdead;
		lmin = k - specbiny2;
		lmax = lmin + specbiny;
		if (lmin < 0) lmin = 0;
		if (lmax > pny) lmax = pny;
		pval = 1. / (lmax - lmin);
		for (l = lmin; l < lmax; l++)
		    psum += parray[l*nout + j];
		for (l = lmin; l < lmax; l++) {
		    if (psum > 0) pval = parray[l*nout + j] / psum;
		    for (m = m_min; m < n; m++) {
			jjj = jj + m;
			data[l*nout + jjj] += pval * frac[m] * weight[ii];
			variance[l*nout + jjj] += pval * frac[m] * vval;
		    }
		}
	    }
	}
     }

    /*
     *  Boxcar extraction.  Use parray to set extraction limits in Y.
     */
     cf_verbose(3, "Doing the boxcar extraction.") ;
     nloop = 0;
     for (j = jmin; j < jmax; j++) {
	for (k = 0; k < pny; k++) {
	   if (parray[k*nout + j] > 1E-4) {
	         (*bkgd_out)[j] += barray[k*nout + j];
	             var_out[j] += variance[k*nout + j];
	      (*weights_out)[j] += data[k*nout + j];
	   }
	}
	(*flux_out)[j] = (*weights_out)[j] - (*bkgd_out)[j];
     }

    /*
     *  Has the user requested optimal extraction?
     */
     cf_verbose(3,"Optimal keyword = %d ", optimal) ;
     if (optimal) {

       /*
	*  If Y centroid is uncertain, abandon optimal extraction.
	*/
	if (ycent_uncertain)
	   cf_verbose(1, "%s uncertain.  Cannot perform optimal extraction.",
		ycent_str);

	else {		/* Proceed with optimal extraction */

 	    cf_verbose(3, "Attempting optimal extraction");

	/*
	 *  Save boxcar versions of flux and variance arrays.
	 */
	for (j = jmin; j < jmax; j++) {
	    flux_box[j] = (*flux_out)[j];
	    var_box[j] = var_out[j];
	}

	/*
	 *  Need array w[j] to scale the variance estimate.
	 *  w[j] is mean scale factor from counts to weights at each wavelength.
	 *  Must calculate for TTAG data; can use TOT_DEAD for HIST data.
	 */
	if (!strncmp(instmode,"T",1)) for (j = jmin; j < jmax; j++) {
	    if ((*counts_out)[j] > 0) {
		for (k = 0; k < pny; k++) w[j] += data[k*nout + j];
		w[j] /= (*counts_out)[j];
	    }
	    else w[j] = tdead;
	}
	else for (j = jmin; j < jmax; j++) w[j] = tdead;

	/*
	 *  Begin big loop.
	 */
	do {
	    /*
	     *  Revise variance estimate
	     *  We can omit w[j] here, as it appears in both the numerator
	     *  and denominator of the flux estimate (below).
	     *  05/17/02 wvd:  If weights_out[j] = 0, set variance = bkgd.
	     *  11/22/05 wvd:  Scale variance estimate by bpmask.
	     */
	     for (j = jmin; j < jmax; j++) {
		if ((*weights_out)[j] > 0)
		   for (k = 0; k < pny; k++)
			variance[k*nout + j] = bpmask[k*nout + j] *
			fabs(parray[k*nout + j] * (*flux_out)[j] + barray[k*nout + j]);
		else
		   for (k = 0; k < pny; k++)
			variance[k*nout + j] = bpmask[k*nout + j] * barray[k*nout + j];
	     }

	    /*
	     *  Optimal extraction
	     */
	     nloop++;
	     for (j = jmin; j < jmax; j++) {
		numerator = denominator = 0.;
		for (k = 0; k < pny; k++) if (variance[k*nout+j] > 0) {
		    numerator += bpmask[k*nout+j] * parray[k*nout + j] * 
                      (data[k*nout + j] - barray[k*nout + j]) / variance[k*nout + j];
		    denominator += bpmask[k*nout+j] * parray[k*nout + j] * 
                       parray[k*nout + j] / variance[k*nout + j];
		}
		flux_old[j] = (*flux_out)[j];
		if (denominator > 0) (*flux_out)[j] = numerator / denominator;
                else (*flux_out)[j] = flux_old[j] = 0. ;
	     }

	    /*
	     *  Repeat loop if any element of flux_out has changed by > 0.01 counts.
	     */
	     for (j = jmin; j < jmax; j++) {
		if (fabs((*flux_out)[j]-flux_old[j]) > 0.01) {
		    cf_verbose(3, "%d\t%d\t%g", nloop, j, (*flux_out)[j]);
		    break;
		}
	     }

	} while (j < jmax && nloop < 50);	/* End of do loop */

	/*
	 *  If nloop < 50, calculate var_out, weights_out, and bkgd_out.
	 *  Note that we must scale var_out by w[j] here.
	 */
	if (nloop < 50) {

	    memset (*bkgd_out, 0, nout*sizeof(float));
	    memset (*weights_out, 0, nout*sizeof(float));
	    memset (var_out, 0, nout*sizeof(float));

	    for (j = jmin; j < jmax; j++) {
		bkgd_num = weights_num = var_num = denominator = 0.;
		for (k = 0; k < pny; k++) if (variance[k*nout+j] > 0) {
		   bkgd_num += bpmask[k*nout+j] * parray[k*nout + j] *
                      barray[k*nout + j] / variance[k*nout + j];
		   weights_num += bpmask[k*nout+j] * parray[k*nout + j] *
                      data[k*nout + j] / variance[k*nout + j];
		   var_num += bpmask[k*nout+j] * parray[k*nout + j];
		   denominator += bpmask[k*nout+j] * parray[k*nout + j] * 
                        parray[k*nout + j] / variance[k*nout + j];
		}
		if (denominator > 0) {
		    (*bkgd_out)[j] = bkgd_num / denominator;
		    (*weights_out)[j] = weights_num / denominator;
		    var_out[j] = var_num / denominator * w[j];
		}
		else {
		    (*bkgd_out)[j] = (*weights_out)[j] = var_out[j] = 0.;
		}
	    }
	}
	/*
	 *  If nloop = 50, use boxcar extraction instead.
	 */
	else {
	    cf_verbose(1, "Optimal extraction failed to converge for aperture %d."
		"  Using boxcar extraction.", aperture);
	    nloop = 0;
	    for (j = jmin; j < jmax; j++) {
		(*flux_out)[j] = flux_box[j];
		var_out[j] = var_box[j];
            }
	}

	}			/* End of if/else test of YCENT */
     }				/* End of optimal extraction block */

    /*
     *  Write number of iterations to file header.
     */
     FITS_update_key(header, TINT, "OPT_EXTR", &nloop, NULL, &status);
     cf_verbose(2, "Number of iterations in optimal extraction = %d", nloop) ;

    /* Convert error array from variance to sigma. */
    for (j = jmin; j < jmax; j ++) {
	if (var_out[j] > 0.) var_out[j] = sqrt(var_out[j]);
	else if ((*bkgd_out)[j] > 0.) var_out[j] = sqrt((*bkgd_out)[j]);
	else var_out[j] = 0.;
    }
    *sigma_out = var_out;

    /*
     *  For HIST data, counts array = weights "uncorrected" for dead time.
     */
     if (!strncmp(instmode,"H",1)) for (j = jmin; j < jmax; j++)
	(*counts_out)[j] = (long) ((*weights_out)[j] / tdead + 0.5);

    /*
     *  Clean up
     */
    free(data);
    free(flux_old);
    free(flux_box);
    free(var_box);
    free(variance);
    free(w);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
