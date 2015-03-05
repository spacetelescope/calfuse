/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Usage:    cf_xcorr exposure_list exp_shift_list
 *
 * Description: Determine pixel shifts between exposures. The
 *              reference spectrum is the one having the longest exposure
 *              time. The input list (ASCII) is the set of exposures to process,
 *              the output list (ASCII) has 3 columns: exposure name, shift and
 *              sigma_shift. Sigma_shift = -1 corresponds to non-detected shift.
 *		Sigma_shift = 1 indicates a valid shift value.
 *              The reference exposure always has shift and sigma_shift of 0.
 *
 * History:     04/07/05        tc      v1.0    First release
 *              04/15/05        tc      v1.1    Use cf_read_col
 *              04/27/05       wvd      v1.2    Use variance arrays properly.
 *						Use QUALITY array to identify
 *						   bad data.
 *						Change spec array to type INT.
 *						Use chi2, not reduced chi2,
 *						   to set error bars.
 *              05/20/05       wvd      v1.3    For continuum spectra, use
 *						single region between 1045 and
 *						1070 A.  For emission-line
 *						targets, use O VI lines.
 *              06/03/05       wvd      v1.4    Delete unused variables.
 *              07/21/05       wvd      v1.5    Adopt the sign convention for 
 *						spectral shifts used by
 *						FUSE_REGISTER.
 *              08/01/05       wvd      v1.6    Store EXPNIGHT of reference
 *						spectrum in ref_expnight.
 *              09/13/05       wvd      v1.7    If called with no arguments,
 *						return calling info, not error.
 *              04/27/06       wvd      v1.8    Scale each spectrum to match
 *						mean of reference spectrum.
 *              05/16/06       wvd      v1.9    Use Lyman beta to align 
 *						background exposures.
 *						Use O VI and C II to align WD's.
 *              03/28/08       wvd      v1.10   Write value of NORM to output.
 *						For HIST data only: if XCORRR
 *						fails, but NORM > 0.5, then
 *						set shift and sigma to 0.
 *              08/15/08       wvd      v1.11   If EXP_STAT = 2, treat as
 *						background observation.
 *
 ******************************************************************************/
 
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

/* Calfuse variables */
static char CF_PRGM_ID[] = "cf_xcorr";
static char CF_VER_NUM[] = "1.11";

/* Parameters */
#define MAXSHIFT     20
#define N_SHIFT	     41

int main(int argc, char *argv[])
{
	/* Variables */
	char *spec_lis, *out_name, src_type[FLEN_VALUE];
	char spec_name[FLEN_FILENAME], ref_spec_name[FLEN_FILENAME];
	char program[FLEN_VALUE]; 
	char instmode[FLEN_CARD] ;
	double *ref_flux, *ref_error, *flux, *error;
	double scale, sum, variance, norm, ref_tot, total;
	float chi_square[N_SHIFT], chi_square_min, chi_square_max;
	float exptime, exptime_max;
	float *wave, w0, wpc, wmin, wdelta;
	int ly_beta=FALSE, exp_stat, objclass, shift, sigma;
	int ind_chi2_min, N, N_tmp, nchi, n_spec, ind_ref_spec=0;
	long expnight, ref_expnight, i, k, n_pix, ref_start, start;
    
	FILE *pt_file, *pt_file_out;
	fitsfile *in_fits;
	int status = 0, hdutype = 0;
	
	/***********************************
	** Enter a timestamp into the log **
	***********************************/
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
	cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);


	/***********************************************************	
	** Check for proper command-line usage and read arguments **
	***********************************************************/
	if (argc != 3) {
		printf("Usage: cf_xcorr exposure_list exp_shift_list\n");
		return(0);
	}

	spec_lis = argv[1];
	out_name = argv[2];
	
	/********************************************	
	** Open input list and create output file  **
	********************************************/
	if ((pt_file = fopen(spec_lis, "r")) == NULL)
	{
		cf_if_error("Unable to open file %s", spec_lis);
	}
	pt_file_out = fopen(out_name, "w");

	/********************************************************************
	** Set the reference spectrum to the one which has maximum EXPTIME **
	********************************************************************/
	n_spec = 0;
	exptime_max = -1;
    
	while (fscanf(pt_file, "%80s", spec_name) != EOF)
	{
		FITS_open_file(&in_fits, spec_name, READONLY, &status);
		FITS_read_key(in_fits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
		FITS_read_key(in_fits, TLONG, "EXPNIGHT", &expnight, NULL, &status);
		FITS_close_file(in_fits, &status);
        	if (exptime > exptime_max)
		{
			exptime_max = exptime;
			ind_ref_spec = n_spec;
			strcpy(ref_spec_name, spec_name);
		}
		n_spec++;
	}
	
	if (n_spec == 0) cf_if_error("At least one spectrum is required");
	
	if (n_spec == 1)
	{
		fprintf(pt_file_out, "%s %2d %2d %5.0f %5ld %c%c%c 1.0\n", ref_spec_name, 0, 0,
			exptime_max, expnight, ref_spec_name[8], ref_spec_name[9], ref_spec_name[10]);
		fclose(pt_file_out);
		return(0);
	}
	printf("n_spec: %d, max exptime: %g, ind_ref: %d, ref_name: %s\n",
		n_spec, exptime_max, ind_ref_spec, ref_spec_name);
		
        /* Read source type and set XCORR limits accordingly. */
	FITS_open_file(&in_fits, ref_spec_name, READONLY, &status);
        FITS_read_key(in_fits, TINT, "EXP_STAT", &exp_stat, NULL, &status);
        FITS_read_key(in_fits, TSTRING, "PRGRM_ID", &program, NULL, &status);
        FITS_read_key(in_fits, TINT, "OBJCLASS", &objclass, NULL, &status);
        FITS_read_key(in_fits, TSTRING, "SRC_TYPE", &src_type, NULL, &status);
	FITS_read_key(in_fits, TSTRING, "INSTMODE", instmode, NULL, &status) ;
	FITS_read_key(in_fits, TFLOAT, "W0", &w0, NULL, &status);
	FITS_read_key(in_fits, TFLOAT, "WPC", &wpc, NULL, &status);
        if (src_type[1] == 'E') {       /* Emission-line source */
			/* Use Lyman beta to align background observations. */
	    if (exp_stat == (int) TEMPORAL_LIMB || 
		objclass == 1 || objclass == 7 || objclass == 90) {
                printf("Assuming background observation.\n");
            	wmin = 1024.;
            	wdelta = 6.;
		ly_beta = TRUE;
	    } else {	/* Use O VI emission for everything else. */
                printf("SRC_TYPE = %s.  Emission-line target.\n", src_type);
            	wmin = 1030.;
            	wdelta = 9.;
	    }
        }
	else if (objclass == 17 || objclass == 29 || objclass == 37) {
			/* Use O VI and C II to align white dwarf spectra. */
            printf("OBJCLASS = %d.  Assuming nearby white dwarf.\n", objclass);
            wmin = 1030.;
            wdelta = 9.;
        }
        else {  /* Use 1045-1070 A region for all other continuum sources. */
            printf("SRC_TYPE = %s.  Assuming continuum target.\n", src_type);
            wmin = 1045.;
            wdelta = 25.;
        }
        start = cf_nlong((wmin - w0) / wpc);
        ref_start = start + MAXSHIFT;
        n_pix = cf_nlong(wdelta / wpc) - (N_SHIFT - 1);

	/* Read wave, flux, and error arrays from reference exposure. */	
	FITS_read_key(in_fits, TLONG, "EXPNIGHT" , &ref_expnight, NULL, &status);
	FITS_movabs_hdu(in_fits, 2, &hdutype, &status);
	if (hdutype != BINARY_TBL) cf_if_error("FITS files must contain BINARY TABLE in HDU #2");
	N = cf_read_col(in_fits, TFLOAT, "WAVE", (void **) &wave);
	N = cf_read_col(in_fits, TDOUBLE, "FLUX", (void **) &ref_flux);
	N = cf_read_col(in_fits, TDOUBLE, "ERROR", (void **) &ref_error);
	FITS_close_file(in_fits, &status);

	/* Compute total flux of reference spectrum in region of interest. */
	ref_tot = 0.;
	for (i = 0; i < n_pix; i++) ref_tot += ref_flux[start+i];
	
	/*****************************************************************
	** Compute the shift of each spectrum relative to the reference **
	******************************************************************/
	n_spec = 0;
	rewind(pt_file);
	
	while (fscanf(pt_file, "%80s", spec_name) != EOF)
	{
		/* Do nothing if the current spectrum is the reference. */
		if (n_spec++ == ind_ref_spec) {
			fprintf(pt_file_out, "%s %3d %2d %5.0f %5ld %c%c%c 1.0\n",
				ref_spec_name, 0, 0, exptime_max, ref_expnight,
				ref_spec_name[8], ref_spec_name[9], ref_spec_name[10]);
			continue;
		}

		/* Read the next spectrum */
		FITS_open_file(&in_fits, spec_name, READONLY, &status);
		FITS_read_key(in_fits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
		FITS_read_key(in_fits, TLONG, "EXPNIGHT", &expnight, NULL, &status);
		FITS_movabs_hdu(in_fits, 2, &hdutype, &status);
		if (hdutype != BINARY_TBL) cf_if_error("FITS files must contain BINARY TABLE in HDU #2");
		N_tmp = cf_read_col(in_fits, TDOUBLE, "FLUX", (void **) &flux);
		N_tmp = cf_read_col(in_fits, TDOUBLE, "ERROR", (void **) &error);
		if (N_tmp != N) cf_if_error("Tables must have the same number of elements");
		FITS_close_file(in_fits, &status);

		/* If EXPTIME < 1, move on to the next data file. */
		if (exptime < 1) {
		    fprintf(pt_file_out, "%s %3d %2d %5.0f %5ld %c%c%c 0.0\n",
			spec_name, 0, -1, exptime, expnight,
			spec_name[8], spec_name[9], spec_name[10]);
		    continue;
		}

		/* Compute total flux of data spectrum in region of interest. */
		norm = 1.;
		if (!ly_beta) {		/* Don't rescale if aligning on airglow. */
		    total = 0.;
		    for (i = 0; i < n_pix; i++) total += flux[start+i];
		    norm = ref_tot / total;
		}
	
        	/* Compute chi-squared for shifts between +/- MAXSHIFT pixels. */
        	for (k = 0; k < N_SHIFT; k++) {
		    nchi = 0;
            	    sum = 0;
            	    for (i = 0; i < n_pix; i++) {
                	variance = ref_error[ref_start + i] * ref_error[ref_start + i] +
                           norm * norm * error[start + i + k] * error[start + i + k];
                        if (variance > 0) {
                            sum += (ref_flux[ref_start + i] - norm * flux[start + i + k]) *
                               (ref_flux[ref_start + i] - norm * flux[start + i + k]) / variance;
                            nchi++;
                        }
                    }
                    scale = (double) nchi / n_pix;
                    chi_square[k] = sum / scale;
                }

		/* Compute min and max values of chi-squared. */
        	ind_chi2_min = 0;
        	chi_square_min = 1E5;
        	chi_square_max = -1E5;
        	for (k = 0; k < N_SHIFT; k++) {
                	if (chi_square[k] < chi_square_min)
                	{
                        	chi_square_min = chi_square[k];
                        	ind_chi2_min = k;
                	}
                	else if (chi_square[k] > chi_square_max)
                        	chi_square_max = chi_square[k];
        	}
                shift = MAXSHIFT - ind_chi2_min;

        	/* If chi-squared changes by less than 20% over range of shifts, bail out.
        	 * Otherwise, return shift corresponding to lowest value of chi-square. */
        	if (chi_square_max / chi_square_min < 1.2) 
			sigma = -1;
        	else 
			sigma = 1;


		/* If we're in HIST mode and sigma = -1 and NORM < 2,
		 * then set shift = sigma = 0. */
		if (!(strncasecmp(instmode, "HIST", 4)) && sigma == -1 && norm < 2.0) 
			shift = sigma = 0;

		/* printf ("shift = %d\tchi_square_max / chi_square_min = %g\n", 
			ind_chi2_min - MAXSHIFT, chi_square_max / chi_square_min); */
		
		fprintf(pt_file_out, "%s %3d %2d %5.0f %5ld %c%c%c %f\n",
			spec_name, shift, sigma, exptime, expnight,
			spec_name[8], spec_name[9], spec_name[10], norm);

		free(flux);
		free(error);
	}
	
	free(wave);
	free(ref_flux);
	free(ref_error);
	
	fclose(pt_file);
	fclose(pt_file_out);

	return(0);
}
