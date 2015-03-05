/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Usage:    cf_nvo all_segments.fit nvo_file.fit
 *
 * Description: Write a National Virtual Observatory (nvo) file. The program
 *              uses the file created by cf_pack, which contains one spectrum
 *              per extension, each corresponding to a detector channel.
 *
 * History:     04/11/05        tc      v1.0    First release
 *		05/16/2005	wvd	v1.1	Check mean flux of each segment.
 *						If preferred channel is missing,
 *						replace it with another.
 *		05/19/2005	wvd	v1.2	Shift each channel to match
 *						LiF1A between 1045 and 1070 A.
 *		05/20/2005	wvd	v1.3	For emission-line sources, 
 *						cross-correlate on O VI lines.
 *		05/23/2005	wvd	v1.4	Don't assume that sides A and B
 *						have the same shift.
 *		06/01/2005	wvd	v1.5	Do assume that sides A and B
 *						have the same shift.
 *		06/03/2005	wvd	v1.6	Delete unused variables.
 *		06/08/2005	wvd	v1.7	Define MAXFLOAT if needed.
 *		06/16/2005	wvd	v1.8	Include values.h
 *		07/06/2005	wvd	v1.9	If FESCENT = FES B, use LiF 2B
 *						as wavelength standard.
 *		07/12/2005	wvd	v1.10	Give up on use of MAXFLOAT.
 *		08/11/2005	wvd	v1.11	Set keyword NEXTEND = 1.
 *		03/22/2006	wvd	v1.12	Allow use of SiC data for 
 *						1000-1100 A region.
 *						Don't cross-correlate segments
 *						with OBSTIME = 0.
 *		05/16/2006	wvd	v1.13	Use Lyman beta to align
 *                                              background exposures.
 *                                              Use O VI and C II to align WD's.
 *		05/19/2006	wvd	v1.14	In copy_spec(), test index to 
 *						prevent extending past the edge.
 *						For BKGD targets, omit flux
 *						comparison when deciding which
 *						regions to use.
 *						Delete index to the file 
 *						extensions from the primary HDU.
 *						Compute mean flux over same
 *						wavelength region for each band.
 *		05/24/2006	wvd	v1.16	For PC targets, don't include
 *						O VI in calculation of mean flux
 *		12/19/2006	wvd	v1.17	Reject segments with 
 *						OBSTIME < 10 seconds when full
 *						exposure is longer than 100 s.
 *		02/19/2008	bot	v1.18	Fixed indexing problem in
 *						copy_spec due to shift and n_copy
 *
 ******************************************************************************/
 
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

#define MAXSHIFT	30
#define N_SHIFT		61
#define NVO_MIN 	900
#define NVO_MAX 	1190

static char CF_PRGM_ID[] = "cf_nvo";
static char CF_VER_NUM[] = "1.18";

static int
compute_shift(fitsfile *pt_fits, int hdu_ref, int hdu, float wmin, float wdelta,
	float wpc, int *shift, int *quality)
{
	char	comb[FLEN_VALUE], ref_comb[FLEN_VALUE];
        int     hdutype, status = 0;
	int	ind_chi2_min, nchi;
        long    i, k, n_rows, ref_start, start, n_pix;
        float   *error, *flux, *wave, *ref_error, *ref_flux, *ref_wave;
	float	chi_square[N_SHIFT], chi_square_min, chi_square_max;
	float	ref_obstime, obstime;
	double	scale, sum, var1, var2;

        /* Don't need to shift reference spectrum. */
        if (hdu_ref == hdu) {
                *shift = 0;
                *quality = TRUE;
                return(0);
        }

	/* Read data for reference and data channels. */
        FITS_movabs_hdu(pt_fits, hdu_ref, &hdutype, &status);
        FITS_read_key(pt_fits, TSTRING, "COMBMETH", &ref_comb, NULL, &status);
        FITS_read_key(pt_fits, TFLOAT, "OBSTIME", &ref_obstime, NULL, &status);
        n_rows = cf_read_col(pt_fits, TFLOAT, "WAVE", (void **) &ref_wave);
        n_rows = cf_read_col(pt_fits, TFLOAT, "FLUX", (void **) &ref_flux);
        n_rows = cf_read_col(pt_fits, TFLOAT, "ERROR", (void **) &ref_error);

        FITS_movabs_hdu(pt_fits, hdu, &hdutype, &status);
        FITS_read_key(pt_fits, TSTRING, "COMBMETH", &comb, NULL, &status);
        FITS_read_key(pt_fits, TFLOAT, "OBSTIME", &obstime, NULL, &status);
        n_rows = cf_read_col(pt_fits, TFLOAT, "WAVE", (void **) &wave);
        n_rows = cf_read_col(pt_fits, TFLOAT, "FLUX", (void **) &flux);
        n_rows = cf_read_col(pt_fits, TFLOAT, "ERROR", (void **) &error);

	/* If either channel was not constructed using cross-correlation, 
	 * return a shift of 0. */
	if (strncmp(ref_comb,"X",1) || strncmp(comb,"X",1)) {
		*shift = 0;
                *quality = FALSE;
		printf("hdu = %d, hdu_ref = %d, COMBMETH != XCORR\n",
			hdu, hdu_ref);
		return(0);
	}

	/* If either channel has OBSTIME = 0, return a shift of 0. */
	if (ref_obstime < 1 || obstime < 1) {
		*shift = 0;
                *quality = FALSE;
		if (ref_obstime < 1)
			printf("hdu_ref = %d, OBSTIME = %f\n", hdu_ref, ref_obstime);
		if (obstime < 1) printf("hdu = %d, OBSTIME = %f\n", hdu, obstime);
		return(0);
	}

	/* Compute chi-squared for shifts between +/- MAXSHIFT pixels. */
        ref_start = cf_nlong((wmin - ref_wave[0]) / wpc) + MAXSHIFT;
        start = cf_nlong((wmin - wave[0]) / wpc);
        n_pix = cf_nlong(wdelta / wpc) - (N_SHIFT - 1);

	for (k = 0; k < N_SHIFT; k++) {
	    nchi = 0;
	    sum = 0;
	    for (i = 0; i < n_pix; i++) {
		var1 = ref_error[ref_start + i] * ref_error[ref_start + i];
		var2 = error[start + i + k] * error[start + i + k];
		if ((var1 + var2) > 0.) {
		    sum += (ref_flux[ref_start + i] - flux[start + i + k]) *
		           (ref_flux[ref_start + i] - flux[start + i + k]) /
			   (var1 + var2);
		    nchi++;
		}
	    }
	    scale = (double) nchi / n_pix;
	    chi_square[k] = sum / scale;
	}

	/* If chi-squared changes by less than 20% over range of shifts, abort.
	 * Otherwise, return shift corresponding to lowest value of chi-square.
	 */
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

        if (chi_square_max / chi_square_min < 1.2) {
                *shift = 0;
                *quality = FALSE;
		printf("hdu = %d, hdu_ref = %d.  Chi-square dist is flat.\n",
		hdu, hdu_ref);
	}
        else {
                *shift = ind_chi2_min - MAXSHIFT;
                *quality = TRUE;
	}

        free(ref_error);
        free(ref_flux);
        free(ref_wave);
        free(error);
        free(flux);
        free(wave);

	return(0);
}


static int
copy_spec(fitsfile *pt_fits, int hdu, float w_min, float w_max, float wpc,
	int shift, float *flux_nvo, float *error_nvo, long n_nvo)
{
	int	hdutype, status = 0;
	long	i, i0, n_rows, nvo_start, start, n_copy, n_copy0;
	float	*error, *flux, w0;

	FITS_movabs_hdu(pt_fits, hdu, &hdutype, &status);
	FITS_read_key(pt_fits, TFLOAT, "W0", &w0, NULL, &status);
	n_rows = cf_read_col(pt_fits, TFLOAT, "FLUX", (void **) &flux);
	n_rows = cf_read_col(pt_fits, TFLOAT, "ERROR", (void **) &error);

	nvo_start = cf_nlong(((double) w_min - NVO_MIN) / wpc);
	start = cf_nlong(((double) w_min - w0) / wpc);
	n_copy = cf_nlong(((double) w_max - w_min) / wpc) + 1;
	if (nvo_start + n_copy > n_nvo) {
	    printf ("HDU = %d, w_min = %.1f, w_max = %.1f, wpc = %5.3f\n",
		hdu, w_min, w_max, wpc);
	    printf ("n_nvo = %ld, nvo_start = %ld, start = %ld, n_copy = %ld\n",
		n_nvo, nvo_start, start, n_copy);
	    cf_if_error("Can't create NVO array.");
	}
	 
	n_copy0 = n_copy;   
	if (start + shift + n_copy > n_rows) {
	    n_copy = n_rows - (start + shift);
	    if (hdu == 3) for (i = n_copy; i < n_copy0; i++) {
		flux_nvo[i+nvo_start]  = 0.0;
		error_nvo[i+nvo_start] = 0.0;
	    }
	}

	i0 = 0;
	if (start + shift < 0) {
	    i0 = -(start + shift);
	    for (i = 0; i < i0; i++) {
		flux_nvo[i+nvo_start]  = 0.0;
		error_nvo[i+nvo_start] = 0.0;
	    }
	}

	for (i = i0; i < n_copy; i++) {
	    flux_nvo[i+nvo_start]  = flux[i+start+shift];
	    error_nvo[i+nvo_start] = error[i+start+shift];
	}
 	   
	free(error);
	free(flux);
	return(0);
}


int main(int argc, char *argv[])
{
	/* Variables */
        char program[FLEN_VALUE];
	double mean_flux[8], sum_flux;
	float *wave, *flux, w0, wpc, hdu_wpc;
	float *wave_nvo, *flux_nvo, *error_nvo;
	float max_obstime, wmin, wdelta;
	long n_rows, n_nvo;
	int ref, i, nextend=1, n_hdus, shift[8], quality[8];
	int bkgd_obs=FALSE, objclass;
	int nkeys, morekeys=0;
	long k, n_elem;
	
	char *ttype[] = { "WAVE", "FLUX", "ERROR" };
	char *tform[] = { "1E", "1E", "1E" };
	char *tunit[] = { "ANGSTROMS", "ERG/CM2/S/A", "ERG/CM2/S/A" };
	
	/* HDU numbers must be consistent with output of cf_pack. */
	char *extname[] = { "1ASIC", "2BSIC", "1ALIF", "2BLIF", "1BSIC", "2ASIC", "1BLIF", "2ALIF" };
	int hdu[]       = { 6,       8,       2,       4,       7,       9,       3      , 5       };
	float w_min[]   = {  997,    1010,    990,     990,     904,     917.5,   1094   , 1090    };
	float w_max[]   = { 1090,    1104,    1080,    1074,    992,     998,     1190   , 1180    };

	/* Regions to use when computing mean flux. */
	float *region;
	float region_array[8][6] = {
		{0, 0, 1030, 1039, 1045, 1070},			/* 1ASIC */
	        {0, 0, 1030, 1039, 1045, 1070},			/* 2BSIC */
	        {0, 0, 1030, 1039, 1045, 1070},			/* 1ALIF */
	        {0, 0, 1030, 1039, 1045, 1070},			/* 2BLIF */
	        {912, 935, 955, 970, 980, 985},			/* 1BSIC */
	        {912, 935, 955, 970, 980, 985},			/* 2ASIC */
	        {1095, 1130, 1140, 1165, 1170, 1180},		/* 1BLIF */
	        {1095, 1130, 1140, 1165, 1170, 1180}		/* 2ALIF */
	};

	char in_name[80], out_name[80];
	char extname_key[80], src_type[FLEN_VALUE], sp_type[FLEN_VALUE], fescent[FLEN_VALUE];
	fitsfile *pt_fits, *pt_fits_out;
	int status = 0, hdutype;

	    
	/***********************************
	** Enter a timestamp into the log **
	***********************************/
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
	cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);


	/***********************************************************	
	** Check for proper command-line usage and read arguments **
	***********************************************************/
	if (argc != 3)
		cf_if_error("Usage: cf_nvo all_segments.fit nvo_file.fit");
	
	strcpy(in_name, argv[1]);
	strcpy(out_name, "!");
	strcat(out_name, argv[2]); /* Force overwrite */
	
	/*******************************************	
	** Open input file, read header keywords. **
	*******************************************/
	FITS_open_file(&pt_fits, in_name, READONLY, &status);
        FITS_read_key(pt_fits, TSTRING, "PRGRM_ID", &program,  NULL, &status);
        FITS_read_key(pt_fits, TINT,    "OBJCLASS", &objclass, NULL, &status);
        FITS_read_key(pt_fits, TSTRING, "SP_TYPE",  &sp_type,  NULL, &status);
        FITS_read_key(pt_fits, TSTRING, "SRC_TYPE", &src_type, NULL, &status);
	FITS_read_key(pt_fits, TFLOAT,  "OBSTIME",  &max_obstime, NULL, &status);
	FITS_read_key(pt_fits, TFLOAT,  "WPC",      &wpc,      NULL, &status);
        FITS_read_key(pt_fits, TSTRING, "FESCENT",  &fescent,  NULL, &status);
	ref = 2;					/* LiF 1A */
	if (!strncmp(fescent, "FES B", 5)) ref = 3;	/* LiF 2B */

	/****************************************************************	
	** If SRC_TYPE = PC, exclude O VI in calculation of mean flux. **
	****************************************************************/
	if (!strncmp(src_type, "PC", 2))
	    for (i = 0; i < 4; i++)
		region_array[i][2] = region_array[i][3] = 0;

	/***************************************	
	** Compute mean flux in each channel. **
	***************************************/
	FITS_get_num_hdus(pt_fits, &n_hdus, &status);
	if (n_hdus != 9)
		cf_if_error("Input file must have 9 extensions (%d found)", n_hdus);
	
	for (i = 0; i < 8; i++) {
		float obstime;
		FITS_movabs_hdu(pt_fits, hdu[i], &hdutype, &status);
		FITS_read_key(pt_fits, TSTRING, "EXTNAME", extname_key, NULL, &status);
		if (strcmp(extname_key, extname[i]))
			cf_if_error("Extension %d's name does not match (%s)",
			hdu[i], extname[i]);
		FITS_read_key(pt_fits, TFLOAT, "WPC", &hdu_wpc, NULL, &status);
		if (hdu_wpc != wpc) 
			cf_if_error("Extension %d's WPC does not match primary HDU",
			hdu[i]);
        	FITS_read_key(pt_fits, TFLOAT, "OBSTIME", &obstime, NULL, &status);
		if (obstime < 1. || (obstime < 10 && max_obstime > 100))
			mean_flux[i] = -1;
		else {
		   n_rows = cf_read_col(pt_fits, TFLOAT, "WAVE", (void **) &wave);
		   n_rows = cf_read_col(pt_fits, TFLOAT, "FLUX", (void **) &flux);

		   sum_flux = n_elem = 0;
		   region = region_array[i];
		   for (k = 0; k < n_rows; k++)
		   {
			if ((wave[k] > region[0] && wave[k] < region[1]) ||
			    (wave[k] > region[2] && wave[k] < region[3]) ||
			    (wave[k] > region[4] && wave[k] < region[5])) {
				sum_flux += flux[k];
				n_elem ++;
			}
		   }
		
		   mean_flux[i] = sum_flux / n_elem;
		   free(wave);
		   free(flux);
		}

		printf("%s\tmean_flux = %g\n", extname_key, mean_flux[i]);
	}
	
	/*******************************************************************
	 * If FESCENT = FES A, compute shifts relative to LiF 1A. 
	 * If FESCENT = FES B, use LiF 2B.
	 * For background targets, compare Lyman beta lines.
	 * For other emission-line targets, compare O VI lines.
	 * For white dwarfs, use O VI and C II.
	 * For other continuum sources, use the region between 1045 and 1070 A
	 * for the four 1000-1100 A channels.
	 ******************************************************************/

	if (src_type[1] == 'E') {	/* Emission-line source */
                        /* Use Lyman beta to align background observations. */
            if (objclass == 7 || (program[0] == 'S' && program[2] == '0' && program[3] == '5')) {
	    	printf("Assuming background observation.\n");
                wmin = 1024.;
                wdelta = 6.;
            } else {    /* Use O VI emission for everything else. */
	    	printf("SRC_TYPE = %s.  Emission-line target.\n", src_type);
                wmin = 1030.;
                wdelta = 9.;
            }
	}
					/* Continuum source */
        else if (objclass == 17 || objclass == 29 || objclass == 37) {
                        /* Use O VI and C II to align white dwarf spectra. */
            printf("OBJCLASS = %d.  Assuming nearby white dwarf.\n", objclass);
            wmin = 1030.;
            wdelta = 9.;
        }
        else {  /* Use 1045-1070 A region for all other continuum sources. */
            printf("SRC_TYPE = %s.  Continuum target.\n", src_type);
            wmin = 1045.;
            wdelta = 25.;
        }

	for (i = 0; i < 4; i++) {
	    if (mean_flux[i] > -1) 
		compute_shift(pt_fits, hdu[ref], hdu[i], wmin, wdelta, wpc, shift+i, quality+i);
	    else
		shift[i] = 0;
	    shift[i+4] = shift[i];
	    printf("hdu = %d, hdu_ref = %d, shift = %d\n", hdu[i], hdu[ref], shift[i]);
	}

	/**************************
	** Set up output arrays. **
	**************************/
	n_nvo = cf_nlong((double) (NVO_MAX - NVO_MIN) / wpc) + 1;
	wave_nvo = (float *) cf_calloc(n_nvo, sizeof(float));
	for (i = 0; i < n_nvo; i++) wave_nvo[i] = NVO_MIN + (double) wpc * i;
	flux_nvo = (float *) cf_calloc(n_nvo, sizeof(float));
	error_nvo = (float *) cf_calloc(n_nvo, sizeof(float));
	
	for (i = 0; i < n_nvo; i++) {
		flux_nvo[i]  = 0.0;
		error_nvo[i] = 0.0;
	}
	
	/**************************************************************
	** In each wave band, copy best data set into output arrays. **
	**************************************************************/
	/* For BKGD targets, omit flux comparison when selecting which
	   regions to use. */
	if (objclass == 7 || !strncmp(sp_type, "BKGD", 4)) bkgd_obs = TRUE;

	/* First fill in the 1010-1104 A region with SiC 2B. */
	copy_spec(pt_fits, hdu[1], w_min[1], w_max[1], wpc, shift[1], flux_nvo, error_nvo, n_nvo);

	/* Replace with SiC 1A if possible.  It has a higher S/N, but extends only to 1090 A. */
	if ((bkgd_obs && mean_flux[0] > -1) || (mean_flux[0] > 0.9 * mean_flux[1])) 
	    copy_spec(pt_fits, hdu[0], w_min[0], w_max[0], wpc, shift[0], flux_nvo, error_nvo, n_nvo);

	/* Next overlay the region between 990 and 1080 A with LiF 1A or LiF 2B, if available. */
	if ((bkgd_obs && mean_flux[2] > -1) || ((mean_flux[2] > 0.9 * mean_flux[3]) &&
	    (mean_flux[2] > 0.7 * mean_flux[0]) && (mean_flux[2] > 0.7 * mean_flux[1])))
	    copy_spec(pt_fits, hdu[2], w_min[2], w_max[2], wpc, shift[2], flux_nvo, error_nvo, n_nvo);
	else if ((bkgd_obs && mean_flux[3] > -1) ||
	    ((mean_flux[3] > 0.7 * mean_flux[0]) && (mean_flux[3] > 0.7 * mean_flux[1])))
	    copy_spec(pt_fits, hdu[3], w_min[3], w_max[3], wpc, shift[3], flux_nvo, error_nvo, n_nvo);
	else 
	    w_max[5] = 1010, w_min[6] = w_min[7] = 1104;
		/* If no LiF data, adopt limits of SiC spectra. */

	/* Use SiC 1B to populate the shortest wavelengths. */
	copy_spec(pt_fits, hdu[4], w_min[4], w_max[4], wpc, shift[4], flux_nvo, error_nvo, n_nvo);

	/* If SiC 2A is good, use it for the main part of 900 - 1000 A. */
	if ((bkgd_obs && mean_flux[5] > -1) || (mean_flux[5] > 0.9 * mean_flux[4]))
	    copy_spec(pt_fits, hdu[5], w_min[5], w_max[5], wpc, shift[5], flux_nvo, error_nvo, n_nvo);
	
	/* Use LiF 1B to populate the longest wavelengths. */
	copy_spec(pt_fits, hdu[6], w_min[6], w_max[6], wpc, shift[6], flux_nvo, error_nvo, n_nvo);

	/* If LiF 2A is good, use it for the main part of 1100 - 1180 A. */
	if ((bkgd_obs && mean_flux[7] > -1) || (mean_flux[7] > 0.9 * mean_flux[6]))
	    copy_spec(pt_fits, hdu[7], w_min[7], w_max[7], wpc, shift[7], flux_nvo, error_nvo, n_nvo);
	
	
	for (i = 0; i < n_nvo; i++) {
	    if (flux_nvo[i] >= 1 || error_nvo[i] >= 1) cf_if_error("Bad values in NVO array.");
	}
	
	/*********************
	** Create NVO file. **
	*********************/
	FITS_create_file(&pt_fits_out, out_name, &status);
	FITS_movabs_hdu(pt_fits, 1, &hdutype, &status);
	FITS_copy_header(pt_fits, pt_fits_out, &status);
	FITS_update_key(pt_fits_out, TINT, "NEXTEND", &nextend, NULL, &status);
	FITS_update_key(pt_fits_out, TSTRING, "FILENAME", (out_name + 1), NULL, &status);
	FITS_update_key(pt_fits_out, TSTRING, "FILETYPE", "NVO SPECTRUM", NULL, &status);
	w0 = NVO_MIN;
	FITS_update_key(pt_fits_out, TFLOAT, "W0", &w0, NULL, &status);

	/* Delete index to file extensions from the primary HDU. */
	fits_get_hdrspace(pt_fits_out, &nkeys, &morekeys, &status);
	for (i = nkeys; i > nkeys-12; i--) fits_delete_record(pt_fits_out, i, &status);
	
	FITS_create_tbl(pt_fits_out, BINARY_TBL, n_nvo, 3, ttype, tform, tunit, "FUSE_SPECTRUM", &status);
	FITS_write_col(pt_fits_out, TFLOAT, 1, 1L, 1L, n_nvo, wave_nvo,  &status);
	FITS_write_col(pt_fits_out, TFLOAT, 2, 1L, 1L, n_nvo, flux_nvo,  &status);
	FITS_write_col(pt_fits_out, TFLOAT, 3, 1L, 1L, n_nvo, error_nvo, &status);
	
	free(wave_nvo);
	free(flux_nvo);
	free(error_nvo);
	
	FITS_close_file(pt_fits_out, &status);
	FITS_close_file(pt_fits, &status);
	
	return(0);
}
