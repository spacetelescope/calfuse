/****************************************************************************
 *
 * Synopsis:    cf_fpa_position(fitsfile *header, long nevents, float *x,
 *                       unsigned char *channel)
 *
 * Description: Modify X position of photon events to account for an offset
 *              of the focal plane assembly (FPA)
 *
 * Arguments:   fitsfile *header        Pointer to the FITS header
 *                                      of the Intermediate Data File
 *              long     npts           Number of photons in the file
 *              float    *x             X positions for each photon event
 *     unsigned char     *channel       Assigned channel of each photon
 *
 * Calls:       None
 *
 * Return:      0 on success 
 *
 * History:	09/04/02 v1.1	wvd	Based on cf_dpix by HSU
 * 		01/14/03 v1.3	wvd	Rename cf_fpa_pos to cf_read_fpa_pos
 *					Move spectrograph optical parameters
 *					 to separate calibration file.
 *					Abandon linear wavelength scale.
 * 		02/13/03 v1.4	wvd	Write FPA shifts to file header.
 * 		04/07/03 v1.5	wvd	Change keywords to FPADXLIF and
 *					FPADXSIC.  Implement cf_verbose.
 *					Comment out cf_proc_check.
 *					Test channel[i] as int, not char.
 *              05/20/03 v1.6   rdr     Added call to cf_proc_check.
 * 		08/21/03 v1.7	wvd	Change channel to unsigned char.
 * 		08/23/05 v1.8	wvd	If either FPA position is out of 
 *					bounds, issue a warning and set
 *					that FPA correction to zero.
 * 		08/01/05 v1.9	wvd	Broaden FPA limits by +/- 1.
 * 		09/06/05 v1.10  wvd	Correct for both X and Z motions
 *					of FPA.
 *					Delete call to cf_read_fpa_pos.c
 * 		12/02/05 v1.11  wvd	Delete unused variables.
 *              04/07/07 v1.12  wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "calfuse.h"

#define PIX1 2500		/* PIX1,2 are just inside the non-linear */
#define PIX2 13500		/*  regions of the detector 		 */

int cf_fpa_position(fitsfile *header, long nevents, float *x,
	unsigned char *channel)
{
    char CF_PRGM_ID[] = "cf_fpa_position";
    char CF_VER_NUM[] = "1.12";

    int   errflg=0, status=0, hdutype, hdunum;
    int   nwave;     /* number of pixels in calibration file */
    char  aperture[FLEN_VALUE];  /* aperture keyword from FITS header */
    char  detector[FLEN_VALUE];  /* detector keyword from FITS header */
    char  wave_file[FLEN_FILENAME]; /* name of wavelength calibration file */
    char  spec_file[FLEN_FILENAME]; /* name of spectrograph calibration file */
    long  i;
			/* Spectrograph optical parameters */
    float alpha, alpha_lif, alpha_sic, sigma_lif, sigma_sic, diam;
    float *wavelength;  /* wavelengths read from calibration file */
    float fpaxlifdata;	/* LiF FPA X position in microns for spectrum */
    float fpaxsicdata;	/* SiC FPA X position in microns for spectrum */
    float fpaxcal;	/* FPA X position in microns for calibration file */
    float dxfpa;	/* fpax_data - fpaxcal */
    float fpazlifdata;	/* LiF FPA Z position in microns for spectrum */
    float fpazsicdata;	/* SiC FPA Z position in microns for spectrum */
    float fpazcal;	/* FPA Z position in microns for calibration file */
    float dzfpa;	/* fpaz_data - fpazcal */
    float dpix, dpix_LiF, dpix_SiC; /* Spectral shift (in pixels)*/

    double w1, w2;	/* wavelengths at which to sample calibration */
    double beta1,beta2;	/* grating exit angle for w1, w2 */
    double sinalpha;	/* sin (spectrograph entrance angle) */
    double cosalpha;	/* cos (spectrograph entrance angle) */
    double sigma;	/* grating groove spacing in Angstroms */
    double cosbeta;	/* cos (mean grating exit angle) */
    double pixscale;	/* mean effective microns/pixel for det segment */
    fitsfile *specfits, *wavefits;

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Check whether routine is appropriate for this data file. */
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Read keywords from data file header. */
    FITS_read_key(header, TSTRING, "DETECTOR", detector, NULL, &status);
    FITS_read_key(header, TSTRING, "APERTURE", aperture, NULL, &status);
    FITS_read_key(header, TFLOAT,  "FPASXPOS", &fpaxsicdata, NULL, &status);
    FITS_read_key(header, TFLOAT,  "FPALXPOS", &fpaxlifdata, NULL, &status);
    FITS_read_key(header, TFLOAT,  "FPASZPOS", &fpazsicdata, NULL, &status);
    FITS_read_key(header, TFLOAT,  "FPALZPOS", &fpazlifdata, NULL, &status);
    FITS_read_key(header, TSTRING, "SPEC_CAL", spec_file, NULL, &status);
    FITS_read_key(header, TSTRING, "WAVE_CAL", wave_file, NULL, &status);

    /* Read spectrograph optical parameters from SPEC_CAL file. */
    FITS_open_file(&specfits, cf_cal_file(spec_file), READONLY, &status);
    FITS_read_key(specfits, TFLOAT, "ALPHALIF", &alpha_lif, NULL, &status);
    FITS_read_key(specfits, TFLOAT, "ALPHASIC", &alpha_sic, NULL, &status);
    FITS_read_key(specfits, TFLOAT, "SIGMALIF", &sigma_lif, NULL, &status);
    FITS_read_key(specfits, TFLOAT, "SIGMASIC", &sigma_sic, NULL, &status);
    FITS_read_key(specfits, TFLOAT, "DIAM", &diam, NULL, &status);
    FITS_close_file(specfits, &status);

    /* Determine the extension number (LiF) from aperture */
    hdunum = 4;					/* Default is LWRS */
    if (!strcmp(aperture, "HIRS")) hdunum = 2;
    if (!strcmp(aperture, "MDRS")) hdunum = 3;

    /* Open wavelength calibration file */
    FITS_open_file(&wavefits, cf_cal_file(wave_file), READONLY, &status);

    /* Loop over the two apertures; first (i==0) LiF, then SiC. */
    /* Determine shifts for LiF and SiC channels.		*/
    for(i=0; i<2; i++) {

	/* Move to appropriate HDU in WAVECAL file. */
	hdunum += i*4;
	FITS_movabs_hdu(wavefits, hdunum, &hdutype, &status);

	/* Get FPA position for which WAVECAL was derived */
	FITS_read_key(wavefits, TFLOAT, "FPACXPOS", &fpaxcal, NULL, &status);
	FITS_read_key(wavefits, TFLOAT, "FPACZPOS", &fpazcal, NULL, &status);

        /* Read the wavelength array from the WAVECAL file. */
	nwave=cf_read_col(wavefits, TFLOAT,"WAVELENGTH", (void **) &wavelength);

	/* Set LiF/SIC-dependent quantities */
	if (i == 1) { 			/* SiC data */
	    dxfpa = fpaxsicdata - fpaxcal;
	    dzfpa = fpazsicdata - fpazcal;
	    sigma = sigma_sic;
	    alpha = alpha_sic;
	} else { 			/* LiF data */
	    dxfpa = fpaxlifdata - fpaxcal;
	    dzfpa = fpazlifdata - fpazcal;
	    sigma = sigma_lif;
	    alpha = alpha_lif;
	}
	/*
	 * Compute mean pixel scale from points just inside the
	 * highly-nonlinear regions at edges of detectors */
        w1 = wavelength[PIX1];
        w2 = wavelength[PIX2];
	sinalpha = sin (alpha * RADIAN);
	cosalpha = cos (alpha * RADIAN);
	beta1 = asin (w1/sigma - sinalpha);
	beta2 = asin (w2/sigma - sinalpha);
	cosbeta = cos ((beta1 + beta2)/2.);  /* mean cosbeta for segment */
	pixscale = diam * fabs(beta1 - beta2) * 1000. / (double)(PIX2-PIX1);

	/* Compute shift of spectrum in pixels due to the FPA offset. */
	if (i == 1)			/* SiC data */
	    dpix = (dxfpa * cosalpha + dzfpa * sinalpha) / cosbeta / pixscale;

	else				/* LiF data */
	    dpix = (dxfpa * cosalpha - dzfpa * sinalpha) / cosbeta / pixscale;

	/* For side 2, multiply the above expressions by -1. */
	if (!strncmp(detector, "2", 1))
	    dpix *= -1;

	/* Write info to log file */
	if (i == 1) {
	   dpix_SiC = dpix;
	   cf_verbose(3, "FPAX data: %8.3f   FPAX calib: %8.3f",
		   fpaxsicdata, fpaxcal);
	   cf_verbose(3, "FPAZ data: %8.3f   FPAZ calib: %8.3f",
		   fpazsicdata, fpazcal);
	   cf_verbose(2, "SiC spectra to be shifted by %8.3f pixels", dpix_SiC);
	}
	else {
	   dpix_LiF = dpix;
	   cf_verbose(3, "FPAX data: %8.3f   FPAX calib: %8.3f",
		   fpaxlifdata, fpaxcal);
	   cf_verbose(3, "FPAZ data: %8.3f   FPAZ calib: %8.3f",
		   fpazlifdata, fpazcal);
	   cf_verbose(2, "LiF spectra to be shifted by %8.3f pixels", dpix_LiF);
	}
  
    } /* End of i loop */

    /* Close the WAVECAL file */
    FITS_close_file(wavefits, &status);

    /* If tabulated FPA position is out of range, set the corresponding
	dpix value to zero and issue a warning. */
    if (fpaxlifdata < -1 || fpaxlifdata > 401) {
	dpix_LiF = 0;
	cf_if_warning("Keyword FPALXPOS is out of bounds.  Setting FPADXLIF = 0.");
    }
    if (fpaxsicdata < -1 || fpaxsicdata > 401) {
	dpix_SiC = 0;
	cf_if_warning("Keyword FPASXPOS is out of bounds.  Setting FPADXSIC = 0.");
    }

    /* Shift only photons with assigned apertures. */
    for(i = 0; i < nevents; i++) {
	switch(channel[i]) {
	case 1: case 2: case 3: case 4:
	    x[i] += dpix_LiF;
	    break;
	case 5: case 6: case 7: case 8:
	    x[i] += dpix_SiC;
	    break;
	default:
	    break;
	}
    }

    /* Write shifts to IDF file header */
    FITS_update_key(header, TFLOAT, "FPADXLIF", &dpix_LiF, NULL, &status);
    FITS_update_key(header, TFLOAT, "FPADXSIC", &dpix_SiC, NULL, &status);

    /* Update processing flags. */
    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");

    return status;
}
