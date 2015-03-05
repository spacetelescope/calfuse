/****************************************************************************
 *	        Johns Hopkins University
 *	        Center For Astrophysical Sciences
 *	        FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_grating_motion(header, nevents, time, x, y, channel,
 *			ntimes, ttime, tsunrise)
 *
 * Description: Reads spectral shift caused by thermal grating motions
 *              and corrects the X and Y locations of each photon event.
 *
 *              fitsfile  *header       Pointer to the FITS header of the
 *                                      Intermediate Data File
 *              long      nevents       Number of photons in the file
 *              int       *time         Time (in seconds) since exposure start
 *              float     *x            Position of the photon (in pixels)
 *              float     *y            Position of the photon (in pixels)
 *     unsigned char      channel       Channel ID of each photon
 *		long	  ntimes	Number of entries in timeline table
 *		float	  ttime		Time array of timeline table
 *		short	  tsunrise	Time since last sunrise
 *
 * Returns:     0 on success
 *
 * History:	09/05/02    1.0  RDR    Begin work; adapted from cf_make_shift
 *		04/21/03    1.3  wvd	Use tsunrise array from timeline table.
 *					Do not assume that ttime is continuous.
 *              05/20/03    1.4  rdr    Added call to cf_proc_check
 *              08/21/03    1.5  wvd    Change channel to unsigned char.
 *              06/25/04    1.6  wvd    Estimate time between sunrises using
 *                                      orbit period in file header.
 *              03/22/05    1.7  wvd    Change tsunrise from float to short.
 *                                      Read orbital period from file header.
 *              04/15/05    1.8  wvd    Close GRAT_CAL file before exiting.
 *					If ORBPERID keyword is not found,
 *					assume that orbit period is 6000 s.
 *		12/30/05    1.9  wvd	Current algorithm corrects for orbital
 *					drift of spectrum, but wavelength zero
 *					point still varies with beta angle and
 *					observation date.  New subroutine
 *					corrects LiF 1 data for these effects.
 *		07/06/06    2.0  wvd	Rewrite to use Dave Sahnow's new
 *					grating-correction file (version 003).
 *		07/12/06    2.1  wvd	Modify to read an arbitrary number of
 *					extensions, each with phase information
 *					in file header keywords.
 *		10/20/06    2.2  wvd	Grating motion also depends on APER_PA.
 *					Now we read it from the file header
 *					and select appropriate coefficients.
 *					If no correction is defined for a 
 *					particular set of parameters, set
 *		11/28/06    2.3  wvd	Grating motion includes a long-term,
 *					non-periodic drift.  Read it from the
 *					first four extensions of GRAT_CAL file.
 *              04/07/07    2.4  wvd	Clean up compiler warnings.
 *
 ****************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "calfuse.h"

static int
read_nonperiodic_shifts(fitsfile *header, fitsfile *shiftfits, 
    long ntimes, float *dxlif, float *dylif, float *dxsic, float *dysic)
{
    char  channel[FLEN_VALUE], comment[FLEN_CARD], detector[FLEN_VALUE];
    int   nrows, tndx;
    int   status=0, hdutype=0, hdunum=0;
    long  i;
    float *mjd, *xshift, *yshift;
    double expstart;
    
    /* Read header keywords from data file. */
    fits_read_key(header, TDOUBLE, "EXPSTART", &expstart, NULL, &status);
    FITS_read_key(header, TSTRING, "DETECTOR", detector, NULL, &status);

    /* First do the LiF channels */
    if (!strncmp(detector,"1",1)) {
         hdunum=2;
    } else if (!strncmp(detector,"2",1)) {
         hdunum=4;
    } else {
         cf_if_error("Cannot parse DETECTOR keyword in data file.");
    }
    cf_verbose(3, "Reading HDU %d of GRAT_CAL file.", hdunum);

    FITS_movabs_hdu(shiftfits, hdunum, &hdutype, &status);
    FITS_read_key(shiftfits, TSTRING, "CHANNEL", channel, NULL, &status);
    sprintf(comment, "LiF%c", detector[0]);
    if (strncmp(channel, comment, 4))
        cf_if_error("Cannot find %s correction info in extension %d",
	comment, hdunum);

    /* Read X and Y corrections as a function of time. */
    nrows = cf_read_col(shiftfits, TFLOAT, "MJD", (void *) &mjd);
    nrows = cf_read_col(shiftfits, TFLOAT, "XSHIFT", (void *) &xshift);
    nrows = cf_read_col(shiftfits, TFLOAT, "YSHIFT", (void *) &yshift);

    /* Select appropriate MJD for this exposure. */
    i = 0;
    while (expstart > mjd[i] && i < nrows) i++;
    tndx = i-1;
    if (tndx < 0) tndx = 0;
    cf_verbose(3,"EXPSTART = %g", expstart);
    cf_verbose(3,"Using LiF corrections for MJD >= %g", mjd[tndx]);
    if (tndx < nrows-1)
        cf_verbose(3,"                      and MJD <  %g", mjd[tndx+1]);
    cf_verbose(3,"xshift = %f, yshift = %f", xshift[tndx], yshift[tndx]);

    /* Copy X and Y corrections to output arrays. */
    for (i = 0; i < ntimes; i++) {
	dxlif[i] = xshift[tndx];
	dylif[i] = yshift[tndx];
    }

    /* Now do the SiC channels */
    hdunum += 1;
    cf_verbose(3, "Reading HDU %d of GRAT_CAL file.", hdunum);
    FITS_movabs_hdu(shiftfits, hdunum, &hdutype, &status);
    FITS_read_key(shiftfits, TSTRING, "CHANNEL", channel, NULL, &status);
    sprintf(comment, "SiC%c", detector[0]);
    if (strncmp(channel, comment, 4))
        cf_if_error("Cannot find %s correction info in extension %d",
	comment, hdunum);

    /* Read X and Y corrections as a function of time. */
    nrows = cf_read_col(shiftfits, TFLOAT, "MJD", (void *) &mjd);
    nrows = cf_read_col(shiftfits, TFLOAT, "XSHIFT", (void *) &xshift);
    nrows = cf_read_col(shiftfits, TFLOAT, "YSHIFT", (void *) &yshift);

    /* Select appropriate MJD for this exposure. */
    i = 0;
    while (expstart > mjd[i] && i < nrows) i++;
    tndx = i-1;
    if (tndx < 0) tndx = 0;
    cf_verbose(3,"Using SiC corrections for MJD >= %g", mjd[tndx]);
    if (tndx < nrows-1)
        cf_verbose(3,"                      and MJD <  %g", mjd[tndx+1]);
    cf_verbose(3,"xshift = %f, yshift = %f", xshift[tndx], yshift[tndx]);

    /* Copy X and Y corrections to output arrays. */
    for (i = 0; i < ntimes; i++) {
	dxsic[i] = xshift[tndx];
	dysic[i] = yshift[tndx];
    }

    return 0;
}


static int
read_periodic_shifts(fitsfile *header, fitsfile *shiftfits, int ncorrection,
    double *lif_mjd0, double *lif_period, double *sic_mjd0, double *sic_period,
    double *lif_x_coef, double *lif_y_coef, double *sic_x_coef, double *sic_y_coef)
{
    char  detector[FLEN_CARD], channel[FLEN_CARD], comment[FLEN_CARD];
    int   status=0, hdutype=0, hdunum=0, anynul=0;
    int   i, nregions, region, typecode;
    long  max_coeff, ncoeff, width;
    float aper_pa, *aperpalo, *aperpahi;
    float *betalo, *betahi, *polelo, *polehi;
    double beta, pole, sunang;
    
    /* Initialize variables. */
    *lif_mjd0 = *lif_period = *sic_mjd0 = *sic_period = 0.;
    for (i = 0; i < 11; i++) {
	lif_x_coef[i] = 0.;
	lif_y_coef[i] = 0.;
	sic_x_coef[i] = 0.;
	sic_y_coef[i] = 0.;
    }

    FITS_read_key(header, TFLOAT,  "APER_PA", &aper_pa, NULL, &status);
    FITS_read_key(header, TDOUBLE, "SUNANGLE", &sunang, NULL, &status);
    FITS_read_key(header, TDOUBLE, "POLEANGL", &pole, NULL, &status);
    FITS_read_key(header, TSTRING, "DETECTOR", detector, NULL, &status);
    while (999.9 > aper_pa && aper_pa > 360.) aper_pa -= 360.;
    while (aper_pa < 0.) aper_pa += 360.;
    beta = 180.0 - sunang;

    cf_verbose(3, "S/C orientation: beta = %5.1f, pole = %5.1f, aper_pa = %5.1f",
	beta, pole, aper_pa);

    /* First do the LiF channels */
    if (!strncmp(detector,"1",1)) {
         hdunum=2;
    } else if (!strncmp(detector,"2",1)) {
         hdunum=4;
    } else {
         cf_if_error("Cannot parse DETECTOR keyword in read_shift_file");
    }
    hdunum += ncorrection * 4;
    cf_verbose(3, "Reading HDU %d of GRAT_CAL file.", hdunum);

    FITS_movabs_hdu(shiftfits, hdunum, &hdutype, &status);
    FITS_read_key(shiftfits, TSTRING, "CHANNEL", channel, NULL, &status);
    sprintf(comment, "LiF%c", detector[0]);
    if (strncmp(channel, comment, 4))
        cf_if_error("Cannot find %s correction info in extension %d",
	comment, hdunum);

    if (ncorrection > 1) {
	FITS_read_key(shiftfits, TDOUBLE, "PERIOD", lif_period, NULL, &status);
	FITS_read_key(shiftfits, TDOUBLE, "MJD_ZERO", lif_mjd0, NULL, &status);
	cf_verbose(3, "%s: period = %g, mjd_zero = %g",
		channel, *lif_period, *lif_mjd0);
    }
    
    nregions = cf_read_col(shiftfits, TFLOAT, "BETALO", (void *) &betalo);
    nregions = cf_read_col(shiftfits, TFLOAT, "BETAHI", (void *) &betahi);
    nregions = cf_read_col(shiftfits, TFLOAT, "POLELO", (void *) &polelo);
    nregions = cf_read_col(shiftfits, TFLOAT, "POLEHI", (void *) &polehi);
    nregions = cf_read_col(shiftfits, TFLOAT, "APERPALO", (void *) &aperpalo);
    nregions = cf_read_col(shiftfits, TFLOAT, "APERPAHI", (void *) &aperpahi);
    FITS_get_coltype(shiftfits, 7, &typecode, &ncoeff, &width, &status);
    max_coeff = ncoeff;

    for (i = 0; i < nregions; i++) 
	if (betalo[i] <= beta && beta < betahi[i] &&
	    polelo[i] <= pole && pole < polehi[i] &&
	    aperpalo[i] <= aper_pa && aper_pa < aperpahi[i]) break;

    if (i == nregions) 
	cf_verbose(3, "LiF correction not defined for beta = %g, pole = %g, "
	"aper_pa = %g", beta, pole, aper_pa);
    else {
	cf_verbose(3, "LIF REGION %2d: beta = %g - %g, pole = %g - %g,", 
	i, betalo[i], betahi[i], polelo[i], polehi[i]);
	cf_verbose(3, "                aper_pa = %g - %g, %d coefficients",
	aperpalo[i], aperpahi[i], ncoeff);

	region = i+1; 		/* CFITSIO counts from 1, not 0. */
	FITS_read_col(shiftfits, TDOUBLE, 7, region, 1, ncoeff, NULL, lif_x_coef,
	&anynul, &status);
	FITS_read_col(shiftfits, TDOUBLE, 8, region, 1, ncoeff, NULL, lif_y_coef,
	&anynul, &status);
    }

    /* Now the SiC channels */
    hdunum += 1;
    cf_verbose(3, "Reading HDU %d of GRAT_CAL file.", hdunum);
    FITS_movabs_hdu(shiftfits, hdunum, &hdutype, &status);
    FITS_read_key(shiftfits, TSTRING, "CHANNEL", channel, NULL, &status);
    sprintf(comment, "SiC%c", detector[0]);
    if (strncmp(channel, comment, 4))
        cf_if_error("Cannot find %s correction info in extension %d", 
	comment, hdunum);

    if (ncorrection > 1) {
	FITS_read_key(shiftfits, TDOUBLE, "PERIOD", sic_period, NULL, &status);
	FITS_read_key(shiftfits, TDOUBLE, "MJD_ZERO", sic_mjd0, NULL, &status);
	cf_verbose(3, "%s: period = %g, mjd_zero = %g",
		channel, *sic_period, *sic_mjd0);
    }
    
    nregions = cf_read_col(shiftfits, TFLOAT, "BETALO", (void *) &betalo);
    nregions = cf_read_col(shiftfits, TFLOAT, "BETAHI", (void *) &betahi);
    nregions = cf_read_col(shiftfits, TFLOAT, "POLELO", (void *) &polelo);
    nregions = cf_read_col(shiftfits, TFLOAT, "POLEHI", (void *) &polehi);
    nregions = cf_read_col(shiftfits, TFLOAT, "APERPALO", (void *) &aperpalo);
    nregions = cf_read_col(shiftfits, TFLOAT, "APERPAHI", (void *) &aperpahi);
    FITS_get_coltype(shiftfits, 7, &typecode, &ncoeff, &width, &status);
    if (max_coeff < ncoeff) max_coeff = ncoeff;

    for (i = 0; i < nregions; i++)
        if (betalo[i] <= beta && beta < betahi[i] &&
            polelo[i] <= pole && pole < polehi[i] &&
	    aperpalo[i] <= aper_pa && aper_pa < aperpahi[i]) break;

    if (i == nregions) 
	cf_verbose(3, "SiC correction not defined for beta = %g, pole = %g, "
	"aper_pa = %g", beta, pole, aper_pa);
    else {
	cf_verbose(3, "SIC REGION %2d: beta = %g - %g, pole = %g - %g,",
	i, betalo[i], betahi[i], polelo[i], polehi[i]);
	cf_verbose(3, "                aper_pa = %g - %g, %d coefficients",
	aperpalo[i], aperpahi[i], ncoeff);

	region = i+1; 		/* CFITSIO counts from 1, not 0. */
	FITS_read_col(shiftfits, TDOUBLE, 7, region, 1, ncoeff, NULL, sic_x_coef,
	    &anynul, &status);
	FITS_read_col(shiftfits, TDOUBLE, 8, region, 1, ncoeff, NULL, sic_y_coef,
	    &anynul, &status);
    }

    cf_verbose(2, "Fourier coefficients of grating-motion correction:");
    cf_verbose(2,"    LiF X      LiF Y      SiC X      SiC Y");
    for (i=0; i<max_coeff; i++) {
	cf_verbose(2, "%10.5f %10.5f %10.5f %10.5f",lif_x_coef[i],
	       lif_y_coef[i], sic_x_coef[i], sic_y_coef[i]);
    }

    return 0;
}


static double
thermal_shift(double theta, double *a)
{
    theta *= 2.0 * PI;
    return      (a[0]+
		 a[1]*sin(1.0*theta)+a[2]*cos(1.0*theta)+
		 a[3]*sin(2.0*theta)+a[4]*cos(2.0*theta)+
		 a[5]*sin(3.0*theta)+a[6]*cos(3.0*theta)+
		 a[7]*sin(4.0*theta)+a[8]*cos(4.0*theta)+
		 a[9]*sin(5.0*theta)+a[10]*cos(5.0*theta));
}


int
cf_grating_motion(fitsfile *header, long nevents, float *time, float *x,
	float *y, unsigned char *channel, long ntimes, float *ttime,
	short *tsunrise)
{
    char CF_PRGM_ID[] = "cf_grating_motion";
    char CF_VER_NUM[] = "2.4";

    char     shift_file[FLEN_CARD];
    int      errflg=0, status=0;
    int      calfvers, num_hdus, ncorrections;
    long     i, j, k;
    float    *dxlif, *dylif, *dxsic, *dysic, period;
    double   lif_mjd0, lif_period, sic_mjd0, sic_period;
    double   expstart, phase, lif_phase, sic_phase;
    double   lif_x_coef[11], lif_y_coef[11], sic_x_coef[11], sic_y_coef[11];
    fitsfile *shiftfits;

    /* Initialize error checking. */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /* Enter a time stamp into the log */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Check whether routine is appropriate for this data file. */
    if ((errflg = cf_proc_check(header, CF_PRGM_ID))) return errflg;

    /* Read exposure start time from IDF header. */
    fits_read_key(header, TDOUBLE, "EXPSTART", &expstart, NULL, &status);

    /* Read orbital period from IDF header. If not found, estimate. */
    fits_read_key(header, TFLOAT, "ORBPERID", &period, NULL, &status);
    if (status) {
	status = 0;
	period = 6000;
	cf_verbose(1, "Keyword ORBPERID not found; assuming 6000 seconds.");
    }
    else
	cf_verbose(2, "Estimated orbital period is %d seconds", cf_nint(period));

    /* Allocate space for shift arrays. */
    dxlif = (float *) cf_calloc(ntimes, sizeof(float));
    dylif = (float *) cf_calloc(ntimes, sizeof(float));
    dxsic = (float *) cf_calloc(ntimes, sizeof(float));
    dysic = (float *) cf_calloc(ntimes, sizeof(float));

    /* Open the GRAT_CAL file, check version number, and count extensions. */
    FITS_read_key(header, TSTRING, "GRAT_CAL", shift_file, NULL, &status);
    FITS_open_file(&shiftfits, cf_cal_file(shift_file), READONLY, &status);
    FITS_get_num_hdus(shiftfits, &num_hdus, &status);
    cf_verbose(3, "GRAT_CAL file %s contains %d HDUs.", shift_file, num_hdus);
    ncorrections = num_hdus / 4;

    /* Exit if the CALFVERS keyword is less than 4. */
    FITS_read_key(shiftfits, TINT, "CALFVERS", &calfvers, NULL, &status);
    if (calfvers < 4) 
        cf_if_error("Grating-calibration file format is out of date.");

    /* First correct for the long-term, non-periodic grating motions. */
    read_nonperiodic_shifts(header, shiftfits, ntimes, dxlif, dylif, dxsic, dysic);

    /* Next correct for periodic motions. */
    for (i = 1; i < ncorrections; i++) {

	/* Read the calibration file.  If there's no correction for this spacecraft
	orientation, the coefficients will be set to zero. */
	read_periodic_shifts(header, shiftfits, i, &lif_mjd0, &lif_period, &sic_mjd0,
	    &sic_period, lif_x_coef, lif_y_coef, sic_x_coef, sic_y_coef);

	/* For orbital motions, the period is just the orbit period. */
	if (i == 1) {
	    for (k=0; k<ntimes; k++) {
		phase = tsunrise[k]/period;
		dxlif[k] += thermal_shift(phase, lif_x_coef);
		dylif[k] += thermal_shift(phase, lif_y_coef);
		dxsic[k] += thermal_shift(phase, sic_x_coef);
		dysic[k] += thermal_shift(phase, sic_y_coef);
	    }
	}
	/* For longer-term motions, read the period and zero point from the file header. */
	else {
	    for (k=0; k<ntimes; k++) {
		lif_phase = ((expstart - lif_mjd0) + ttime[k]/(24.*3600.))/lif_period;
		sic_phase = ((expstart - sic_mjd0) + ttime[k]/(24.*3600.))/sic_period;
		dxlif[k] += thermal_shift(lif_phase, lif_x_coef);
		dylif[k] += thermal_shift(lif_phase, lif_y_coef);
		dxsic[k] += thermal_shift(sic_phase, sic_x_coef);
		dysic[k] += thermal_shift(sic_phase, sic_y_coef);
	    }
	}
    }

    /* Close the GRAT_CAL file. */
    FITS_close_file(shiftfits, &status);

    /* Apply grating-motion correction. */
    for (j=k=0; j<nevents; j++) {
        while(ttime[k+1] - FRAME_TOLERANCE < time[j] && k+1 < ntimes) k++;
	if (channel[j] > 0 && channel[j] < 4) {
	    x[j] += dxlif[k];
	    y[j] += dylif[k];
	}
	else if (channel[j] > 4) {
	    x[j] += dxsic[k];
	    y[j] += dysic[k];
	} 
    }
    
    /* Release memory. */
    free(dxlif);
    free(dylif);
    free(dxsic);
    free(dysic);

    cf_proc_update(header, CF_PRGM_ID, "COMPLETE");
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");

    return status;
}
