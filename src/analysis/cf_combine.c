/*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_combine [-ahk] [-v level] file_list output_idf_file
 *
 * Description: This is a utility program that combines a series of
 *		FUSE spectra in the form of FITS binary tables as produced
 *		by the CalFUSE pipeline.  The files are combined by taking
 *		the exposure-weighted mean of the fluxes and flux errors
 *		and summing the counts.
 *		It is assumed that the error arrays contain standard
 *		deviations, not variances.
 *
 *
 * History:     01/04/00        jwk     started work
 *	        05/01/00 v1.2	jwk	handle SHORT qual arrays
 *	        09/11/00 v1.3	jwk	remove include <sunmath.h>
 *	        10/11/00 v1.4   jwk	update FILENAME keyword in output file
 *              10/27/00 v1.5   peb     #ifdef ieee_retrospective
 *              04/18/01 v1.6   wvd     Change declaration of timeval from
 *                                      long to time_t
 *              10/29/02 v1.7   wvd     Correct error in history line.
 *	        05/19/03 v1.8	wvd	Two updates from jwk:
 *					Close each input file in initial loop;
 *                                      without this cfitsio crashes upon
 *                                      opening the 60th file.
 *	                                Weight each pixel by exposure time
 *                                      after excluding 'bad' pixels.
 *              12/15/03 v1.9	bjg	Changed to work on CalFUSE 3 spectrum
 *                                      data format
 *	        01/12/04 v1.10	bjg	npix_1 is now correct
 *                                      copy header of 1st file instead of 2nd
 *                                      Updated the name of some FITS cols that
 *                                      were still using 2.4 format
 *              03/09/04 v1.11	wvd     Change POTHOLE to QUALITY throughout.
 *              04/05/04   -    bjg     Declarations for cf_wrspec7 and 
 *                                      cf_wrspec_cf2
 *                                      Remove unused variables
 *              06/07/04 v1.12  bjg     Update some keywords (plantime,     
 *                                      rawtime, exp_saa, nbadevnt ...)
 *              04/07/05 v1.13  tc      Read 1D|2D spectrum format. Clean unused
 *                                      variables. Allows one to shift a 2D
 *                                      spectrum (2nd column in the filelist).
 *                                      Use calloc rather than malloc. The
 *                                      quality value is now the nearest short.
 *                                      Update ROOTNAME keyword. Change EXP to
 *                                      OBS if "-k" (optional). Proper use of
 *                                      fclose() and rewind()
 *              04/16/05 v1.14  tc      Use cf_read_col()
 *              05/02/05 v1.15  wvd	Write NSPEC, SPEC###, WOFF### to 
 *					output file header (primary HDU).
 *              06/03/05 v1.16  wvd	Delete unused variables.
 *              07/06/05 v1.17  wvd	Use fgets to read input file.
 *					Change pix_shift to type int.
 *              07/21/05 v1.18  wvd	Adopt the sign convention for spectral
 *					shifts used by FUSE_REGISTER.
 *              08/26/05 v1.19  wvd	Without arguments, don't return an
 *					error.  Just print help and exit.
 *              08/29/05 v1.20  wvd	If input list contains only file 
 *					names, assume shifts are zero.
 *              01/31/06 v1.21  wvd	Revise creation date of output file.
 *              05/19/06 v1.22  wvd	Use same options as idf_combine.
 *					Ignore spectra with EXP_STAT != 0.
 *					Override with -a option.
 *              05/22/06 v1.23  wvd	Rewrite program using a single while
 *					loop.
 *              05/23/06 v1.24  wvd	Properly normalize quality array.
 *              05/26/06 v1.25  wvd	Include unistd.h
 *              07/12/06 v1.26  wvd	Clean up i/o.
 *              08/14/06 v1.27  wvd	Count SPEC### and WAVE### using ngood,
 *					not nfiles.
 *              05/04/07 v1.28  wvd	If the input file list contains a
 *					blank line, skip to the next line.
 *              08/27/07 v1.29  bot	Changed TINT to TLONG l.412
 *					Added L to long l.125 to 136, 351, 378
 *              04/04/08 v1.30  wvd	If -a flag is set, include files with
 *					positive values of EXP_STAT, but not
 *					those with negative values.
 *		07/25/08 v1.31  wvd	Set output EXP_STAT keyword to lowest
 *					non-negative value among input files.
 *
 ****************************************************************************/

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "calfuse.h"

void cf_wrspec7(fitsfile *outfits, long npts, float *wave, float *flux, 
		       float *error, long *counts, float *weights, 
		       float *bkgd,short *quality);
void cf_wrspec_cf2(fitsfile *outfits, int npts, float *wave,
			  float *spec, float *errs, short *qual,
			  float *counts, float *cntserr, int ncol, 
			  int areaflag);

static char CF_PRGM_ID[] = "cf_combine";
static char CF_VER_NUM[] = "1.31";

int main(int argc, char *argv[])
{
    char  *timestr;
    time_t  timeval;
    char  flist[80], infile[80], line[FLEN_CARD], outfile[80], string0[100];
    char  comment[FLEN_CARD];
    char  errstr[80];
    int	  nargs, nfiles, ngood;
    int   status=0;
    int	  hdutype_1, hdutype_o;
    int	  hdunum_1, hdunum_o;
    long  npix_1;
    int	  tfields_1, delta = 0;
    long  j;
    short *quality1;
    float *wave1, *flux1, *error1, *weights1, *bkgd1;
    
    double *flux_o, *error_o, *quality_o;
    float *weights_o, *bkgd_o;
    long *counts_o, *counts1;
    float *exptime_p;
    float f1, shift, wpc;
    
    float exptime_1,  exptime_o = 0.;
    float rawtime_1,  rawtime_o = 0.;
    long plantime_1,  plantime_o = 0L;
    long nbadevnt_1,  nbadevnt_o = 0L;
    long nbadpha_1 ,  nbadpha_o = 0L;
    long exp_bad_1 ,  exp_bad_o = 0L;
    long exp_brst_1,  exp_brst_o = 0L;
    long exp_hv_1,    exp_hv_o = 0L;
    long exp_jitr_1,  exp_jitr_o = 0L;
    long exp_lim_1 ,  exp_lim_o = 0L;
    long exp_saa_1 ,  exp_saa_o = 0L;
    long expnight_1,  expnight_o = 0L;

    long nevents_1, nevents_o = 0L;
    long exp_stat, exp_stat_out=999;

    double	mjd_start, mjd_end;
    double	mjd_start_o = 99999., mjd_end_o = 0.;
    fitsfile *infits_1, *outfits;
    FILE  *ifp;
    char keyword[FLEN_KEYWORD], root_name[32];
	int k_flag = 0;
	
	float *temp_float;
	short *temp_short;
	long *temp_long;

	/* Function prototypes*/
	void getnshift_float(float*, float*, long, int);
	void getnshift_short(short*, short*, long, int);
	void getnshift_long(long*, long*, long, int);	

  int optc;
  int ignore_exp_stat=TRUE;

  char opts[] = "ahkv:";
  char usage[] =
    "Usage:\n"
    "  cf_combine [-ahk]  [-v level] file_list output_idf_file\n";
  char option[] =
    "Options:\n"
    "  -h:  this help message\n"
    "  -v:  verbosity level (=1; 0 is silent)\n"
    "  -a:  ignore EXP_STAT keyword\n"
    "  -k:  Change EXP keywords to OBS\n";

  verbose_level = 1;

  /* Check number of options and arguments */
  while ((optc = getopt(argc, argv, opts)) != -1) {
    switch(optc) {

    case 'h':
      printf("%s\n%s", usage, option);
      return 0;
    case 'a':
      ignore_exp_stat=FALSE;
      break;
    case 'v':
      verbose_level = atoi(optarg);
      break;
    case 'k':
      k_flag=1;
      break;
    }
  }

    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Initialize error checking */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

  if (argc < optind+2) {
    printf("%s", usage);
    return 1;
  }

    strcpy (flist, argv[optind]);
    strcpy (outfile, argv[optind+1]);

    if ((ifp = fopen (flist, "r")) == NULL)
	   {
	   sprintf (errstr, "Unable to open file %s\n", flist);
	   cf_if_error (errstr);
	   }

  if (ignore_exp_stat)
    cf_verbose(1,"Will include only files with EXP_STAT=0") ;
  else
    cf_verbose(1,"Will include only files with EXP_STAT >= 0") ;

    nfiles = ngood = 0;

    /* BIG LOOP THROUGH INPUT FILES */
    while (fgets(line, FLEN_CARD, ifp) != NULL) {
	nargs = sscanf(line, "%s %d", infile, &delta);
	if (nargs == EOF) continue;
	if (nargs == 1) delta = 0;
	nfiles++;
	cf_verbose (3, "Input file loop: nfiles = %d, file = %s", nfiles, infile);

	/* Open input file */
	FITS_open_file(&infits_1, infile, READONLY, &status);
	FITS_read_key(infits_1, TLONG, "EXP_STAT", &exp_stat, NULL, &status);
	FITS_read_key(infits_1, TDOUBLE, "EXPSTART", &mjd_start, NULL, &status);
	FITS_read_key(infits_1, TDOUBLE, "EXPEND"  , &mjd_end  , NULL, &status);
	FITS_read_key(infits_1, TFLOAT, "EXPTIME" , &exptime_1 , NULL, &status);
	FITS_read_key(infits_1, TLONG, "NEVENTS", &nevents_1, NULL, &status);
	FITS_read_key(infits_1, TFLOAT, "RAWTIME" , &rawtime_1 , NULL, &status);
	FITS_read_key(infits_1, TLONG, "PLANTIME", &plantime_1, NULL, &status);
	nevents_o += nevents_1;
	rawtime_o  += rawtime_1;
	plantime_o += plantime_1;

	if (mjd_start < mjd_start_o) mjd_start_o = mjd_start;
	if (mjd_end > mjd_end_o) mjd_end_o = mjd_end;
	
	/* Create output file, copy the header from the input file */
	if (nfiles==1)
	  {
	    /* open output file */
	    remove(outfile);	/* In case the old file name exists */    
	    FITS_create_file (&outfits, outfile, &status);
	    FITS_copy_header(infits_1, outfits, &status);
	    FITS_write_date(outfits, &status);

	    FITS_read_key(infits_1, TFLOAT, "WPC" , &wpc , NULL, &status);
		
	    /* Now go to the second HDU and figure out how long the arrays are. */
	    hdunum_1 = 2;
	    FITS_movabs_hdu(infits_1, hdunum_1, &hdutype_1, &status);
	    FITS_read_key(infits_1, TINT, "TFIELDS", &tfields_1, NULL, &status);

	    if (tfields_1 != 7)
	      {
		sprintf (errstr, "Expecting 7 cols, found %d\n", tfields_1);
		cf_if_error (errstr);
	      }
		
	    npix_1 = cf_read_col(infits_1, TFLOAT, "WAVE", (void **) &wave1);

	    /* Allocate data arrays (calloc initializes to 0 by default). */
	    exptime_p = (float *) cf_calloc(npix_1, sizeof(float));
	    
	    flux1 = (float *) cf_calloc(npix_1, sizeof(float));
	    error1 = (float *) cf_calloc(npix_1, sizeof(float));
	    counts1 = (long *) cf_calloc(npix_1, sizeof(long));
	    weights1 = (float *) cf_calloc(npix_1, sizeof(float));
	    bkgd1 = (float *) cf_calloc(npix_1, sizeof(float));
	    quality1 = (short *) cf_calloc(npix_1, sizeof(short));
		
	    flux_o = (double *) cf_calloc(npix_1, sizeof(double));
	    error_o = (double *) cf_calloc(npix_1, sizeof(double));
	    counts_o = (long *) cf_calloc(npix_1, sizeof(long));
	    weights_o = (float *) cf_calloc(npix_1, sizeof(float));
	    bkgd_o = (float *) cf_calloc(npix_1, sizeof(float));
	    quality_o = (double *) cf_calloc(npix_1, sizeof(double));
		
            hdunum_1 = 1;
            FITS_movabs_hdu(infits_1, hdunum_1, &hdutype_1, &status);

	  }  /* end of if(nfiles==1) */
	
	/* Set EXP_STAT to lowest non-negative value in input files. */
	if (exp_stat_out > exp_stat && exp_stat >= 0) exp_stat_out = exp_stat;

	/* Skip files with EXP_TIME = 0 or EXP_STAT != 0. */
	if (ignore_exp_stat && exp_stat != 0) {
	    cf_if_warning("File %s rejected. EXP_STAT not equal to zero. Use -a flag to override.",
		infile) ;
	    FITS_close_file(infits_1,&status);
	    continue;
	}

	/* Never include files with negative values of EXP_STAT. */
	else if (exp_stat < 0) {
	    cf_if_warning("File %s rejected. EXP_STAT less than zero.", infile) ;
	    FITS_close_file(infits_1,&status);
	    continue;
	}

	if (exptime_1 < 1.) {
	    cf_verbose (3, "Rejecting file %s: EXPTIME = %f", infile, exptime_1);
	    FITS_close_file(infits_1,&status);
	    continue;
	}
	ngood++;
	
	FITS_read_key(infits_1, TLONG, "NBADEVNT", &nbadevnt_1, NULL, &status);
	FITS_read_key(infits_1, TLONG, "NBADPHA" , &nbadpha_1 , NULL, &status);
	FITS_read_key(infits_1, TLONG, "EXP_BAD" , &exp_bad_1 , NULL, &status);
	FITS_read_key(infits_1, TLONG, "EXP_BRST", &exp_brst_1, NULL, &status);
	FITS_read_key(infits_1, TLONG, "EXP_HV",   &exp_hv_1,   NULL, &status);
	FITS_read_key(infits_1, TLONG, "EXP_JITR", &exp_jitr_1, NULL, &status);
	FITS_read_key(infits_1, TLONG, "EXP_LIM" , &exp_lim_1 , NULL, &status);
	FITS_read_key(infits_1, TLONG, "EXP_SAA" , &exp_saa_1 , NULL, &status);
	FITS_read_key(infits_1, TLONG, "EXPNIGHT", &expnight_1, NULL, &status);

	exptime_o  += exptime_1;
	nbadevnt_o += nbadevnt_1;
	nbadpha_o  += nbadpha_1; 
	exp_bad_o  += exp_bad_1;
	exp_saa_o  += exp_saa_1;
	exp_lim_o  += exp_lim_1;
	exp_brst_o += exp_brst_1;
	exp_hv_o   += exp_hv_1;
	exp_jitr_o += exp_jitr_1;
	expnight_o += expnight_1;

	/* Read the data from input file (HDU 2) */
	hdunum_1 = 2;
	FITS_movabs_hdu(infits_1, hdunum_1, &hdutype_1, &status);
	
	(void) cf_read_col(infits_1, TFLOAT, "FLUX", (void **) &temp_float);
	getnshift_float(flux1, temp_float, npix_1, delta);
	free(temp_float);
	
	(void) cf_read_col(infits_1, TFLOAT, "ERROR", (void **) &temp_float);
	getnshift_float(error1, temp_float, npix_1, delta);
	free(temp_float);
	
	(void) cf_read_col(infits_1, TSHORT, "QUALITY", (void **) &temp_short);
	getnshift_short(quality1, temp_short, npix_1, delta);
	free(temp_short);
			
	(void) cf_read_col(infits_1, TLONG, "COUNTS", (void **) &temp_long);
	getnshift_long(counts1, temp_long, npix_1, delta);
	free(temp_long);
			
	(void) cf_read_col(infits_1, TFLOAT, "WEIGHTS", (void **) &temp_float);
	getnshift_float(weights1, temp_float, npix_1, delta);
	free(temp_float);
			
	(void) cf_read_col(infits_1, TFLOAT, "BKGD", (void **) &temp_float);
	getnshift_float(bkgd1, temp_float, npix_1, delta);
	free(temp_float);

        /* Average fluxes, weighting by exposure time; add counts. */
        /* Include only pixels with non-zero quality flags. */
        for (j=0L; j<npix_1; j++) {
	    if (quality1[j] > 0) {
	       f1 = exptime_1;
	       exptime_p[j] += f1;
	       flux_o[j]    += f1 * flux1[j];
	       error_o[j]   += f1 * f1 * (double) error1[j]*error1[j];
	       counts_o[j]  += counts1[j];
	       weights_o[j] += weights1[j];
	       bkgd_o[j]    += bkgd1[j];
	       quality_o[j] += f1 * (double) quality1[j];
	    }
	}
	
	FITS_close_file(infits_1, &status);

	/* Write input file name and spectral shift to output file header. */
	sprintf(keyword, "SPEC%03d", ngood);
	FITS_update_key(outfits, TSTRING, keyword, infile, NULL, &status);
	sprintf(keyword, "WOFF%03d", ngood);
	shift = delta * wpc;
	FITS_update_key(outfits, TFLOAT, keyword, &shift, "[A]", &status);


      }  /* end of loop over input files */

    /* Convert variances back to sigma. Normalize flux, error, and quality by EXPTIME. 
	Convert from doubles to floats. */
    for (j=0L; j<npix_1; j++) if (exptime_p[j] > 0.) {
	flux1[j] = flux_o[j] / exptime_p[j];
	error1[j] = sqrt(error_o[j]) / exptime_p[j];
	quality1[j] = cf_nint(quality_o[j] / exptime_p[j]);
    }
    else {
	flux1[j] = error1[j] = 0.;
	quality1[j] = 0;
    }
    
	/* Write the output file (in multi-row format). */
    cf_wrspec7 (outfits, npix_1, wave1, flux1, error1, counts_o, weights_o, bkgd_o, quality1);
    
    free(exptime_p);
    free(wave1);
    free(flux1);
    free(error1);
    free(bkgd1);
    free(counts1);
    free(weights1);
    free(quality1);
    free(flux_o);
    free(error_o);
    free(quality_o);
    free(bkgd_o);
    free(counts_o);
    free(weights_o);
    
    hdunum_o = 1;
    FITS_movabs_hdu (outfits, hdunum_o, &hdutype_o, &status);

	    /* update output header keywords*/

	    FITS_update_key(outfits, TLONG, "EXP_STAT",  &exp_stat_out, NULL, &status);
	    FITS_update_key(outfits, TDOUBLE, "EXPSTART", &mjd_start_o, NULL, &status);
	    FITS_update_key(outfits, TDOUBLE, "EXPEND", &mjd_end_o  , NULL, &status);
	    FITS_update_key(outfits, TFLOAT, "EXPTIME", &exptime_o , NULL, &status);
	    FITS_update_key(outfits, TLONG, "NEVENTS",  &nevents_o, NULL, &status);
	    FITS_update_key(outfits, TFLOAT, "RAWTIME", &rawtime_o , NULL, &status);
	    FITS_update_key(outfits, TLONG, "PLANTIME", &plantime_o, NULL, &status);
	    FITS_update_key(outfits, TLONG, "NBADEVNT", &nbadevnt_o, NULL, &status);
	    FITS_update_key(outfits, TLONG, "NBADPHA" , &nbadpha_o , NULL, &status);
	    FITS_update_key(outfits, TLONG, "EXP_BAD" , &exp_bad_o , NULL, &status);
	    FITS_update_key(outfits, TLONG, "EXP_BRST", &exp_brst_o, NULL, &status);
	    FITS_update_key(outfits, TLONG, "EXP_HV",   &exp_hv_o, NULL, &status);
	    FITS_update_key(outfits, TLONG, "EXP_JITR", &exp_jitr_o, NULL, &status);
	    FITS_update_key(outfits, TLONG, "EXP_LIM" , &exp_lim_o , NULL, &status);
	    FITS_update_key(outfits, TLONG, "EXP_SAA" , &exp_saa_o , NULL, &status);
	    FITS_update_key(outfits, TLONG, "EXPNIGHT", &expnight_o, NULL, &status);

	    /* update FILENAME keyword in output file */
	    FITS_update_key(outfits, TSTRING, "FILENAME", outfile, NULL, &status);
	    
	    strcpy(string0, "COMBINED SPECTRA");
	    FITS_update_key(outfits, TSTRING, "FILETYPE", string0, NULL, &status);

	    sprintf(string0, "000");
	    FITS_update_key(outfits, TSTRING, "EXP_ID", string0, NULL, &status);
		
	    FITS_read_key(outfits, TSTRING, "ROOTNAME", root_name, NULL, &status);
	    sprintf(string0, "%.8s000", root_name);
	    FITS_update_key(outfits, TSTRING, "ROOTNAME", string0, NULL, &status);

	    FITS_update_key(outfits, TINT, "NSPEC", &ngood,
		"Number of combined spectral files", &status);

	if (k_flag)
	{
	    cf_if_fits_error(fits_modify_name(outfits, "EXP_STAT", "OBS_STAT", &status));
	    cf_if_fits_error(fits_modify_name(outfits, "EXPSTART", "OBSSTART", &status));
	    cf_if_fits_error(fits_modify_name(outfits, "EXPEND", "OBSEND", &status));
	    cf_if_fits_error(fits_modify_name(outfits, "EXPTIME", "OBSTIME", &status));
	    cf_if_fits_error(fits_modify_name(outfits, "EXPNIGHT", "OBSNIGHT", &status));
	}
	else {
	    /* Add HISTORY lines to output file */
	    if (ignore_exp_stat) {
	    sprintf (comment, "cf_combine:  Including only files with EXP_STAT=0");
	    FITS_write_history (outfits, comment, &status);
	    }

	    else {
	    sprintf (comment, "cf_combine:  Including only files with EXP_STAT >= 0");
	    FITS_write_history (outfits, comment, &status);
	    }

	    sprintf (comment, "cf_combine:  File list = %s", flist);
	    FITS_write_history (outfits, comment, &status);
    
	    time (&timeval);
	    timestr = ctime (&timeval);
	    sprintf (comment, "cf_combine:  Completed on %.24s", timestr);
	    FITS_write_history (outfits, comment, &status);
	}
	
    /* Close the files. */
	fclose(ifp);
    FITS_close_file(outfits, &status);
    
    /* Enter a timestamp into the log. */
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
#ifdef SOLARIS
    ieee_retrospective(stdout);
#endif
    return 0;
}


void getnshift_float (float *out_vect, float *in_vect, long n, int delta)
{
	long k, ind;
	
	for (k = 0; k < n; k++)	
	{
		ind = k - delta;
		if((ind >= 0) && (ind < n))
			out_vect[k] = in_vect[ind];
		else
			out_vect[k] = 0;
	}
}

void getnshift_short (short *out_vect, short *in_vect, long n, int delta)
{
	long k, ind;
	
	for (k = 0; k < n; k++)	
	{
		ind = k - delta;
		if((ind >= 0) && (ind < n))
			out_vect[k] = in_vect[ind];
		else
			out_vect[k] = 0;
	}
}

void getnshift_long (long *out_vect, long *in_vect, long n, int delta)
{
	long k, ind;
	
	for (k = 0; k < n; k++)	
	{
		ind = k - delta;
		if((ind >= 0) && (ind < n))
			out_vect[k] = in_vect[ind];
		else
			out_vect[k] = 0;
	}
}
