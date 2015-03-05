/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences 
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_extract_spectra options intermediate_file  
 *
 * Description: This module computes a background model, reads the various
 *		calibration files, and performs either a standard or optimal
 *		extraction for the target spectrum in each of the LiF and
 *		SiC channels.
 *
 *		Note: some of the routines called by this module modify values
 *		in the IDF file header.  We want those changes to appear in
 *		the final, extracted spectral files, but not in the archived
 *		IDF, so we make a copy in memory of the IDF header and pass it
 *		to the various routines.  We delete it at the end.
 *
 * Arguments:   input_file              Intermediate data file
 *
 * Calibration files:                   
 *
 * Returns:     0 on successful completion
 *
 * Calls:       cf_apply_filters, cf_scale_bkgd, cf_make_wave_array,
 *		cf_rebin_and_flux_calibrate, cf_rebin_probability_array,
 *              cf_standard_or_optimal_extraction, cf_optimal_extraction
 *              cf_write_extracted_spectrum
 *
 * History:     02/13/03   peb    1.1   Begin work
 *		03/02/03   wvd    1.3	Add cf_standard_or_optimal_extraction()
 *              03/10/03   peb    1.4   Added verbose_level, unistd.h for
 *                                      getopt portability, and remove debug
 *                                      print statements
 *		04/08/03   wvd    1.6	Copy header of IDF into memory
 *					so that subsequent routines can
 *					modify it without changing IDF.
 *					Delete this virtual header when done.
 *					Read arrays from timeline table and
 *					pass to cf_apply_filters.
 *              05/06/03   rdr    1.7   Read from column ERGCM2 rather than
 *                                      ERGCM2S
 *              05/07/03   wvd    1.8   Change ergcm2 to ergcm2 throughout
 *              05/28/03   rdr    1.9   Read pothole mask
 *              06/09/03   rdr    1.10  Incorporate the tscreen flag
 *		06/11/03   wvd    1.11  Pass datatype to cf_read_col
 *		08/25/03   wvd    1.12 	Change coltype from string to int in
 *					cf_read_col.
 *		09/29/03   wvd    1.13 	Move standard_or_optimal so that it
 *					runs just once.  Don't check return
 *					values from subroutines.
 *		10/08/03   wvd    1.14 	Change counts_out array to type long.
 *		10/16/03   wvd    1.15 	If EXPTIME = 0, generate a pair of
 *					null-valued spectra and exit.
 *		10/31/03   wvd    1.16 	Change channel to unsigned char.
 *		12/12/03   wvd    1.18 	Clean up i/o.
 *		03/16/04   wvd    1.19 	Don't pass wave_out to optimal
 *					extraction subroutine.
 *					Change pycent from int to float.
 *		03/24/04   wvd    1.20 	Eliminate got_no_data(); fold into
 *					main routine.  If either EXPTIME = 0
 *					or NEVENTS = 0, write null spectrum.
 *		04/26/04   wvd    1.21 	Replace cf_rebin_and_flux_calibrate
 *					background with cf_rebin_background.
 *					cf_optimal_extraction now returns
 *					flux_out and sigma_out in units of
 *					counts; must flux-calibrate each.
 *		06/10/04   wvd    1.22 	Don't pass time array from timeline
 *					table to cf_apply_filters.
 *              06/11/04   peb    1.23  Add the -r option to override the
 *                                      default rootname.
 *		01/28/05   wvd    1.24	If EXP_STAT < 0, generate a pair of
 *                                      null-valued spectra and exit.
 *		03/09/05   wvd    1.25	Change cf_ttag_bkgd to cf_scale_bkgd,
 *					pass weights array to it.
 *		04/12/05   wvd    1.27	Add -n flag, which forces extraction
 *					of night-only spectrum.  Its argument
 *					is the name of a night-only BPM file.
 *		06/15/05   wvd    1.28  BUG FIX: program always read the
 *					point-source probability arrays from
 *					the WGTS_CAL file.  Now uses value
 *					of extended to determine which HDU
 *					to read.
 *		11/23/05   wvd    1.29  Add flag to disable optimal 
 *					extraction from the command line.
 *		01/27/06   wvd    1.30  Add flag to force use of optimal 
 *					extraction from the command line.
 *		05/19/06   wvd    1.31	Don't discard spectra with non-zero
 *					values of EXP_STAT.
 *		06/12/06   wvd    1.32	Change cf_if_warning to cf_verbose.
 *
 ****************************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_extract_spectra";
static char CF_VER_NUM[] = "1.32";

int
main(int argc, char *argv[])
{
    char  *bpm_file=NULL, *outrootname = NULL;
    unsigned char *channel=NULL, *timeflags=NULL, *loc_flags=NULL;
    unsigned char *time_status=NULL;
    int   extended, status=0, optc, j, tscreen=1;
    int   binx, biny, bnx[2], bny[2], bymin[2], aper[2];
    int   night_only=FALSE, optimal=FALSE;
    int   force_optimal=FALSE, override_optimal=FALSE;
    long  nevents=0, ntimes=0, ngood=0, dtime=0, ntime=0, *good_index=NULL;
    float *weight=NULL, *x=NULL, *y=NULL, *lambda=NULL;
    float exptime, *bpmask=NULL ;
    float *bimage[2]={NULL, NULL};

    fitsfile *memp, *header;
    
    char opts[]  = "hn:ofr:sv:";
    char usage[] = 
	"Usage:\n"
	"  cf_extract_spectra [-hsof] [-n bpm_filename] [-r rootname] [-v level]"
		" idf_file\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
	"  -n:  extract night-only spectrum using given bad-pixel map\n"
	"  -o:  disable optimal extraction\n"
	"  -f:  force optimal extraction\n"
        "  -r:  override the default rootname\n"
        "  -s:  do not perform screening on time flags\n" 
        "  -v:  verbosity level (=1; 0 is silent)\n" ;

    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
        case 'n':
            bpm_file = optarg;
	    night_only = TRUE;
            break;
        case 'o':
	    override_optimal = TRUE;
            break;
        case 'f':
	    force_optimal = TRUE;
            break;
        case 'r':
            outrootname = optarg;
            break;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
        case 's':
	  tscreen = 0 ;
          break;
	}
    }
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    if (argc <= optind) {
        printf("%s", usage);
	return -1;
    }
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    /* Open input data file. */
    FITS_open_file(&header, argv[optind], READONLY, &status);
    FITS_read_key(header, TFLOAT, "EXPTIME", &exptime, NULL, &status);

    FITS_movabs_hdu(header, 2, NULL, &status);
    nevents = cf_read_col(header, TFLOAT, "WEIGHT",  (void **) &weight);
    nevents = cf_read_col(header, TFLOAT, "X",       (void **) &x);
    nevents = cf_read_col(header, TFLOAT, "Y",       (void **) &y);
    nevents = cf_read_col(header, TBYTE,  "CHANNEL", (void **) &channel);
    nevents = cf_read_col(header, TBYTE,  "TIMEFLGS", (void **) &timeflags);
    nevents = cf_read_col(header, TBYTE,  "LOC_FLGS", (void **) &loc_flags);
    nevents = cf_read_col(header, TFLOAT, "LAMBDA",  (void **) &lambda);

    FITS_movabs_hdu(header, 4, NULL, &status);
    ntimes = cf_read_col(header, TBYTE, "STATUS_FLAGS", (void **) &time_status);

   /* 
    *  Copy header of IDF into memory, close IDF,
    *  and pass copy of header to subsequent routines.
    */
    FITS_movabs_hdu(header, 1, NULL, &status);
    FITS_create_file(&memp, "mem://", &status);
    FITS_copy_hdu(header, memp, 0, &status);
    FITS_close_file(header, &status);
    header = memp;

    /* If night-only spectrum is requested, modify header accordingly. */
    if (night_only) {
	cf_verbose(3, "Night-only spectrum requested from command line.");
	FITS_read_key(header,  TFLOAT, "EXPNIGHT", &exptime, NULL, &status);
	FITS_update_key(header, TFLOAT, "EXPTIME",  &exptime, NULL, &status);
	FITS_update_key(header, TSTRING, "DAYNIGHT", "NIGHT", NULL, &status);
	FITS_update_key(header, TSTRING, "BPM_CAL", bpm_file, NULL, &status);
    }

    extended = cf_source_aper(header, aper);

    /* Determine whether to attempt optimal extraction. */
    if (override_optimal) {
	cf_verbose (1, "User set -o flag.  No optimal extraction.");
	optimal=FALSE;
    }
    else if (force_optimal) {
	cf_verbose (1, "User set -f flag.  Forcing use of optimal extraction.");
	optimal=TRUE;
    }
    else cf_standard_or_optimal_extraction(header, &optimal);

    /* Make list of all photons that satisfy requested screenings. */
    cf_apply_filters(header, tscreen, nevents, timeflags, loc_flags,
	ntimes, time_status, &dtime, &ntime, &ngood, &good_index);

    /* Generate model backgrounds for LiF and SiC channels. */
    cf_scale_bkgd(header, nevents, x, y, weight, channel, timeflags,
	loc_flags, ngood, good_index, dtime, ntime, &binx, &biny,
	&bnx[0], &bny[0], &bymin[0], &bimage[0],
	&bnx[1], &bny[1], &bymin[1], &bimage[1]);

    /* Run once for each of LiF and SiC target channels. */
    for (j=0; j<2; j++) {
	int   pny, valid_spectrum=TRUE;
	float pycent, wpc;
        long *counts_out=NULL, i, nout=0, nphotons=0;
	float *wave_out=NULL, *bkgd_out=NULL, *flux_out=NULL, *sigma_out=NULL;
	float *weights_out=NULL, *barray=NULL, *parray=NULL;
        short *bpix_out=NULL;
	unsigned char *channel_temp=NULL;

	/* Generate output wavelength array. */
	cf_make_wave_array(header, aper[j], &nout, &wave_out);

	/* If no data, generate null-valued output arrays. */
	for (i = 0; i < ngood; i++) {
	    if (channel[good_index[i]] == aper[j]) {
		nphotons++;
		break;
	    }
	}
	if (exptime < 1. || nphotons < 1) {
            bkgd_out = (float *) cf_calloc(nout, sizeof(float));
            bpix_out = (short *) cf_calloc(nout, sizeof(short));
            counts_out = (long *) cf_calloc(nout, sizeof(long));
            flux_out = (float *) cf_calloc(nout, sizeof(float));
            sigma_out = (float *) cf_calloc(nout, sizeof(float));
            weights_out = (float *) cf_calloc(nout, sizeof(float));
	    valid_spectrum=FALSE;	
	    cf_verbose(1, "No data for aperture %d.  Generating "
		"null-valued output spectrum.", aper[j]);
	}
	else {			/* Extract target spectrum. */

	/* Bin background model to match output wavelength array. */
	    cf_rebin_background(header, aper[j], nout, wave_out,
		binx, biny, bnx[j], bny[j], bimage[j], &barray);

	/* Use QUAL_CAL files to create a bad-pixel mask for this exposure. */
            cf_make_mask(header, aper[j], nout, wave_out, bny[j], bymin[j],
		&bpmask);

	/* Bin the weights/probability array to match output wave array. */
	    cf_rebin_probability_array(header, extended, aper[j], nout,
		wave_out, &pny, &pycent, &parray);

	/* Use optimal extraction if possible, standard extraction otherwise. */
	    cf_optimal_extraction(header, optimal, aper[j],
		weight, y, channel, lambda, ngood, good_index,
		barray, bpmask, pny, pycent, parray, nout, wave_out, &flux_out,
		&sigma_out, &counts_out, &weights_out, &bkgd_out, &bpix_out);

	/* Optimal-extraction routine returns flux_out and sigma_out in units
	   of counts.  Must flux-calibrate both arrays. */
	    channel_temp =
		(unsigned char *) cf_malloc(sizeof(unsigned char) * nout);
	    memset (channel_temp, aper[j], sizeof(unsigned char) * nout);
	    FITS_read_key(header, TFLOAT, "WPC", &wpc, NULL, &status);
	    FITS_update_key(header, TSTRING, "FLUX_COR", "PERFORM", NULL,
		&status) ;
	    cf_convert_to_ergs(header, nout, flux_out, flux_out, channel_temp,
		wave_out);
	    FITS_update_key(header, TSTRING, "FLUX_COR", "PERFORM", NULL,
		&status) ;
	    cf_convert_to_ergs(header, nout, sigma_out, sigma_out, channel_temp,
		wave_out);
	    free(channel_temp);

	/* Divide flux and errors by EXPTIME * WPC to get units of erg/cm2/s/A. */
	    for (i = 0; i < nout; i++) {
		flux_out[i]  /= exptime * wpc;
		sigma_out[i] /= exptime * wpc;
	    }
	}

	/* Write extracted spectrum to output file. */
	cf_write_extracted_spectrum(header, aper[j], valid_spectrum,
		nout, wave_out, flux_out, sigma_out, counts_out,
		weights_out, bkgd_out, bpix_out, outrootname);

	free(counts_out);
	free(bkgd_out);
	free(weights_out);
	free(sigma_out);
	free(flux_out);
	free(wave_out);
        free(bpix_out) ;
	free(parray);
	free(barray);
    }

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    FITS_delete_file(header, &status);

    free(bimage[1]);
    free(bimage[0]);
    free(lambda);
    free(loc_flags);
    free(timeflags);
    free(channel);
    free(y);
    free(x);
    free(weight);
    free(time_status);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
