/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences 
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_remove_motions options intermediate_file     
 *
 * Description: Associates a photon with an aperture and corrects for grating,
 *              mirror, spacecraft (jitter) and focal plane assembly (FPA)
 *              motion.  The position of the Y centroid is also determined
 *		both as a function of time (before correction for s/c motion)
 *		and for the entire exposure (after motion correction).
 *              The resulting coordinates represent the data obtained by a
 *		motionless instrument. 
 *
 *              By default, all corrections are performed.  Command line
 *              options allow one or more corrections to be omitted.
 *
 * Arguments:   input_file              Intermediate data file
 *
 * Calibration files:                   
 *
 * Returns:     0 on successful completion
 *
 * Calls:       cf_find_spectra, cf_calculate_ycent_motion,
 *              cf_grating_motion, cf_fpa_position, cf_mirror_motion,
 *              cf_satellite_jitter, cf_calculate_y_centroid
 *
 * History:     09/30/02   peb    1.1   Begin work
 *		12/04/02   wvd    1.5   Install Y centroid routines.
 *                                      Modify call to cf_satellite_jitter
 *		12/12/02   wvd    1.6   Change infits -> header throughout
 *              03/04/03   peb    1.10  Minor code changes to remove gcc -Wall
 *                                      messages: unused j in main.
 *                                      Added unistd.h, so getopt is portable
 *		03/05/03   wvd    1.11  Added cf_target_count_rate
 *              03/10/03   peb    1.12  Added verbose_level, unistd.h for
 *                                      getopt portability, and changed
 *                                      locflags to unsigned char *.
 *              03/13/03   rdr    1.13  corrected small error in tabulating
 *                                      count rate_sic
 *              04/17/03   wvd	  1.17  Add cf_find_spectra for initial
 *					assignment of photons to
 *					target apertures.  Add last_call to
 *					cf_identify_channel, last_call and
 *					weight to cf_calculate_y_centroid
 *		04/21/03   wvd    1.18  Pass tsunrise to cf_grating_motion,
 *					tsunset to cf_mirror_motion,
 *					channel to cf_satellite_jitter.
 *		06/02/03   wvd    1.20  Don't modify count-rate arrays
 *					for HIST data.
 *		06/11/03   wvd    1.21  Pass datatype to cf_read_col and
 *					cf_write_col.
 *		08/01/03   wvd    1.23  Add cf_error_init after calls to
 *					subroutines.
 *		08/06/03   wvd    1.24  Change channel to unsigned char.
 *					Remove GTI's from cf_satellite_jitter.
 *		08/25/03   wvd    1.25  Change coltype from string to int in
 *					cf_read_col and cf_write_col.
 *		09/18/03   wvd    1.26  Pass locflags to cf_identify_channel
 *		10/21/03   wvd    1.27  For HIST data, don't write blank
 *					LiF and SiC count-rate arrays to IDF.
 *		10/27/03   wvd    1.28  Pass pha array to cf_find_spectra.
 *					Delete last_call from argument list
 *					of cf_calculate_y_centroid.
 *		10/31/03   wvd    1.29  Treat channel as unsigned char
 *					throughout.
 *              04/06/04   bjg    1.30  Remove unused variables
 *		05/03/04   wvd    1.31  Improve documentation.
 *					Move first call of cf_identify_channel
 *					from cf_find_spectra to this routine.
 *		06/02/04   wvd    1.32  Copy/move cf_set_photon_flags(),
 *					cf_set_good_time_intervals(), and
 *					cf_modify_hist_times() from 
 *					cf_screen_photons.  Run them only if
 *					cf_satellite_jitter rejects some time
 *					AND if they were already run before.
 *					Pass timeflags to cf_find_spectra()
 *					and cf_calculate_y_centroid().
 *					Get exp_jitr from cf_satellite_jitter.
 *		07/21/04   wvd    1.33  Add include file <string.h>
 *		03/10/05   wvd    1.34  Change use of variable pad.  Now add
 *					10 to user_pad on first call to 
 *					cf_identify_channel, 0 on second call.
 *		03/22/05   wvd    1.35  Change TIME_SUNRISE and TIME_SUNSET
 *					from floats to shorts.
 *		05/20/05   wvd    1.36  Clean up i/o.
 *		09/09/05   wvd    1.37  Add -j to list of allowed options.
 *		11/29/05   wvd    1.38  Separate cf_satellite_jitter into two
 *					routines.  Screening is performed 
 *					earlier, so we need not copy flags.
 *		08/30/06   wvd    1.39  Add -a to list of allowed options.
 *		12/29/06   wvd    1.40  Pass weights array and initial LiF
 *					and SiC counts arrays to 
 *					cf_target_count_rate.  Make arrays
 *					rate_lif and rate_sic floats.  If 
 *					they are too large to store as shorts, 
 *					store as unsigned shorts.
 *		02/02/07   wvd    1.41  If LiF or SiC count rate is too large
 *					for an unsigned int, set to zero.
 *
 ****************************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_remove_motions";
static char CF_VER_NUM[] = "1.41";

int
main(int argc,  char *argv[])
{
    unsigned char *channel=NULL, *timeflags=NULL, *locflags=NULL,
	*statflag=NULL;
    int  grating_motion=1, fpa_position=1, mirror_motion=1, sat_jitter=1;
    int  last_call, set_times=FALSE, set_gtis=FALSE, airglow_centroid=FALSE;
    int  user_pad=0, status=0, set_rate=FALSE, optc;
    long i, nevents, ntimes;
    float *time=NULL, *ttime=NULL, *xfarf=NULL, *yfarf=NULL, *x, *y;
    float *weight=NULL, *ycent_lif, *ycent_sic;
    float max_rate_lif=0, max_rate_sic=0, *rate_lif=NULL, *rate_sic=NULL;
    float tscale=1, tzero=32768; 
    short *tsunrise=NULL, *tsunset=NULL;
    GTI  gti;
    fitsfile *header;
    
    char opts[]  = "hagfmjp:v:";
    char usage[] = 
	"Usage:\n"
	"  cf_remove_motions [-hagfmj] [-p pad] [-v level] idf_file\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
	"  -v:  verbose level (default is 1; 0 is silent)\n"
	"  -a:  use airglow lines to determine spectral centroid\n"
	"  -p:  aperture padding in pixels (default is 10)\n"
	"  -g:  no grating motion correction\n"
	"  -f:  no fpa position correction\n"
	"  -m:  no mirror motion correction\n"
	"  -j:  no satellite jitter correction\n";

    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	case 'p':
	    user_pad = atoi(optarg);
	    break;
	case 'g':
	    grating_motion = 0;
	    break;
	case 'f':
	    fpa_position = 0;
	    break;
	case 'a':
	    airglow_centroid = TRUE;
	    break;
	case 'm':
	    mirror_motion = 0;
	    break;
	case 'j':
	    sat_jitter = 0;
	    break;
	}
    }
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    if (argc <= optind) {
        printf("%s", usage);
        cf_if_error("Incorrect number of command-line arguments");
    }
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    FITS_open_file(&header, argv[optind], READWRITE, &status);
 
    FITS_movabs_hdu(header, 2, NULL, &status);
    nevents = cf_read_col(header, TFLOAT, "TIME",   (void **) &time);
    nevents = cf_read_col(header, TFLOAT, "WEIGHT", (void **) &weight);
    nevents = cf_read_col(header, TFLOAT, "XFARF",  (void **) &xfarf);
    nevents = cf_read_col(header, TFLOAT, "YFARF",  (void **) &yfarf);
    nevents = cf_read_col(header, TBYTE,  "TIMEFLGS", (void **) &timeflags);
    nevents = cf_read_col(header, TBYTE,  "LOC_FLGS", (void **) &locflags);
    channel = cf_calloc(nevents, sizeof(unsigned char));

    FITS_movabs_hdu(header, 4, NULL, &status) ;
    ntimes = cf_read_col(header, TFLOAT, "TIME", (void **) &ttime) ;
    ntimes = cf_read_col(header, TBYTE,  "STATUS_FLAGS", (void **) &statflag);
    ntimes = cf_read_col(header, TSHORT, "TIME_SUNRISE", (void **) &tsunrise) ;
    ntimes = cf_read_col(header, TSHORT, "TIME_SUNSET",  (void **) &tsunset) ;
    ntimes = cf_read_col(header, TFLOAT, "LIF_CNT_RATE",  (void **) &rate_lif) ;
    ntimes = cf_read_col(header, TFLOAT, "SIC_CNT_RATE",  (void **) &rate_sic) ;
    ycent_lif=(float *) cf_calloc(ntimes, sizeof(float)) ;
    ycent_sic=(float *) cf_calloc(ntimes, sizeof(float)) ;

    FITS_movabs_hdu(header, 1, NULL, &status);
    x = xfarf;
    y = yfarf;

    /* Find all six spectra and compute centroids. */
    cf_find_spectra(header, nevents, weight, x, y, channel, timeflags,
        locflags, airglow_centroid);

    /* Initial assignment of channel numbers to photons */
    cf_identify_channel(header, nevents, xfarf, yfarf, channel, locflags,
	user_pad+10, last_call=FALSE);

    /* Track Y centroids of moving spectra. */
    cf_calculate_ycent_motion(header, nevents, time, y, channel, locflags,
	ntimes, ttime, ycent_lif, ycent_sic);

    /* Correct for detector motion during exposure. */
    if (grating_motion)
	cf_grating_motion(header, nevents, time, x, y, channel,
		ntimes, ttime, tsunrise);

    /* Shift spectra in X to account for position of the FPA. */
    if (fpa_position)
	cf_fpa_position(header, nevents, x, channel);

    /* Correct for mirror motion during exposure. */
    if (mirror_motion)
	cf_mirror_motion(header, nevents, time, x, y, channel,
		ntimes, ttime, tsunset);

    /* Correct for spacecraft motion during exposure. */
    if (sat_jitter)
	cf_satellite_jitter(header, nevents, time, x, y, channel,
		ntimes, ttime, statflag);
 
    /* Compute centroids of motion-corrected spectra. */
    cf_calculate_y_centroid(header, nevents, weight, x, y, channel, 
	timeflags, locflags);

    /* Final assignment of channel numbers to photons */
    cf_identify_channel(header, nevents, x, y, channel, locflags,
	user_pad, last_call=TRUE); 

    /* Compute count rates for TTAG target spectra. */
    set_rate = !(cf_target_count_rate(header, nevents, time, weight,
	channel, locflags, ntimes, ttime, rate_lif, rate_sic)); 

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    cf_verbose(3, "Writing X, Y, CHANNEL to IDF");
    FITS_movabs_hdu(header, 2, NULL, &status);
    if (set_times)
        cf_write_col(header, TFLOAT, "TIME", (void *) time, nevents);
    cf_write_col(header, TFLOAT, "X",      (void *) x,   nevents);
    cf_write_col(header, TFLOAT, "Y",      (void *) y,   nevents);
    cf_write_col(header, TBYTE, "CHANNEL", (void *) channel, nevents);

    if (set_gtis) {
	cf_verbose(3, "Writing good-time intervals to IDF");
	FITS_movabs_hdu(header, 3, NULL, &status);
	cf_write_col(header, TDOUBLE, "START", (void *) gti.start, gti.ntimes);
	cf_write_col(header, TDOUBLE, "STOP", (void *) gti.stop, gti.ntimes);
	free(gti.stop);
	free(gti.start);
    }

    FITS_movabs_hdu(header, 4, NULL, &status);
    cf_verbose(3, "Writing timeline information to IDF");
    if (set_rate) {
	for (i = 0; i < ntimes; i++) {
	    if (rate_lif[i] > 65535) rate_lif[i] = 0;
	    if (rate_sic[i] > 65535) rate_sic[i] = 0;
	    if (max_rate_lif < rate_lif[i]) max_rate_lif = rate_lif[i];
	    if (max_rate_sic < rate_sic[i]) max_rate_sic = rate_sic[i];
	}
	cf_verbose(2, "max_rate_lif = %.0f", max_rate_lif);
	cf_verbose(2, "max_rate_sic = %.0f", max_rate_sic);
	if (max_rate_lif > 32767) {
	    FITS_update_key(header, TFLOAT, "TSCAL10", &tscale, NULL, &status);
	    FITS_update_key(header, TFLOAT, "TZERO10", &tzero, NULL, &status);
	    fits_set_tscale(header, 10, (double)tscale, (double)tzero, &status);
	}
	if (max_rate_sic > 32767) {
	    FITS_update_key(header, TFLOAT, "TSCAL11", &tscale, NULL, &status);
	    FITS_update_key(header, TFLOAT, "TZERO11", &tzero, NULL, &status);
	    fits_set_tscale(header, 11, (double)tscale, (double)tzero, &status);
	}
	FITS_write_col(header, TFLOAT, 10, 1, 1, ntimes, rate_lif, &status);
	FITS_write_col(header, TFLOAT, 11, 1, 1, ntimes, rate_sic, &status);
    }
    cf_write_col(header, TFLOAT, "YCENT_LIF", (void *) ycent_lif, ntimes);
    cf_write_col(header, TFLOAT, "YCENT_SIC", (void *) ycent_sic, ntimes);
    cf_verbose(3, "IDF updates complete");

    FITS_close_file(header, &status);

    free(time);
    free(weight);
    free(xfarf);
    free(yfarf);
    free(locflags);
    free(timeflags);
    free(channel);
    free(ttime);
    free(rate_lif);
    free(rate_sic);
    free(ycent_lif);
    free(ycent_sic);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return status;
}
