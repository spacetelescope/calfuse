/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_convert_to_farf options intermediate_file
 *
 * Description: Corrects the RAW data for IDS and electronic deadtime
 *              and for thermal, count rate, geometric, and PHA distortions.
 *              The resulting data are in the Flight Aligned Reference
 *              Frame (FARF), the location of an event on an ideal
 *              detector.  The input file is modified in place.
 *
 *              By default all corrections are performed.  Command line
 *              options allow one or more corrections to be omitted.
 *
 * Arguments:   input_file              Intermediate data file
 *
 * Calibration files:                   
 *
 * Returns:     0 on successful completion
 *
 * Calls:       cf_ids_dead_time, cf_electronics_dead_time,
 *              cf_thermal_distort, cf_count_rate_y_distort,
 *              cf_geometric_distort, cf_pha_x_distort
 *
 * History:     08/08/02   peb    1.1   Begin work
 *              10/27/02   peb    1.3   Added IDS and electronics dead time
 *                                      corrections and time-dependent
 *                                      correction for count_rate_y_distort.
 *                                      Changed the help option to -h.
 *                                      Corrected function description.
 *              11/12/02   peb    1.4   Added argument to cf_thermal_distort,
 *                                      cf_count_rate_y_distort, and
 *                                      cf_geometric_distort calls, so
 *                                      stim-pulse flag can be set.
 *                                      Added function to flag events in
 *                                      inactive regions.
 *		12/09/02   wvd    1.5   Set keyword TOT_DEAD to mean value
 *					for exposure.
 *		01/28/03   wvd    1.6   Added cf_check_digitizer
 *              03/10/03   peb    1.7   Added verbose_level option and unistd.h
 *                                      so getopt.h will be portable.
 *              03/10/03   peb    1.8   Changed locflgs to unsigned char *
 *              03/10/03   peb    1.9   Changed pha to unsigned char *
 *		06/11/03   wvd    1.10  Pass datatype to cf_read_col and
 *					cf_write_col.
 *              06/16/03   rdr    1.11  Correct datatype in reading LOC_FLGS
 *		08/01/03   wvd    1.12  Add cf_apply_dead_time to properly
 *					correct both TTAG and HIST data.
 *		08/04/03   wvd    1.13  Convert count-rate arrays to shorts.
 *		08/25/03   wvd    1.14  Change coltype from string to int in
 *					cf_read_col and cf_write_col.
 *		11/26/03   wvd    1.15  Read aic_rate and fec_rate as floats.
 *		02/09/04   wvd    1.16  Don't pass aic_rate to
 *					cf_apply_dead_time().
 *              02/27/04   rdr    1.17  Added weight to cf_thermal_distortion
 *              03/02/05   wvd    1.18	Must assign channel numbers and pulse-
 *					height values to HIST data before walk
 *					correction.
 *					Read XFARF and YFARF from XRAW and YRAW.
 *		05/20/05   wvd    1.19  Clean up i/o.
 *		11/02/06   wvd    1.20  Add call to cf_time_xy_distort.
 *              12/29/06   wvd    1.21  Move cf_screen_fifo_overflow from
 *                                      cf_screen_photons.  Rename it
 *					cf_fifo_dead_time.  Set all dead-time
 *					keywords to 1.0 by default.
 *              02/08/08   wvd    1.22  For HIST data, write modified PHA
 *					values to the IDF.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "calfuse.h"

int
main(int argc,  char *argv[])
{
    char CF_PRGM_ID[] = "cf_convert_to_farf";
    char CF_VER_NUM[] = "1.22";

    unsigned char *pha=NULL, *locflags=NULL;
    char instmode[FLEN_VALUE], aperture[FLEN_VALUE];
    int  check_digitizer=1, ids_deadtime=1, elec_deadtime=1, fifo_deadtime=1;
    int  thermal_corr=1;
    int  countrate_corr=1, geometric_corr=1, pha_corr=1, active_flag=1;
    int  temporal_corr=1, status=0, hdutype, last_call, optc, pad;
    long i, nevents, nseconds;
    float *elec_dtc=NULL, *ids_dtc=NULL, unity = 1.0;
    float *time=NULL, *xfarf=NULL, *yfarf=NULL, *weight=NULL;
    float *timeline=NULL; 
    float *aic_rate=NULL, *fec_rate=NULL;
    fitsfile *header;
    
    char opts[]  = "hceiftrsgpav:";
    char usage[] = 
	"Usage:\n"
	"  cf_convert_to_farf [-hceiftrsgpa] [-v level] idf_file\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
        "  -v:  verbosity level (default is 1; 0 is silent)\n"
	"  -c:  no check of digitizer values\n"
	"  -e:  no electronics deadtime correction\n"
	"  -i:  no IDS deadtime correction\n"
	"  -f:  no FIFO deadtime correction\n"
	"  -t:  no thermal distortion correction\n"
	"  -r:  no count-rate distortion correction\n"
	"  -s:  no time-dependent distortion correction\n"
	"  -g:  no geometric distortion correction\n"
	"  -p:  no PHA distortion correction\n"
	"  -a:  do not flag events beyond active region\n";

    verbose_level=1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	case 'c':
	    check_digitizer = 0;
	    break;
	case 'i':
	    ids_deadtime = 0;
	    break;
	case 'e':
	    elec_deadtime = 0;
	    break;
	case 'f':
	    fifo_deadtime = 0;
	    break;
	case 't':
	    thermal_corr = 0;
	    break;
	case 'r':
	    countrate_corr = 0;
	    break;
	case 's':
	    temporal_corr = 0;
	    break;
	case 'g':
	    geometric_corr = 0;
	    break;
	case 'p':
	    pha_corr = 0;
	    break;
	case 'a':
	    active_flag = 0;
	    break;
	}
    }
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    if (argc <= optind) {
        printf("%s", usage);
        cf_if_error("Incorrect number of program arguments");
    }
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    FITS_open_file(&header, argv[optind], READWRITE, &status);

    /* Read a couple of keywords from the file header. */
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    FITS_read_key(header, TSTRING, "APERTURE", aperture, NULL, &status);

    FITS_movabs_hdu(header, 2, &hdutype, &status);
    nevents = cf_read_col(header, TFLOAT, "TIME",   (void **) &time);
    nevents = cf_read_col(header, TBYTE,  "PHA",    (void **) &pha);
    nevents = cf_read_col(header, TFLOAT, "XRAW",   (void **) &xfarf);
    nevents = cf_read_col(header, TFLOAT, "YRAW",   (void **) &yfarf);
    nevents = cf_read_col(header, TFLOAT, "WEIGHT", (void **) &weight);
    nevents = cf_read_col(header, TBYTE,  "LOC_FLGS", (void **) &locflags);

    FITS_movabs_hdu(header, 4, &hdutype, &status);
    nseconds = cf_read_col(header, TFLOAT, "TIME",  (void **) &timeline);
    nseconds = cf_read_col(header, TFLOAT, "FEC_CNT_RATE", (void **) &fec_rate);
    nseconds = cf_read_col(header, TFLOAT, "AIC_CNT_RATE", (void **) &aic_rate);
    FITS_movabs_hdu(header, 1, &hdutype, &status);

    /* Initialize dead-time correction arrays and header keywords. */
    FITS_update_key (header, TFLOAT, "DET_DEAD", &unity, NULL, &status);
    FITS_update_key (header, TFLOAT, "IDS_DEAD", &unity, NULL, &status);
    FITS_update_key (header, TFLOAT, "TOT_DEAD", &unity, NULL, &status);
    elec_dtc = (float *) cf_malloc(sizeof(float) * nseconds);
    ids_dtc  = (float *) cf_malloc(sizeof(float) * nseconds);
    for (i = 0; i < nseconds; i++) elec_dtc[i] = ids_dtc[i] = 1.0;

    if (check_digitizer)
	cf_check_digitizer(header);

    if (elec_deadtime)
	cf_electronics_dead_time(header, nseconds, fec_rate, elec_dtc);

    if (ids_deadtime)
	cf_ids_dead_time(header, nseconds, aic_rate, ids_dtc);

    if (fifo_deadtime)
        cf_fifo_dead_time(header, nevents, time,
                nseconds, timeline, aic_rate, ids_dtc);

    if (ids_deadtime || elec_deadtime || fifo_deadtime)
	cf_apply_dead_time(header, nevents, time, weight,
			 nseconds, timeline, ids_dtc, elec_dtc);

    if (thermal_corr)
	cf_thermal_distort(header, nevents, xfarf, yfarf, weight, locflags);

    if (countrate_corr)
	cf_count_rate_y_distort(header, nevents, time, yfarf, locflags,
				nseconds, timeline, fec_rate);
    if (temporal_corr)
	cf_time_xy_distort(header, nevents, xfarf, yfarf, locflags);

    if (geometric_corr)
	cf_geometric_distort(header, nevents, xfarf, yfarf, locflags);

    /* Must assign PHA values to HIST data before applying walk correction. */
    if (pha_corr && !strncmp(instmode, "HIST", 4)) {
	unsigned char *channel;
	channel = cf_calloc(nevents, sizeof(unsigned char));
	if (!strncmp(aperture, "RFPT", 4)) pad = 10;
	else pad = 100;
	cf_identify_channel(header, nevents, xfarf, yfarf, channel, locflags,
	    pad, last_call=FALSE);
	cf_modify_hist_pha(header, nevents, pha, channel);
	free(channel);
    }

    if (pha_corr)
	cf_pha_x_distort(header, nevents, pha, xfarf, locflags);

    if (active_flag)
	cf_active_region(header, nevents, xfarf, yfarf, locflags);

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    FITS_movabs_hdu(header, 2, &hdutype, &status);
    if (pha_corr && !strncmp(instmode, "HIST", 4))
	cf_write_col(header, TBYTE,  "PHA", (void *) pha, nevents);
    cf_write_col(header, TFLOAT, "WEIGHT", (void *) weight, nevents);
    cf_write_col(header, TFLOAT, "XFARF",  (void *) xfarf,  nevents);
    cf_write_col(header, TFLOAT, "YFARF",  (void *) yfarf,  nevents);
    cf_write_col(header, TBYTE,  "LOC_FLGS", (void *) locflags, nevents);

    FITS_close_file(header, &status);

    free(fec_rate);
    free(aic_rate);
    free(timeline);
    free(locflags);
    free(weight);
    free(yfarf);
    free(xfarf);
    free(pha);
    free(time);
    free(elec_dtc);
    free(ids_dtc);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
