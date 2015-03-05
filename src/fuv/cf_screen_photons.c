/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_screen_photons options intermediate_file
 *
 * Description: Filters events in the intermediate data file (IDF) for limb
 *              angle, SAA crossing, high voltage changes, bursts, and PHA.
 *              Sets the event flags and GTIs.  Data are modified in place.
 *
 *              By default all corrections are performed.  Command line
 *              options allow one or more corrections to be omitted.
 *
 * Arguments:   input_file              FARF-corrected intermediate data file
 *
 * Calibration files:                   
 *
 * Returns:     0 on successful completion
 *
 * Calls:       cf_screen_limb_angle, cf_screen_saa, cf_screen_high_voltage,
 *              cf_screen_burst, cf_screen_jitter,
 *		cf_set_user_gtis, cf_set_photon_flags, cf_modify_hist_times, 
 *		cf_screen_pulse_height, cf_screen_airglow, cf_screen_bad_pixels
 *
 * History:     11/05/02   1.1   peb    Begin work
 *              11/13/02   1.2   peb    Added screening functions to program
 *              11/14/02   1.3   peb    Added burst-screening function
 *              11/18/02   1.4   peb    Corrected logic dealing with reading,
 *                                      writing, and freeing of GTI struct.
 *                                      Corrected some grammatical errors.
 *              12/10/02   1.5   rdr    changed call to cf_screen_bursts
 *              12/18/02   1.6   rdr    updated names of some columns in
 *                                      the calls to the timeline extension
 *              12/20/02   1.7   wvd    Change flags to unsigned char
 *              03/10/03   1.9   peb    Added verbose_level and added unistd.h
 *                                      for getopt portability
 *              05/10/03   1.10  wvd    Pass locflag to cf_set_photon_flags
 *              05/30/03   1.11  wvd    Pass weight to cf_set_photon_flags
 *              06/11/03   1.12  wvd    Treat HV as array of shorts.
 *		                        Pass datatype to cf_read_col and
 *					cf_write_col.
 *		08/25/03   1.13  wvd	Change coltype from string to int in
 *					cf_read_col and cf_write_col.
 *		09/10/03   1.14  wvd	Change background array to type short.
 *					Add cf_set_user_gtis().
 *		10/02/03   1.15  wvd	Pass locflag array to
 *					cf_screen_pulse_height; write to IDF
 *		02/10/04   1.16  wvd	Add cf_screen_fifo_overflow()
 *		06/02/04   1.17  wvd	Add cf_modify_hist_times()
 *		06/03/04   1.18  wvd	Read XFARF and YFARF, not X and Y.
 *					If INSTMODE = TTAG, don't bother to
 *					call cf_modify_hist_times().
 *		02/02/05   1.19  wvd	Pass rate_aic to cf_screen_burst
 *		03/10/05   1.20  wvd	Add cf_screen_airglow()
 *		03/30/05   1.21  wvd	Don't change TIME_COR to SKIPPED
 *					for TTAG data.  Leave as OMIT.
 *		05/20/05   1.22  wvd	Clean up i/o.
 *		06/03/05   1.23  wvd	Fix typo in CF_VER_NUM.
 *		09/09/05   1.24  wvd	Update list of allowed options.
 *		11/09/05   1.25  wvd	Change argument list for
 *					cf_screen_fifo_overflow.
 *		11/22/05   1.26  wvd	Add cf_screen_jitter and
 *					cf_screen_bad_pixels.
 *		01/24/06   1.27  wvd	Move cf_screen_airglow before
 *					cf_screen_burst
 *              12/29/06   1.28  wvd    Move cf_screen_fifo_overflow to
 *                                      cf_convert_to_farf.
 *              07/18/08   1.29  wvd    Write STATUS_FLAGS to IDF before
 *					calling cf_set_photon_flags, which
 *					may change the array.
 *
 ****************************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

int
main(int argc, char *argv[])
{
    char CF_PRGM_ID[] = "cf_screen_photons";
    char CF_VER_NUM[] = "1.29";

    char instmode[FLEN_VALUE];
    unsigned char *pha=NULL, *timeflag=NULL, *statflag=NULL, *locflag=NULL;
    int  limb_screening=1, saa_screening=1, voltage_screening=1;
    int  burst_screening=1, set_flags=1, airglow_screening=1, pha_screening=1;
    int  modify_times=1, set_gtis=1, set_user_gtis=1;
    int  bad_pixel_screening=1, set_times=FALSE, jitter_screening=1;
    int  status=0, optc;
    long nevents, nseconds;
    short *background=NULL, *rate_lif=NULL, *rate_sic=NULL, *voltage=NULL;
    float *time=NULL, *x=NULL, *y=NULL, *weight=NULL;
    float *timeline=NULL, *limb=NULL, *longitude=NULL, *latitude=NULL;
    float *rate_aic=NULL;
    GTI  gti;
    fitsfile *header;
    
    char opts[]  = "halstbjufimdpv:";
    char usage[] = 
	"Usage:\n"
	"  cf_screen_photons [-halstbjufimdp] [-v level] idf_file\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
	"  -v:  verbosity level (default is 1; 0 is silent)\n"
	"  -a:  no airglow screening\n"
	"  -l:  no limb screening\n"
	"  -s:  no SAA screening\n"
	"  -t:  no high voltage screening\n"
	"  -b:  no burst screening\n"
	"  -j:  no jitter screening\n"
	"  -u:  no user-defined good-time intervals\n"
	"  -f:  do not set photon event flags\n"
        "  -i:  do not update GTI arrays\n"
	"  -m:  do not modify HIST photon times\n"
	"  -d:  no bad-pixel screening\n"
	"  -p:  no PHA screening\n";

    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	case 'l':
	    limb_screening = 0;
	    break;
	case 't':
	    voltage_screening = 0;
	    break;
	case 's':
	    saa_screening = 0;
	    break;
	case 'b':
	    burst_screening = 0;
	    break;
	case 'j':
	    jitter_screening = 0;
	    break;
	case 'u':
	    set_user_gtis = 0;
	    break;
	case 'f':
	    set_flags = 0;
	    break;
        case 'i':
            set_gtis = 0;
            break;
	case 'm':
	    modify_times = 0;
	    break;
	case 'a':
	    airglow_screening = 0;
	    break;
	case 'd':
	    bad_pixel_screening = 0;
	    break;
	case 'p':
	    pha_screening = 0;
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

    FITS_movabs_hdu(header, 2, NULL, &status);
    nevents = cf_read_col(header, TFLOAT, "TIME",   (void **) &time);
    nevents = cf_read_col(header, TBYTE,  "PHA",    (void **) &pha);
    nevents = cf_read_col(header, TFLOAT, "WEIGHT", (void **) &weight);
    nevents = cf_read_col(header, TFLOAT, "XFARF",  (void **) &x);
    nevents = cf_read_col(header, TFLOAT, "YFARF",  (void **) &y);
    nevents = cf_read_col(header, TBYTE,  "TIMEFLGS", (void **) &timeflag);
    nevents = cf_read_col(header, TBYTE,  "LOC_FLGS", (void **) &locflag);

    if (burst_screening) {
	FITS_movabs_hdu(header, 3, NULL, &status);
	gti.ntimes = cf_read_col(header, TDOUBLE, "START",(void **) &gti.start);
	gti.ntimes = cf_read_col(header, TDOUBLE, "STOP", (void **) &gti.stop);
    }
    FITS_movabs_hdu(header, 4, NULL, &status);
    nseconds = cf_read_col(header, TFLOAT, "TIME",         (void **) &timeline);
    nseconds = cf_read_col(header, TBYTE,  "STATUS_FLAGS", (void **) &statflag);
    nseconds = cf_read_col(header, TFLOAT, "LIMB_ANGLE",   (void **) &limb);
    nseconds = cf_read_col(header, TFLOAT, "LONGITUDE",   (void **) &longitude);
    nseconds = cf_read_col(header, TFLOAT, "LATITUDE",     (void **) &latitude);
    nseconds = cf_read_col(header, TSHORT, "HIGH_VOLTAGE", (void **) &voltage);
    nseconds = cf_read_col(header, TSHORT, "LIF_CNT_RATE", (void **) &rate_lif);
    nseconds = cf_read_col(header, TSHORT, "SIC_CNT_RATE", (void **) &rate_sic);
    nseconds = cf_read_col(header, TFLOAT, "AIC_CNT_RATE", (void **) &rate_aic);
    background = (short *) cf_calloc(nseconds, sizeof(short)) ;

    /* If INSTMODE = TTAG, don't bother to call cf_modify_hist_times. */
    FITS_movabs_hdu(header, 1, NULL, &status);
    FITS_read_key(header, TSTRING, "INSTMODE", instmode, NULL, &status);
    if (!strncmp(instmode, "TTAG", 4)) modify_times = FALSE;

    if (airglow_screening)
	cf_screen_airglow(header, nevents, x, y, locflag);

    if (limb_screening)
	cf_screen_limb_angle(header, nseconds, statflag, limb);

    if (saa_screening)
	cf_screen_saa(header, nseconds, statflag, longitude, latitude);

    if (voltage_screening)
	cf_screen_high_voltage(header, nseconds, statflag, voltage);

    if (burst_screening) {
	cf_screen_burst(header, nevents, time, x, y, locflag, &gti,
			nseconds, timeline, statflag, rate_aic, background);
	free(gti.stop);
	free(gti.start);
    }

    if (jitter_screening)
	cf_screen_jitter(header, nseconds, timeline, statflag);

    if (set_user_gtis)
	cf_set_user_gtis(header, nseconds, timeline, statflag);

    /* Write STATUS_FLAGS array to IDF before calling cf_set_photon_flags. */
    FITS_movabs_hdu(header, 4, NULL, &status);
    cf_write_col(header, TBYTE, "STATUS_FLAGS", (void *) statflag, nseconds);
    FITS_movabs_hdu(header, 1, NULL, &status);

    if (set_flags)
	cf_set_photon_flags(header, nevents, time, weight, timeflag,
		locflag, nseconds, timeline, statflag);

    if (set_gtis)
        cf_set_good_time_intervals(header, nseconds, timeline,
                statflag, &gti);

    if (modify_times)
	set_times = !(cf_modify_hist_times(header, nevents, time, &gti));

    if (bad_pixel_screening)
	cf_screen_bad_pixels(header, nevents, x, y, locflag);

    if (pha_screening)
	cf_screen_pulse_height(header, nevents, pha, locflag);

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    FITS_movabs_hdu(header, 2, NULL, &status);
    if (set_times)
	cf_write_col(header, TFLOAT, "TIME", (void *) time, nevents);
    cf_write_col(header, TBYTE, "TIMEFLGS", (void *) timeflag, nevents);
    cf_write_col(header, TBYTE, "LOC_FLGS", (void *) locflag, nevents);

    if (set_gtis) {
        FITS_movabs_hdu(header, 3, NULL, &status);
        cf_write_col(header, TDOUBLE, "START", (void *) gti.start, gti.ntimes);
        cf_write_col(header, TDOUBLE, "STOP", (void *) gti.stop, gti.ntimes);
        free(gti.stop);
        free(gti.start);
    }

    if (burst_screening) {
	FITS_movabs_hdu(header, 4, NULL, &status);
        cf_write_col(header, TSHORT, "BKGD_CNT_RATE", (void *) background,
	nseconds);
    }

    FITS_close_file(header, &status);

    free(background);
    free(voltage);
    free(latitude);
    free(longitude);
    free(limb);
    free(statflag);
    free(timeline);
    free(timeflag);
    free(pha);
    free(y);
    free(x);
    free(time);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
