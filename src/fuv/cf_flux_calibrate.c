/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_flux_calibrate options intermediate_file
 *
 * Description: Flux calibrate each photon assigned to a channel.
 *
 *		** NEXT STEPS COMMENTED OUT FOR NOW **
 *		For each photon in the target aperture,
 *		apply astigmatism correction to XFARF array and use the
 *		resulting coordinate system to apply flat-field and
 *		worm corrections to photon weights.
 *
 *              By default all corrections are performed.  Command line
 *              options allow one or more corrections to be omitted.
 *
 * Arguments:   input_file              Intermediate data file
 *
 * Calibration files:                   AEFF_CAL
 *
 * Returns:     0 on successful completion
 *
 * Calls:       cf_effective_area, cf_astig_farf, cf_flat_field,
 *		cf_worm_correction
 *
 * History:     12/06/02   wvd    1.1   Begin work
 *		02/12/03   wvd    1.2   Read ERGCM2S from IDF
 *					Move cf_convert_to_ergs to end of
 *                                      module
 *		03/11/03   wvd    1.3   Change channel to type char
 *		03/12/03   wvd    1.4   If ERGCM2S[i] is undefined,
 *					set channel[i] = 0
 *              03/19/03   peb    1.5   Added verbose_level option.
 *              05/06/03   rdr    1.6   Read from column ERGCM2 rather
 *                                      than ERGCM2S
 *              05/07/03   wvd    1.7   Change ergcm2s to ergcm2
 *		06/11/03   wvd    1.8   Pass datatype to cf_read_col and
 *					cf_write_col.
 *					Comment out call to cf_astig_farf().
 *		08/25/03   wvd    1.9   Change coltype from string to int in
 *					cf_read_col and cf_write_col.
 *		09/17/03   wvd    1.10  Add some documentation.
 *		09/17/03   wvd    1.11  Change channel to unsigned char.
 *		05/20/05   wvd    1.12  Clean up i/o.
 *
 ****************************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_flux_calibrate";
static char CF_VER_NUM[] = "1.12";

int
main(int argc,  char *argv[])
{
    unsigned char *channel=NULL;
    int  effective_area=1, astig_farf=1, flat_field=1, worm_correction=1;
    int  status=0, hdutype, optc;
    long nevents, nseconds;
    float *time=NULL, *ttime=NULL, *xfarf=NULL, *yfarf=NULL, *weight=NULL;
    float *lambda=NULL, *ycentl=NULL, *ycents=NULL, *ergcm2=NULL;
    fitsfile *header;
    
    char opts[]  = "heafwv:";
    char usage[] = 
	"Usage:\n"
	"  cf_flux_calibrate [-heafw] [-v level] idf_file\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
	"  -v:  verbosity level (default is 1; 0 is silent)\n"
	"  -e:  no effective-area calibration\n"
	"  -a:  no astigmatism correction of XFARF array\n"
	"  -f:  no flat-field correction\n"
	"  -w:  no worm correction\n";

    verbose_level = 1;
 
    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	case 'e':
	    effective_area = 0;
	    break;
	case 'a':
	    astig_farf = 0;
	    break;
	case 'f':
	    flat_field = 0;
	    break;
	case 'w':
	    worm_correction = 0;
	    break;
	}
    }
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    if (argc != optind+1) {
        printf("%s", usage);
        cf_if_error("Incorrect number of command-line arguments");
    }
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");

    FITS_open_file(&header, argv[optind], READWRITE, &status);

    FITS_movabs_hdu(header, 2, &hdutype, &status);
    nevents = cf_read_col(header, TFLOAT, "TIME",   (void **) &time);
    nevents = cf_read_col(header, TFLOAT, "WEIGHT",  (void **) &weight);
    nevents = cf_read_col(header, TFLOAT, "XFARF",  (void **) &xfarf);
    nevents = cf_read_col(header, TFLOAT, "YFARF",  (void **) &yfarf);
    nevents = cf_read_col(header, TBYTE,  "CHANNEL", (void **) &channel);
    nevents = cf_read_col(header, TFLOAT, "LAMBDA",  (void **) &lambda);
    nevents = cf_read_col(header, TFLOAT, "ERGCM2",  (void **) &ergcm2);

    FITS_movabs_hdu(header, 4, &hdutype, &status) ;
    nseconds = cf_read_col(header, TFLOAT, "TIME", (void **) &ttime) ;
    nseconds = cf_read_col(header, TFLOAT, "YCENT_LIF", (void **) &ycentl) ;
    nseconds = cf_read_col(header, TFLOAT, "YCENT_SIC", (void **) &ycents) ;

    FITS_movabs_hdu(header, 1, &hdutype, &status);

   /***************************************************************************
    * Since we currently have no flat-field or worm correction,
    * these subroutines are commented out.
    *
    if (astig_farf)
	cf_astig_farf(header, nevents, xfarf, yfarf, channel, time,
		nseconds, ttime, ycentl, ycents);
    if (flat_field)
	cf_flat_field(header, nevents, time, weight, xfarf, yfarf, channel,
		nseconds, ttime, ycentl, ycents);

    if (worm_correction)
	cf_worm_correction(header, nevents, time, weight, xfarf, yfarf, channel,
		nseconds, ttime, ycentl, ycents);
    ***************************************************************************/

    if (effective_area)
	cf_convert_to_ergs(header, nevents, weight, ergcm2, channel, lambda);

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    FITS_movabs_hdu(header, 2, &hdutype, &status);
    cf_write_col(header, TFLOAT, "WEIGHT",  (void *) weight, nevents);
    cf_write_col(header, TBYTE,  "CHANNEL", (void *) channel, nevents);
    cf_write_col(header, TFLOAT, "ERGCM2",  (void *) ergcm2, nevents);

    FITS_close_file(header, &status);

    free(time);
    free(weight);
    free(xfarf);
    free(yfarf);
    free(channel);
    free(lambda);
    free(ergcm2);
    free(ttime);
    free(ycentl);
    free(ycents);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return 0;
}
