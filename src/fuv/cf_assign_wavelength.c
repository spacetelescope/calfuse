/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_assign_wavelength options intermediate_data_file
 *
 * Description: Assign a wavelength to each photon, correcting for
 *              spacecraft motion and instrumental astigmatism.  
 *		The final wavelength scale is heliocentric.
 *
 *              By default all corrections are performed.  Command-line
 *              options allow one or more corrections to be omitted.
 *
 * Arguments:   input_file              FARF-corrected intermediate data file
 *
 * Calibration files:                   ASTG_CAL, WAVE_CAL
 *
 * Returns:     0 on successful completion
 *
 * Calls:       cf_astigmatism, cf_dispersion, cf_doppler_and_heliocentric
 *
 * History:     11/05/02   1.1   peb    Begin work
 *		12/16/02   1.2   wvd    Install subroutines.
 *		02/12/03   1.3   wvd    Change ORBITAL_VELOCITY to ORBITAL_VEL
 *              03/10/03   1.4   peb    Added verbose_level, changed channel to
 *                                      char *, and added unistd.h for getopt
 *                                      portability
 *		03/12/03   1.5   wvd    If lambda[k] is undefined,
 *					set channel[k] = 0 in IDF
 *		03/23/03   1.6   wvd    If -w option is set, write astigmatism-
 *					corrected X values to IDF.
 *		06/11/03   1.7   wvd    Pass datatype to cf_read_col and
 *					cf_write_col.
 *		08/25/03   1.8   wvd	Change coltype from string to int in
 *					cf_read_col and cf_write_col.
 *		09/16/03   1.9   wvd	Write comment lines to IDF header if
 *					-w option is requested.
 *		10/31/03   1.10  wvd	Change channel to unsigned char.
 *		12/03/03   1.11  wvd	Initialize astig_return to 1.
 *              04/06/04   1.12  bjg    Function returns EXIT_SUCCESS
 *                                      Fix option[]
 *		06/21/04   1.13  wvd	If -w flag is set, run the program
 *					even if it's been run before.
 *		05/20/05   1.14  wvd	Clean up i/o.
 *		05/15/06   1.15  wvd	Divide cf_astigmatism_and_dispersion
 *					into two routines.  Note that wave-
 *					length calibration is performed even
 *					if astigmatism correction is not.
 *
 ****************************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include "calfuse.h"

int
main(int argc,  char *argv[])
{
    char CF_PRGM_ID[] = "cf_assign_wavelength";
    char CF_VER_NUM[] = "1.15";

    char comment[FLEN_CARD], datestr[FLEN_CARD];
    unsigned char *channel=NULL;
    int  astig_corr=1, doppler_corr=1, status=0, optc, timeref;
    int  astig_return=1, wavecal = FALSE;
    long nevents, nseconds;
    float *time=NULL, *x=NULL, *y=NULL, *lambda=NULL, *timeline=NULL;
    float *velocity=NULL;
    fitsfile *header;
    
    char opts[]  = "hadwv:";
    char usage[] = 
	"Usage:\n"
	"  cf_assign_wavelength [-hadw] [-v level] idf_file\n";
    char option[] =
	"Options:\n"
	"  -h:  this help message\n"
	"  -v:  verbose mode (default is 1; 0 is silent)\n"
	"  -a:  no astigmatism correction\n"
	"  -d:  no doppler correction\n"
	"  -w:  write astigmatism-corrected X values to IDF\n";

    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	case 'a':
	    astig_corr = 0;
	    break;
	case 'd':
	    doppler_corr = 0;
	    break;
	case 'w':
	    wavecal = TRUE;
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
    nevents = cf_read_col(header, TFLOAT, "X",      (void **) &x);
    nevents = cf_read_col(header, TFLOAT, "Y",      (void **) &y);
    nevents = cf_read_col(header, TBYTE,  "CHANNEL",(void **) &channel);
    nevents = cf_read_col(header, TFLOAT, "LAMBDA", (void **) &lambda);

    FITS_movabs_hdu(header, 4, NULL, &status);
    nseconds = cf_read_col(header, TFLOAT, "TIME",        (void **) &timeline);
    nseconds = cf_read_col(header, TFLOAT, "ORBITAL_VEL", (void **) &velocity);
    FITS_movabs_hdu(header, 1, NULL, &status);

    if (wavecal) {
        FITS_update_key(header, TSTRING, "WAVE_COR", "PERFORM", NULL, &status);
        FITS_update_key(header, TSTRING, "DOPP_COR", "PERFORM", NULL, &status);
    }

    /* Correct X coordinate for 2-D astigmatism. */
    if (astig_corr)
	astig_return = cf_astigmatism(header, nevents, x, y, channel);

    /* Assign wavelength to each photon event according to channel and X */
    cf_dispersion(header, nevents, x, channel, lambda);

    /* Correct for spacecraft motion and shift to heliocentric wavelength scale. */
    if (doppler_corr)
	cf_doppler_and_heliocentric(header, nevents, time, channel,
	    lambda, nseconds, timeline, velocity);

    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    FITS_movabs_hdu(header, 2, NULL, &status);
    cf_write_col(header, TBYTE,  "CHANNEL", (void *) channel, nevents);
    cf_write_col(header, TFLOAT, "LAMBDA",  (void *) lambda, nevents);

    /* If -w option is requested and astigmatism correction was applied,
	write astigmatism-corrected X array to IDF. */
    if (wavecal && !astig_return) {
	cf_write_col(header, TFLOAT, "X", (void *) x, nevents);
        FITS_movabs_hdu(header, 1, NULL, &status);
	sprintf(comment, "Writing astigmatism-corrected X array to IDF.");
        cf_if_warning(comment);
        FITS_write_comment(header, "  ", &status);
        FITS_write_comment(header, comment, &status);
        fits_get_system_time(datestr, &timeref, &status);
        sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
        FITS_write_comment(header, comment, &status);
        FITS_write_comment(header, "  ", &status);
    }

    FITS_close_file(header, &status);

    free(velocity);
    free(timeline);
    free(lambda);
    free(channel);
    free(y);
    free(x);
    free(time);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done processing");
    return EXIT_SUCCESS;
}
