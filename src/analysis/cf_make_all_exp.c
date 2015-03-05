/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis: cf_make_all_exp output_filename input1 input2 input3 input4
 *
 * Description: For a single exposure, copies header of channel with largest	
 *              value of EXPTIME into an otherwise empty "all" file.
 *		Input files should represent each of the four FUSE channels.
 *
 *
 * History:	08/11/05	wvd	v1.0	Based on cf_quicklook.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "calfuse.h"

static char CF_PRGM_ID[]= "cf_make_all_exp";
static char CF_VER_NUM[]= "1.0";

int main(int argc,char *argv[]){
  
	char	*corrected_filename, *string_pointer;
	char	comment[FLEN_CARD], datestr[FLEN_CARD];
	int	i, nextend=0, ref_file=2, timeref;
	int	status=0;
	float	exptime, max_exptime;

	fitsfile *infits,*outfits; 

	if (argc != 6) {
		printf("Incorrect number of arguments.\n");
		printf("Calling sequence:  cf_make_all_exp output_filename "
			"input1 input2 input3 input4\n");
		exit(1);
	}

	/* Initialize error checking. */
	cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");

	/* Determine the channel with the longest exposure time */

	max_exptime = 0;
	for (i = 0; i < 4; i++) {
		FITS_open_file(&infits, argv[2+i], READONLY, &status);
		FITS_read_key(infits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
		if (exptime >= max_exptime) {
        		ref_file = 2 + i;
    			max_exptime = exptime;
		}
		FITS_close_file(infits, &status);
	}

	/* Copy header from the reference file to the exposure-level "all" file. */

	FITS_open_file(&infits, argv[ref_file], READONLY, &status);
	FITS_create_file(&outfits, argv[1], &status);
	FITS_copy_hdu(infits, outfits, 0, &status);

	/* Modify a few header keywords. */
	string_pointer = strrchr(argv[1],'/');
	if (string_pointer == NULL) corrected_filename = argv[1];
	else corrected_filename = &(string_pointer[1]);
	FITS_update_key(outfits, TINT, "NEXTEND", &nextend, NULL, &status);
	FITS_write_date(outfits, &status);
	FITS_update_key(outfits, TSTRING, "FILENAME", corrected_filename, NULL, &status);
	FITS_update_key(outfits, TSTRING, "FILETYPE", "EXPOSURE-LEVEL STATUS FILE", NULL, &status);

	/* Add a note warning that this file contains no data. */
        FITS_write_comment(outfits, "  ", &status);
        FITS_write_comment(outfits,
                "This file contains no data.", &status);
        FITS_write_comment(outfits,
                "Its header keywords are used by the DADS cataloging software.", &status);
        fits_get_system_time(datestr, &timeref, &status);
        sprintf(comment, "CalFUSE v%s   %.10s", CALFUSE_VERSION, datestr);
        FITS_write_comment(outfits, comment, &status);
        FITS_write_comment(outfits, "  ", &status);


	FITS_close_file(infits, &status);
	FITS_close_file(outfits,&status);

	return EXIT_SUCCESS;
}
