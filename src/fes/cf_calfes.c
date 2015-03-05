/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:    cf_calfes rootname
 *                   where rootname is the exposure root name
 *
 * Description: 
 *      This main program calls all the routines below to apply 
 *      FES processing steps.  It
 *         1) Opens the input and output FES FITS files;
 *         2) Loops through all images present:
 *              a) Calls the procedures listed below to process all the images;
 *              b) Writes the current processed image out to the 
 *                    current output HDU.
 *         3) Closes the input and output files.
 *
 * Arguments:  rootname exposure_root_name
 *
 * Returns:    0 upon successful completion
 *
 * History:    07/08/98    gak    calfes_design.070898 design documented
 *             07/15/98    mlr    begin work
 *             04/12/99    mlr    modified to only write each image only
 *                                once - now passing the image in addition
 *                                to the fits hdu to each subroutine.
 *             04/19/99    mlr    finished modifications to utilize
 *                                libcf and FITSIO.h(error handling).
 *	       08/23/04	   wvd    Change cf_velang_calc to cf_velang.
 *
 ******************************************************************************/

#include <string.h>
#include "calfuse.h"
#include "cf_calfes.h"

#define CF_PRGM_ID      "cf_calfes"
#define CF_VER_NUM	"1.4"
#define RAW_FES "RAW FES           "
#define CAL_FES "CALIBRATED FES    "

int main(int argc, char *argv[])
{
    fitsfile   *infits, *outfits;
    int        i, status;
    int        hdu, hdutype, num_of_hdus;
    int        naxes1, naxes2;
    long       bitpix;
    float      bscale, bzero;
    double     mjd_start, mjd_end;

    char       infilename[30], outfilename[30];
    char       file_type[20];
    float      *image;

        /* Enter a timestamp into the log. */
    	cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin Processing");
        cf_error_init(CF_PRGM_ID,CF_VER_NUM, stdout);

        if (argc != 2)
           cf_if_error("Usage: cf_calfes filerootname ");

        /* determine file names of the input file and the output file */
	strcpy(infilename, argv[1]);
        strcat(infilename, "raw.fit");

        strcpy(outfilename, argv[1]);
        strcat(outfilename, "cal.fit");

        status=0;

        /* open the input FES FITS file  and verify that it is a RAW FES file*/
	FITS_open_file(&infits, infilename, READONLY, &status);

        FITS_read_key(infits, TSTRING, "FILETYPE", &file_type, NULL, &status);

        if (strncmp(file_type, RAW_FES, 7) != 0)
                cf_if_error("The input file is not a RAW FES file.");

        /* open the output file for calibrated FES data, copy the primary
           header from the input FES file to the output FES file, then 
           modify keywords as needed. 
        */
	FITS_create_file(&outfits, outfilename, &status);
        FITS_copy_header(infits, outfits, &status);

        bitpix = FLOAT_IMG;
        FITS_update_key(outfits, TLONG,   "BITPIX",   &bitpix, NULL, &status);
        FITS_update_key(outfits, TSTRING, "FILETYPE", CAL_FES, NULL, &status);

        read_tle(outfits);

        /* Determine how many extensions there are in this file.
           There is 1 FES image in each extension.  So this tells us how
           many images there are in this file and controls the loop that
           processes each one.
        */

        FITS_read_key(infits, TINT, "NEXTEND ", &num_of_hdus, NULL, &status);

        printf("This file %s has %d FES image(s).\n", infilename, num_of_hdus);
        hdutype = -1;

        for (i=1 ; i <= num_of_hdus; i++)
        {
            hdu = i + 1;

	    FITS_create_hdu(outfits, &status);
	    FITS_movabs_hdu(infits, hdu, &hdutype, &status);
            FITS_copy_header(infits, outfits, &status);

            bitpix = FLOAT_IMG;
            bzero = 0.0; 
            bscale = 1.0;
            FITS_update_key(outfits, TLONG,  "BITPIX", &bitpix, NULL, &status);
            FITS_update_key(outfits, TFLOAT, "BZERO",  &bzero,  NULL, &status);
            FITS_update_key(outfits, TFLOAT, "BSCALE", &bscale, NULL, &status);

            FITS_read_key(outfits,TDOUBLE,"EXPSTART",&mjd_start,NULL, &status);
            FITS_read_key(outfits,TDOUBLE,"EXPEND",  &mjd_end,  NULL, &status);
            cf_velang(outfits, (mjd_start+mjd_end)/2.0 );
            cf_min_limbang(outfits, mjd_start, mjd_end);
 
             status =+ cf_fes_read(infits, &image, &naxes1, &naxes2);
             status =+ cf_fes_init(outfits);
             status =+ cf_fes_write(outfits, hdu, image); 

             status =+ cf_fes_cal(FES_MASK, outfits, &image, naxes1, naxes2);
             status =+ cf_fes_cal(FES_BIAS, outfits, &image, naxes1, naxes2);
             status =+ cf_fes_cal(FES_FLAT, outfits, &image, naxes1, naxes2);
             status =+ cf_fes_write(outfits, hdu, image); 


             free(image);

             if (status != 0)
        	cf_if_error("cf_calfes failed.");
	}

        /* close the FES FITS files */
	FITS_close_file(infits, &status);
	FITS_close_file(outfits, &status);


        /* Enter a timestamp into the log. */
        cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done Processing");

	return(status);
}
