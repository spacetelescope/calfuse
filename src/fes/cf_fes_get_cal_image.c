 /*******************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *******************************************************************************
 *
 * Synopsis:int cf_fes_get_cal_image(char cal_file_name[], float **cal_image, 
 *                                                  int naxis1, int naxis2)
 *
 * Description: open the specified calibration file for an FES image
 *              loop through each extension in the calibration file
 *             until you find the one that has the same size naxis1 and
 *              naxis2 as the input image.
 *              read the image in that extention into *cal_image
 *              close the calibration file
 *
 * Arguments:	char cal_file_name[], 
 *              float **cal_image, 
 *              int naxis1, int naxis2
 *
 * Returns:	none
 *
 * History:	04/14/99        mlr     started work
 *              04/16/99        mlr     modified to use libcf and FITSIO.h
 *
 ******************************************************************************/

#include "calfuse.h"

#define CF_PRGM_ID      "cf_fes_get_cal_image"
#define CF_VER_NUM      "1.4"

int cf_fes_get_cal_image(char cal_file_name[], float **cal_image, 
                                  int naxis1, int naxis2)
{
      fitsfile *fes_cal_file;
      int status;
      int i, hdu, hdutype, num_exts;
      int cal_axis1, cal_axis2;
      char buffer[FLEN_CARD];
 
      cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin getting Cal image");

      status = 0;

      /* open the input FES CALIBRATION file */

      FITS_open_file(&fes_cal_file,cf_cal_file(cal_file_name),READONLY,&status);

      FITS_read_key(fes_cal_file, TINT, "NEXTEND", &num_exts, buffer, &status);

      for (i = 1; i<= num_exts; i++)
      {
 	  hdu = i + 1;
          FITS_movabs_hdu(fes_cal_file, hdu, &hdutype, &status);
          cf_fes_read(fes_cal_file, cal_image, &cal_axis1, &cal_axis2);
          if ((cal_axis1 == naxis1) && (cal_axis2 == naxis2))
                   break;
      }
      FITS_close_file(fes_cal_file, &status);

      if (i>num_exts) 
      { 
          *cal_image=NULL;
          cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "No Cal image of requested size");
          return(-1);
      }
       
      cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done getting Cal image");
      return(status);
}
