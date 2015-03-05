/*****************************************************************************
 *	        Johns Hopkins University
 *	        Center For Astrophysical Sciences 
 *	        FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_hist_init input_file output_file 
 *
 * Description: Reads a raw hist file and generates a list of the pixels
 *              that contain counts. This list has the same format 
 *              as a TTAG list, so analysis can be carried out using
 *		many of the modules in the TTAG pipeline.
 *
 *              The output file will have the same format as a TTAG 
 *              file, with the following definitions:
 *
 *              TIME            FLOAT   all equal to the midpoint of the exp.
 *              XRAW            INT     raw x position 
 *              YRAW            INT     raw y position (midpoint of 8 pixel bin)
 *              PHA             BYTE    all 20 - indicating a good point
 *              WEIGHT          FLOAT   number of photons falling in that bin
 *              XFARF           FLOAT   x position in geometrically
 *                                      corrected frame
 *              YFARF           FLOAT   y position in geometrically
 *                                      corrected frame
 *              X               FLOAT   final x position (after correcting
 *                                      for motions)
 *              Y               FLOAT   final y position (after correcting
 *                                      for motions)
 *              CHANNEL         CHAR    aperture ID for the photon
 *              TIMEFLGS        CHAR    temporal flags - various bits set 
 *              LOC_FLGS        CHAR    location flags - various bits set 
 *              LAMBDA          FLOAT   wavelength for each photon
 *		ERGCM2		FLOAT	flux calibration for each photon
 *
 *              This routine calls cf_fuv_init to populate the calibration
 *              entries in the header.
 *
 *              In addition, the routine creates a timeline table containing
 *              the orbital information, count rates, and high-voltage
 *              information to be used by later pipeline modules.
 * 
 *
 * Returns:	Zero upon successful completion.
 *
 * History:	05/09/03   v1.1   rdr	Adapted from cf_ttag_init
 *		06/04/03   v1.4   wvd	Add call to cf_proc_check
 *              06/10/03   v1.5   rdr   Scan for hot pixels and remove 
 *		06/11/03   v1.6   wvd	Compress data with TSCALE and TZERO
 *		06/11/03   v1.7   wvd	Changed calfusettag.h to calfuse.h
 *		07/17/03   v1.8   wvd   Move various subroutines to libcf.
 *		07/22/03   v1.9   wvd   Read all extensions of raw hist file.
 *					Test for hot pixels only in active
 *					region and only if img[j]>15.
 *		07/28/03   v1.10  wvd   Close ELEC_CAL.
 *		10/26/03   v1.11  wvd   Increase morekeys to 38.
 *		11/05/03   v1.12  wvd   Added stdlib.h and unistd.h.
 *              12/05/03          bjg   Modified cf_make_pixel_list
 *                                      to unbin data in y.
 *              12/10/03   v1.13  bjg   Unbin is now optional 
 *                                      Update keyword SPECBINY
 *              03/03/04   v1.14  rdr   Exit if no photons found
 *              03/05/04   v1.15  rdr   Warning if data is missing
 *              03/09/04   v1.16  rdr   Fix error in correcting hot pixels
 *              03/26/04   v1.17  rdr   Continue when there is no data
 *              03/30/04   v1.20  wvd   Clean up i/o.
 *              03/31/04   v1.21  bjg   functions return EXIT_SUCCESS instead
 *                                      of zero. 
 *                                      >>>Include math.h for fmod<<<
 *                                      Type conversion for
 *                                      double fmod(double, double)
 *              04/06/04   v1.22  bjg   Remove unused variables
 *                                      Change format to match arg type
 *                                      in printf
 *                                      Remove static keyword in struct key 
 *                                      definition
 *                                      Add braces in TSCAL and TZERO
 *                                      definitions
 *                                      Test for fill data
 *              10/15/04   v1.23  bjg   Fill Data Photon gets flag
 *                                      LOCATION_FILL only if it is
 *                                      in Active Aperture
 *              12/24/2004 v1.24  bjg   npix_max now takes binx and biny into
 *                                      account
 *		01/24/2005 v1.25  wvd   Modify cf_make_pixel_list to ignore
 *					duplicate data, which can sneak in
 *					with the stim pulses.
 *		02/11/2005 v1.26  wvd   Hot pixels and fill data cause the
 *					header keyword NEVENTS to be too 
 *					large, so we correct/exclude these
 *					points and update NEVENTS.
 *		03/10/2005 v1.27  wvd   Move flagging of airglow photons to
 *					cf_screen_airglow.
 *					Initialize XFARF and YFARF with zero.
 *		05/20/2005 v1.28  wvd   Clean up i/o.
 *		05/31/2006 v1.29  wvd   OPUS has been updated, so we can set
 *					morekeys to 0.
 *		06/12/2006 v1.30  wvd   Call cf_verbose rather than 
 *					cf_if_warning when an extension 
 *					contains no photons.
 *		02/23/2007 v1.31  wvd   If a raw hist file contains no 
 *					extensions, set EXP_STAT to -1.
 *					Make active-area limits inclusive.
 *		11/07/2008 v1.32  wvd   Clean up i/o. 
 *
 ****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "calfuse.h" 

struct key {
        char keyword[FLEN_KEYWORD];
        float value;
        };


/*****************************************************************************/

static int
cf_make_pixel_list(fitsfile *infits, fitsfile *outfits, char unbin) {

/*****************************************************************************
 *
 * Description: procedure to take input data from the raw HIST data file and 
 *              generate a table with 14 columns. This table lists
 *              all pixels that contain counts. The
 *              table is designed to mimic a TTAG list, so that it can be used
 *              with various modules in the TTAG pipeline
 *
 ****************************************************************************/

    fitsfile *elecfits; 
    char  fmt_byte[FLEN_CARD], fmt_float[FLEN_CARD], fmt_short[FLEN_CARD];
    char  *inpha, *dumarrb;
    char  elecfile[FLEN_CARD];
    char  keyword[FLEN_CARD], card[FLEN_CARD];
    unsigned char *dumarrbu ;
    int   status=0, anynull=0, hdutype;
    int	  colval;
    int   tfields = 14;      /*  output table will have 14 columns */
    int   active_l, active_r ;
    long  i, j, k, nevents, nfdp;
    long  frow=1, felem=1, nrows=1; /* output table will have 1 row */
    float sum_weights;
    float *intime, *dumarrf;

    int   sia_index[512], sia_temp[512];  /* Indices of SIA rectangles */
    int   sia_width = 2048, sia_height = 16;  /* Dimensions of SIA recs */
    int   xx, yy;		/* Pixel coordinates (unbinned) */

    int x0, y0, xstart, ystart ;
    int nx, ny, yndx, binx, biny ;
    int nhot, nhot_total=0;
    short *xpix, *ypix, *img, *fill_data_pix, nmask, nullval;
    long npix, npix_max, npts ;
    float *wtpix, exptime, mtime, nave;
    

    char extname[]="TTAG DATA";   /* name of the extension containing 
                                        the binary table holding the data */

    char *ttype[]={"TIME", "XRAW", "YRAW", "PHA", "WEIGHT", "XFARF", 
		   "YFARF",  "X", "Y", "CHANNEL", "TIMEFLGS",
		   "LOC_FLGS", "LAMBDA", "ERGCM2" };

    char *tform[14];  /* will asign values when we find out 
                         the number of elements in the data set */


    char *tunit[]={"SECONDS", "PIXELS", "PIXELS", " ", " ", "PIXELS",
		   "PIXELS", "PIXELS", "PIXELS", "UNITLESS",
		   "UNITLESS", "UNITLESS", "ANGSTROMS","ERG CM^-2"};

    struct key tscal[] = {{"TSCAL6", 0.25},  {"TSCAL7", 0.1},
                          {"TSCAL8", 0.25},  {"TSCAL9", 0.1}};
    struct key tzero[] = {{"TZERO6", 8192.}, {"TZERO7", 0.},
                          {"TZERO8", 8192.}, {"TZERO9", 0.}};

    /* Initialize the sia_index array. */
    (void) memset((void *)sia_index, 0, sizeof(int)*512);
 
    /*********************************************************************
             Read information to be used in the screening analysis.
    **********************************************************************/

    /* Read the limits to the active area of the detector */
    FITS_read_key(outfits, TSTRING, "ELEC_CAL", elecfile, NULL, &status);
    FITS_open_file(&elecfits, cf_cal_file(elecfile), READONLY, &status);
    FITS_read_key(elecfits, TINT, "ACTIVE_L", &active_l, NULL,
		  &status);
    FITS_read_key(elecfits, TINT, "ACTIVE_R", &active_r, NULL,
		  &status);
    FITS_close_file(elecfits, &status);
    cf_verbose(3, "Limits to the active area: left=%d, right=%d",
	   active_l, active_r);

    /* Determine the exposure time */
    FITS_movabs_hdu(infits, 1, &hdutype, &status);
    FITS_read_key(infits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    cf_verbose(3, "exptime = %d ", cf_nint(exptime)) ;

    /* Determine the binning factors in x and y */
    FITS_read_key(infits, TINT, "SPECBINX", &binx, NULL, &status);
    FITS_read_key(infits, TINT, "SPECBINY", &biny, NULL, &status);
    cf_verbose(3, "binx=%d, biny=%d \n", binx, biny) ;


   /***********************************************************************
              Get information from the raw data file.
   ************************************************************************/

    /* Set up arrays to contain the positions of each pixel with counts 
        (xpix and ypix) and the number of counts in that pixel (wtpix) */
	if (unbin) npix_max = (NXMAX * NYMAX) / binx ;
	else npix_max = (NXMAX * NYMAX) / (binx * biny) ;
        npix = -1 ;
        xpix =  (short *) cf_calloc(npix_max, sizeof(short)) ;
        ypix =  (short *) cf_calloc(npix_max, sizeof(short)) ;
        wtpix = (float *) cf_calloc(npix_max, sizeof(float)) ;
	fill_data_pix = (short *) cf_calloc(npix_max, sizeof(short)) ;
         
    /*
     *  Loop over all extensions in the input file, and place contents
     *  of each extension into outbuff.  Keep trying to read a new HDU
     *  in the input file until we hit the EOF.
     *
     *  Note that the next line cannot use FITS_movrel_hdu because
     *  that routine automatically exits upon error.
     */
    for (i = 1; !(fits_movrel_hdu(infits, 1, &hdutype, &status)); i++) {
        FITS_read_key(infits, TINT, "NAXIS1",  &nx, NULL, &status);
        FITS_read_key(infits, TINT, "NAXIS2",  &ny, NULL, &status);
        FITS_read_key(infits, TINT, "XORIGIN", &x0, NULL, &status);
        FITS_read_key(infits, TINT, "YORIGIN", &y0, NULL, &status);

	npts = nx * ny ;
        cf_verbose(3, "Reading extension %d ", i) ;
        cf_verbose(3, "nx=%d, ny=%d, npts=%d ", nx, ny, npts) ;

	/* Initialize the sia_temp array */
	(void) memset((void *)sia_temp, 0, sizeof(int)*512);
 
        /*
         *  Allocate space for and read the input image
         */
        img = cf_malloc(sizeof(short) * npts);
        FITS_read_img(infits, TSHORT, 1L, npts, &nullval,
                      img, &anynull, &status);

	/* Specify the starting positions in x and y.
	    Place the point in the middle of the bin */
	xstart = x0 * binx + binx / 2;
        ystart = y0 * biny + biny / 2;
	cf_verbose(3, "x0=%d, y0=%d, xstart=%d, ystart=%d",
		x0, y0, xstart, ystart) ;
      
	/* Go through the image and tabulate the pixels with counts. */
	nhot=0 ;
        for (j=0; j<npts; j++) 
          if (img[j] > 0) {

	    /* If this pixel falls in an SIA rectangle that has already
		been written to the output data list, skip it and move to
		the next SIA rectangle. */

		xx = (x0 + (j % nx)) * binx;
		yy = (y0 + (j / nx)) * biny;
		k = xx / sia_width + yy / sia_height * 8;
		if (k < 0 || k >= 512) continue;
		if (sia_index[k] == 1) {
			j += sia_width / binx - 1;
			continue;
		}
		sia_temp[k] = 1;

	    /* Check for a hot pixel by comparing the value in the pixel with
               those on either side.  Ignore pixels outside of detector active
               area or with fewer than 16 counts. */

            if ((img[j] > 15) && (img[j] != FILL_DATA) && 
		(xx >= active_l) && (xx <= active_r)) {
                nave = (img[j-2] + img[j-1] + img[j+1] + img[j+2]) / 4. ;
		  if (nave > 32000.) nave=32000. ;
                  if (nave < 1.) nave=1. ;
	        if (img[j] > 8*nave) {
		  nhot++ ;
		  cf_verbose(3,"hot pixel at point %d: val=%d, ave=%f", 
                       j, img[j], nave) ;
		  /* Repair hot pixel */
                  nmask=1 ;
		  while(nmask < 8*nave && nmask < 16000) nmask*=2 ; 
		  nmask-=1 ;
		  img[j] = (img[j] & nmask) ;
                  cf_verbose(3,"corrected value = %d ",img[j]) ;
		}
	    }

	    yndx = (int) ( j / nx ) ;

            if (unbin){
	      for (k=0;k<biny;k++){
		npix++ ;
		xpix[npix]  = (short) (xstart + (j - (yndx * nx)) * binx) ;
		ypix[npix]  = (short) (ystart + yndx*biny + k) ;
		wtpix[npix] = (float) ( ((float) img[j]) / ((float) biny) ) ;

		fill_data_pix[npix]=0;
		if (img[j] == FILL_DATA) {
		  wtpix[npix]=0.0;
		  fill_data_pix[npix]=1;
		}


	      }
	    }
	    else {
	      npix++ ;
	      xpix[npix]  = (short) (xstart + (j - (yndx * nx)) * binx) ;
	      ypix[npix]  = (short) (ystart + yndx*biny) ;
	      wtpix[npix] = (float) img[j];
	      
	      fill_data_pix[npix]=0;
	      if (img[j] == FILL_DATA) {
		wtpix[npix]=0.0;
		fill_data_pix[npix]=1;
	      }
	    }

	   
	  }

	nhot_total += nhot;
        cf_verbose(2, "%d hot pixels found in extn %d  ",nhot, i) ;
        cf_verbose(3, "total number of pixels with counts = %d \n", npix) ;
        if (nhot > 10) cf_if_warning("More than 10 hot pixels found.") ;

        free(img) ;

	/* Copy new SIA indices to sia_index array. */
	for (k = 0; k < 512; k++)
		if (sia_temp[k]) sia_index[k] = 1;

    }  /* Finished loop over all extensions in infits */
    /*
     *  Status upon completion of loop through extensions should be
     *  "end of file."
     *  If it is, reset value of status.  If not, there was an
     *  unexpected error.  Print and exit.
     */
    if( status == END_OF_FILE ) {  /* !!! This might cause a problem */
        status = 0;
    } else  {
        cf_if_fits_error(status);
    }

    cf_verbose(3, "Finished reading data ") ;

    if (npix < 0) {
	npix = 1 ;
        xpix[0]=0 ;
        ypix[0]=0 ;
        wtpix[0]=1. ;
	fill_data_pix[0]=0;
    }

   /*************************************************************************
           Output the analysis to a table
   ***************************************************************************/  

    /* Generate the tform array */
    sprintf(fmt_byte,  "%ldB", npix);
    sprintf(fmt_float, "%ldE", npix);
    sprintf(fmt_short, "%ldI", npix);

    tform[0] = fmt_float;
    tform[1] = fmt_short;
    tform[2] = fmt_short;
    tform[3] = fmt_byte;
    tform[4] = fmt_float;
    for (i=5; i<9; i++)
        tform[i] = fmt_short;
    for ( ; i<12; i++)
        tform[i] = fmt_byte;
    tform[12] = fmt_float;
    tform[13] = fmt_float;

   /* Append a new empty binary table to the output file */

    cf_verbose(3, "Creating table") ;
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);
    FITS_create_tbl(outfits, BINARY_TBL, nrows, tfields, ttype, tform,
                        tunit, extname, &status);

    /* Write TSCALE and TZERO entries to header */
    for (i=0; i<4; i++) {
        sprintf(keyword, "TUNIT%ld", i+6);
        if (fits_read_keyword(outfits, keyword, card, NULL, &status))
		cf_if_fits_error(status);
        FITS_insert_key_flt(outfits, tscal[i].keyword, tscal[i].value, -2,
		NULL, &status);
        FITS_insert_key_flt(outfits, tzero[i].keyword, tzero[i].value, -4,
		NULL, &status);
    }
   
    /*   Allocate space for the new output columns */
  
    cf_verbose(3, "Allocating space for arrays ") ;
    dumarrf  = (float *) cf_malloc(sizeof(float) * npix);
    dumarrb  = (char *) cf_malloc(sizeof(char) * npix);
    dumarrbu  = (unsigned char *) cf_malloc(sizeof(unsigned char) * npix);
 
    cf_verbose(3, "Filling columns") ;  
    /* Specify a time array - assign the same time to all pixels */
    intime = (float *) cf_calloc(npix, sizeof(float) ) ;
    mtime = exptime/2. ;
    for (i=0; i<npix; i++) intime[i]=mtime ;
    FITS_write_col(outfits, TFLOAT, 1, frow, felem, npix, intime, &status);

    /* Output X and Y positions  */
    FITS_write_col(outfits, TSHORT, 2, frow, felem, npix, xpix, &status);
    FITS_write_col(outfits, TSHORT, 3, frow, felem, npix, ypix, &status);

    /* Now do the PHA column - assume 20 */
    inpha = (char *) cf_calloc(npix, sizeof(char) ) ;
    for (i=0; i<npix; i++) inpha[i] = 20 ;
    FITS_write_col(outfits, TBYTE, 4, frow, felem, npix, inpha, &status);

    /* Fill weights column with counts per pixel */ 
    FITS_write_col(outfits, TFLOAT, 5, frow, felem, npix, wtpix, &status);

    /* Fill XFARF, YFARF, X and Y arrays with 0 */
    for (i=0; i<npix; i++) dumarrf[i] = 0.;
    for (colval=6; colval<10; colval++) {
        FITS_write_col(outfits, TFLOAT, colval, frow, felem, npix,
		dumarrf, &status);
    }

    /* Initialize the channel and temporal flag arrays with 0 - byte */
    for (i=0;i<npix; i++)
	dumarrb[i] = 0;
    FITS_write_col(outfits, TBYTE, 10, frow, felem, npix, dumarrb, &status);
    for (i=0;i<npix; i++)
	dumarrbu[i] = 0;
    FITS_write_col(outfits, TBYTE, 11, frow, felem, npix, dumarrbu, &status);

    /* Flag events that lie outside of the active area or
	in fill-data regions. */
    for (i=0 ; i<npix ; i++) {
	dumarrbu[i] = 0 ;

	/* Flag photons that lie outside the active area. */
	if (xpix[i] < active_l || xpix[i] > active_r) 
	    dumarrbu[i] = dumarrbu[i] | LOCATION_SHLD ;

	/* If inside, is it fill data? */
	else if (fill_data_pix[i])
	    dumarrbu[i] = dumarrbu[i] | LOCATION_FILL ; 
    }

    FITS_movabs_hdu(outfits, 2, &hdutype, &status);
    FITS_write_col(outfits, TBYTE, 12, frow, felem, npix,
		       dumarrbu, &status);

    /* Fill lambda and energy arrays with 0 */
    FITS_write_col(outfits, TFLOAT, 13, frow, felem, npix, dumarrf, &status);
    FITS_write_col(outfits, TFLOAT, 14, frow, felem, npix, dumarrf, &status);

    /* If file contains hot pixels or fill data, recalculate NEVENTS. */
    nfdp = 0;
    sum_weights = 0.;
    for (i=0; i<npix; i++) {
	nfdp += fill_data_pix[i];
	sum_weights += wtpix[i];
    }
    if (nfdp > 0 || nhot_total > 0) {
	nevents = cf_nlong(sum_weights);
	FITS_movabs_hdu(outfits, 1, &hdutype, &status);
	FITS_update_key(outfits, TLONG, "NEVENTS", &nevents, NULL, &status);
	cf_verbose(1, "Image contains bad pixels.  Updating NEVENTS = %ld",
	nevents);
    }

    free(intime);
    free(xpix);
    free(ypix);
    free(wtpix) ;
    free(inpha);
    free(dumarrf);
    free(dumarrb);

    cf_verbose(3, "Leaving cf_make_pixel_list") ;

    return EXIT_SUCCESS ;
}


/******************************************************************************
 *
 *                   CF_MAKE_GTI_TABLE
 *
 *      Procedure to generate a table of GTI values and put it into the second
 *      extension of the output file. There is assumed to be only one GTI, 
 *      starting at 0 and ending at the exposure time.
 *
 ******************************************************************************/

int cf_make_gti_table (fitsfile *outfits) {

    int status=0, nrows=1, hdutype, frow=1, felem=1  ;
    float exptime ;
    double start[1]={0.}, stop[1] ;

    char extname[]="GTI";
    int tfields=2;	/* output table will have 2 columns */
    char *ttype[]={"START", "STOP"};

    char *tform[]={"1D","1D"}; 
    char *tunit[]={"seconds", "seconds"};

    cf_verbose(3, "Entering cf_make_gti_table ") ;

    /* read the exposure time from the output file */
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);
    FITS_read_key(outfits, TFLOAT, "EXPTIME", &exptime, NULL, &status);
    stop[0] = (double) exptime ;


    /* Append a new empty binary table to the output file */
    FITS_movabs_hdu(outfits, 2, &hdutype, &status);
    FITS_create_tbl(outfits, BINARY_TBL, nrows, tfields, ttype, tform,
			tunit, extname, &status);

    /* Write the data to the table */
    FITS_write_col(outfits, TDOUBLE, 1, frow, felem, nrows, start, &status);
    FITS_write_col(outfits, TDOUBLE, 2, frow, felem, nrows, stop, &status);

    return EXIT_SUCCESS ;
}

/*****************************************************************************/


int main(int argc, char *argv[])
{
    char CF_PRGM_ID[] = "cf_hist_init"; 
    char CF_VER_NUM[] = "1.32";

    fitsfile *infits, *outfits;
   
    int   morekeys=0;		/* Change to zero when OPUS is updated. */
    int   status=0, optc;
   
    int one=1, nextn=0;
    long nphot ;
    char unbin=0;

    char opts[] = "hsv:";
    char usage[] =
	"Usage:\n"
	"  cf_hist_init [-h] [-s] [-v level] rawhistfile idf_file\n";
    char option[] = 
	"Options:\n"
	"  -h:  this help message\n"
        "  -s:  to unbin the data in Y\n"
        "  -v:  verbosity level (default is 1; 0 is silent)\n";
      

    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 's':
	    unbin=1;
	    break;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	}
    }
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    if (argc != optind+2) {
        printf("%s", usage);
        cf_if_error("Incorrect number of arguments");
    }
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin execution");

    FITS_open_file(&infits, argv[optind++], READONLY, &status);

    /*  Exit if data are not in HIST mode. */
    if ( cf_proc_check(infits, CF_PRGM_ID) ) {
        FITS_close_file(infits, &status);
        return (-1);
    }

    /* Issue warning message if data are missing. */
    FITS_read_key(infits, TLONG, "NEVENTS", &nphot, NULL, &status); 
    FITS_read_key(infits, TINT, "NEXTEND", &nextn, NULL, &status); 
    cf_verbose(3,"Number of extensions to the file = %d ", nextn) ;
    if (nphot < 1) cf_verbose(1, "%s contains no photons.", argv[optind-1]);
    if (nextn < 1) cf_if_warning("%s contains no extensions.", argv[optind-1]); 
    else if (nextn < 4) cf_verbose(1, "%s contains only %d extensions.",
	argv[optind-1], nextn) ; 

    /*  Create output file. */
    FITS_create_file(&outfits, argv[optind], &status);
    /*  Copy primary image header from the input file. */
    FITS_copy_hdu(infits, outfits, morekeys, &status);

    /* Update header keywords. */
    FITS_update_key(outfits, TSTRING, "FILENAME", argv[optind], NULL, &status);
    FITS_update_key(outfits, TSTRING, "FILETYPE",
                 "INTERMEDIATE DATA FILE", NULL, &status);

    if (unbin) FITS_update_key(outfits, TINT, "SPECBINY", &one, NULL, &status);

    /* If input file contains no extensions, set EXP_STAT to -1. */
    if (nextn < 1) {
	char comment[] = "Raw file contains no data.";
	int eflag = -1; 
	FITS_update_key(outfits, TINT, "EXP_STAT", &eflag, comment, &status);
    }

    /* Populate the calibration keywords in the primary header. */
    cf_fuv_init(outfits);

    /* Fill in the first extension with the photon information */
    cf_make_pixel_list(infits, outfits, unbin) ;

    /* Generate a good-times table and put it in the second extension. */
    cf_make_gti_table(outfits);
  
    /*  Append a table containing the orbital timeline */
    cf_timeline(outfits);
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /*  Update processing flags. */
    cf_proc_update(outfits, CF_PRGM_ID, "COMPLETE");

    FITS_close_file(infits, &status);
    FITS_close_file(outfits, &status);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");

    return EXIT_SUCCESS ;
}
