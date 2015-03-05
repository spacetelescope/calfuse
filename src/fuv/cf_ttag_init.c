/*****************************************************************************
 *	        Johns Hopkins University
 *	        Center For Astrophysical Sciences 
 *	        FUSE
 *****************************************************************************
 *
 * Synopsis:	cf_ttag_init input_file output_file 
 *
 * Description: Reads a raw ttag file, which is a FITS binary table
 *              consisting of four columns:
 *
 *              TIME            FLOAT 
 *              X               INT
 *              Y               INT
 *              PHA             BYTE
 *
 *              This routine will convert the raw ttag file to another binary
 *              table with 13 columns the values for most of these columns
 *              will be added by various modules of the pipeline:
 *
 *              TIME            FLOAT   time of arrival for photon
 *              XRAW            INT     raw x position
 *              YRAW            INT     raw y position
 *              PHA             BYTE    pha value
 *              WEIGHT          FLOAT   effective area and flux calibration
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
 *              This routine calls cf_fuv_init to set up the calibration
 *              entries in the header.
 *
 *              In addition, the routine creates a timeline containing
 *              orbital information, count rates, and high-voltage
 *              information which is used in later pipeline modules.
 *
 *		Turns out that OPUS good-time intervals may exclude SAA
 *		passes that remain in the photon-event list.  This means
 *		that there may be photons in the photon-event list whose
 *		times are not listed in the timeline table.  To catch
 *		them, we set all photon status flags to TEMPORAL_HV.
 * 
 *
 * Returns:	Zero upon successful completion
 *
 * History:	03/30/98	gak	Begin work
 *              03/18/99        emm     Modified header.
 *              04/06/99        emm     Added check on nrows to be no
 *                                        larger than nevents.
 *              04/15/99    barrett     Deleted FF and DY columns,
 *                                        added FITS wrapper functions,
 *                                        tidied code.
 *              04/22/99    emm         Program now returns int.
 *              04/23/99    emm         Fixed two bugs causing segmentation
 *                                      faults (pha was being written as
 *                                      long, at definition of outdx was
 *                                      given as 1B instead of 1E.) Added
 *                                      check on number of arguments.
 *              04/29/99    emm         Added check on nrows to make sure
 *                                      it is smaller than nevents.
 *              04/30/99    peb         Added cf_malloc functions.
 *              06/07/99    peb         Added reporting of version number.
 *              18/07/99 1.12    jm     Populate APER_ACT keyword
 *              10/28/99 1.14   emm     Added DY column.
 *              12/02/99 1.15   emm     Added corrections to some keywords
 *                                      removed APER_ACT keyword update
 *              12/21/00 1.16   peb     Fixed type mismatches in cfitsio
 *                                      calls
 *              05/16/01 1.17   rdr     added night/day flag column (dnflg)
 *              08/20/01 1.18   wvd	Check keywords SPECBINX and SPECBINY;
 *					if 0, set to 1
 *              10/19/01 1.19   wvd	Move correction of SPECBINX, etc.
 *					from cf_ttag_image.
 *              04/03/02 1.20   wvd	Clean up keywords in photon-event list.
 *              08/02/02 1.21   rdr     - changed the format of the binary
 *                                        tables to improve speed
 *                                      - added new columns to the table
 *              09/15/02 1.22   rdr     added extension containing timeline
 *              11/07/02 1.23   peb     Modified columns for event and
 *                                      timeline flags.  Removed temporary
 *                                      keywords.  Changed include files to
 *                                      use calfusettag.h
 *              11/13/02 1.24   rdr     Added ycentl and ycents columns to
 *					timeline.
 *                                      Added section to flag geocoronal photons
 *                                      and photons in the inactive part of the 
 *                                      detector. Moved processing of the 
 *                                      photon list into a separate subroutine.
 *                                      Fixed problem of program crashing when
 *                                      needed columns are not in the hskp file
 *              12/06/02 1.25   rdr     removed variables which were not being 
 *					used
 *		12/12/02 1.5	wvd	Installed program.
 *              12/18/02 1.6    rdr     use AIRG calibration file for airglow
 *                                      locations
 *              12/20/02 1.7    rdr     - check for negative count rates and
 *                                        correct
 *                                      - convert status flags to unsigned char
 *                                      - flag photons near start and end of
 *                                        GTI
 *              01/24/03 1.8    rdr     - correct situation where hskp file
 *                                        has blank columns
 *              02/11/03 1.9    rdr     Added ERGCM2S column to photon table
 *              03/04/03 1.10   rdr     removed YFARFC column and deleted some
 *                                      unused variables
 *              03/19/03 1.12   peb     Added command-line options and
 *                                      verbose_level parameter
 *                                      Changed all printf calls to cf_verbose
 *                                      calls with verbose level == 1
 *              03/25/03 1.13   peb     Changed a cf_verbose to cf_timestamp.
 *              04/01/03 1.15   wvd     Call cf_error_init after cf_timeline.
 *              04/04/03 1.16   wvd     Update FILETYPE keyword. Set verbose
 *					level of each call to cf_verbose.
 *              04/07/03 1.17   wvd     Update FILENAME keyword.
 *              04/07/03 1.18   wvd     Add FPADXLIF and FPADXSIC to header.
 *              04/09/03 1.19   rdr     Correct error in determining the name
 *                                      of the housekeeping file
 *              04/29/03 1.21   wvd	Allow possiblity that housekeeping
 *                                      filename is in lower-case letters.
 *              05/06/03 1.22   rdr     Add daynight keyword, bkg_max and
 *                                      bkg_min keywords, and changed the 
 *                                      column ergcm2s to ergcm2
 *              05/10/03 1.23   wvd	Modify internal copy of GTI's.
 *					Make number of timeline entries more
 *				 	nearly equal to EXPTIME.  Discard
 *					photons with bad time markers.
 *              06/03/03 1.24   wvd	Add call to cf_proc_check.
 *					Keep trying to get GTI's right.
 *              06/09/03 1.25   wvd	Compress data with TSCALE and TZERO.
 *		06/11/03 1.26   wvd     Change calfusettag.h to calfuse.h
 *              07/07/03 1.27   rdr     Correct initialization of dn_flag_last
 *                                      when searching for time_sunrise and
 *                                      time_sunset 
 *		07/17/03 1.28   wvd     Move various subroutines to library.
 *		07/23/03 1.29   wvd     Close ELEC_CAL file.
 *					Delete goto statement.
 *		09/15/03 1.31   wvd     Set temporal flags to TEMPORAL_HV.
 *		10/02/03 1.32   wvd     Delete reference to LOCATION_GTI flag
 *		10/26/03 1.33   wvd     Increase to morekeys to 38.
 *		11/05/03 1.34   wvd     Added stdlib.h and unistd.h
 *              03/03/04 1.35   rdr     Exit if no photons found.
 *              04/06/04 1.37   bjg     Change formats to match arg types 
 *                                      in printf
 *                                      Functions return EXIT_SUCCESS
 *                                      Remove static keyword in struct key 
 *                                      definition
 *                                      Add braces in TSCAL and TZERO
 *                                      definitions
 *                                      When no events are present in
 *                                      the raw file, create a photon list
 *                                      with a dummy event.
 *		03/10/05 1.38   wvd	Move airglow screening to 
 *					cf_screen_airglow.
 *					Initialize XFARF and YFARF to zero.
 *		05/20/05 1.39   wvd	Clean up i/o.
 *		05/31/06 1.40   wvd     OPUS has been updated, so we can set
 *					morekeys to 0.
 *		06/12/06 1.41   wvd     Change cf_if_warning to cf_verbose.
 *		08/21/07 1.42   wvd     Add -e flag to set HKEXISTS to YES.
 *
 ****************************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include "calfuse.h" 

struct key {
	char keyword[FLEN_KEYWORD];
	float value;
	};

static int
cf_init_null_photon_list(fitsfile *outfits) {
       
    int   status=0;

    long i;

    char  keyword[FLEN_CARD], card[FLEN_CARD];
      
    char extname[]="TTAG DATA";   /* Name of this extension */

    char *ttype[]={"TIME", "XRAW", "YRAW", "PHA", "WEIGHT", "XFARF", 
		   "YFARF",  "X", "Y", "CHANNEL", "TIMEFLGS",
		   "LOC_FLGS", "LAMBDA", "ERGCM2" };

    char *tform[14]={"1E","1I","1I","1B","1E","1I","1I",
		     "1I","1I","1B","1B","1B","1E","1E"};  

    char *tunit[]={"SECONDS", "PIXELS", "PIXELS", "UNITLESS", "UNITLESS",
		   "PIXELS", "PIXELS", "PIXELS", "PIXELS", "UNITLESS",
		   "UNITLESS", "UNITLESS", "ANGSTROMS", "ERG CM^-2"};

    struct key tscal[] = {{"TSCAL6", 0.25},  {"TSCAL7", 0.1}, 
                          {"TSCAL8", 0.25},  {"TSCAL9", 0.1}};
    struct key tzero[] = {{"TZERO6", 8192.}, {"TZERO7", 0.}, 
                          {"TZERO8", 8192.}, {"TZERO9", 0.}};

    float   hdu2_time[1]      =  {0.0};
    short   hdu2_xraw[1]      =  {0};
    short   hdu2_yraw[1]      =  {0};
    char    hdu2_pha[1]       =  {0};
    float   hdu2_weight[1]    =  {0.0};
    float   hdu2_xfarf[1]     =  {0.0};
    float   hdu2_yfarf[1]     =  {0.0}; 
    float   hdu2_x[1]         =  {0.0};
    float   hdu2_y[1]         =  {0.0};
    char    hdu2_channel[1]   =  {0};
    char    hdu2_timeflgs[1]  =  {TEMPORAL_HV};
    char    hdu2_loc_flgs[1]  =  {LOCATION_SHLD};
    float   hdu2_lambda[1]    =  {0.0};
    float   hdu2_ergcm2s[1]   =  {0.0};

    cf_verbose(3, "Entering cf_init_null_photon_list.");

    /* Append a new empty binary table to the output file */
    FITS_create_tbl(outfits, BINARY_TBL, 1, 14, ttype, tform,
			tunit, extname, &status);

    /* Write TSCALE and TZERO entries to header */
    for (i=0; i<4; i++) {
	sprintf(keyword, "TUNIT%ld", i+6);
	if (fits_read_keyword(outfits, keyword, card, NULL, &status))
		cf_if_fits_error(status);
	FITS_insert_key_flt(outfits, tscal[i].keyword, tscal[i].value,
		-2, NULL, &status);
	FITS_insert_key_flt(outfits, tzero[i].keyword, tzero[i].value,
		-4, NULL, &status);
    }



  
    FITS_write_col(outfits, TFLOAT,  1, 1, 1, 1, hdu2_time     , &status);
    FITS_write_col(outfits, TSHORT,  2, 1, 1, 1, hdu2_xraw     , &status);
    FITS_write_col(outfits, TSHORT,  3, 1, 1, 1, hdu2_yraw     , &status);
    FITS_write_col(outfits, TBYTE ,  4, 1, 1, 1, hdu2_pha      , &status);
    FITS_write_col(outfits, TFLOAT,  5, 1, 1, 1, hdu2_weight   , &status);
    FITS_write_col(outfits, TFLOAT,  6, 1, 1, 1, hdu2_xfarf    , &status);
    FITS_write_col(outfits, TFLOAT,  7, 1, 1, 1, hdu2_yfarf    , &status);
    FITS_write_col(outfits, TFLOAT,  8, 1, 1, 1, hdu2_x        , &status);
    FITS_write_col(outfits, TFLOAT,  9, 1, 1, 1, hdu2_y        , &status);
    FITS_write_col(outfits, TBYTE , 10, 1, 1, 1, hdu2_channel  , &status);
    FITS_write_col(outfits, TBYTE , 11, 1, 1, 1, hdu2_timeflgs , &status);
    FITS_write_col(outfits, TBYTE , 12, 1, 1, 1, hdu2_loc_flgs , &status);
    FITS_write_col(outfits, TFLOAT, 13, 1, 1, 1, hdu2_lambda   , &status);
    FITS_write_col(outfits, TFLOAT, 14, 1, 1, 1, hdu2_ergcm2s  , &status);

    cf_verbose(3, "Exiting cf_init_null_photon_list.");
    return EXIT_SUCCESS;
}

/*****************************************************************************/

static int
cf_init_photon_list(fitsfile *infits, fitsfile *outfits) {

/*****************************************************************************
 *
 * Description: procedure to take input data from the raw data file and 
 *              generate a table with 14 columns. These columns
 *              will be modified by later pipeline modules.
 *
 * ARGUMENTS:
 *              infits        fitsfile  name of the input raw data file
 *              outfits       fitsfile  name of the output intermediate data
 *                                      file
 * INPUT DATA: 
 *              TIME            FLOAT   time of arrival for each photon
 *              XRAW            INT     x position of photon
 *              YRAW            INT     y position of photon
 *              PHA             BYTE    pha value of the photon
 *
 * OUTPUT DATA: 
 *              TIME            FLOAT   time of arrival for photon
 *              XRAW            INT     raw x position
 *              YRAW            INT     raw y position
 *              PHA             BYTE    pha value
 *              WEIGHT          FLOAT   effective area and flux calibration
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
 *		ERGCM2		FLOAT	flux cal for each photon
 *
 ****************************************************************************/

    fitsfile *elecfits; 
    char  fmt_byte[FLEN_CARD], fmt_float[FLEN_CARD], fmt_short[FLEN_CARD];
    char  *inpha;
    char  elecfile[FLEN_FILENAME];
    char  keyword[FLEN_CARD], card[FLEN_CARD];
    unsigned char *channel, *location, *temporal;
    short *inx, *iny ;
    int   status=0, intnull=0, anynull, hdutype;
    int	  tcol, xcol, ycol, phacol, colval;
    int   tfields = 14;      /*  output table will have 14 columns */
    int   active_l, active_r;
    long  i;
    long  frow=1, felem=1, nevents, nrows=1; /* output table will have 1 row */
    float *intime, *dumarrf;

    char extname[]="TTAG DATA";   /* Name of this extension */

    char *ttype[]={"TIME", "XRAW", "YRAW", "PHA", "WEIGHT", "XFARF", 
		   "YFARF",  "X", "Y", "CHANNEL", "TIMEFLGS",
		   "LOC_FLGS", "LAMBDA", "ERGCM2" };

    char *tform[14];  /* We'll assign values when we know 
                         the number of elements in the data set. */


    char *tunit[]={"SECONDS", "PIXELS", "PIXELS", "UNITLESS", "UNITLESS",
		   "PIXELS", "PIXELS", "PIXELS", "PIXELS", "UNITLESS",
		   "UNITLESS", "UNITLESS", "ANGSTROMS", "ERG CM^-2"};

    struct key tscal[] = {{"TSCAL6", 0.25},  {"TSCAL7", 0.1}, 
                          {"TSCAL8", 0.25},  {"TSCAL9", 0.1}};
    struct key tzero[] = {{"TZERO6", 8192.}, {"TZERO7", 0.}, 
                          {"TZERO8", 8192.}, {"TZERO9", 0.}};

    cf_verbose(3, "Entering cf_init_photon_list.");

    /*  Move to the extension of the input file containing the data. */
    FITS_movabs_hdu(infits, 2, &hdutype, &status);

    /* Read the limits to the active area of the detector. */
    FITS_read_key(outfits, TSTRING, "ELEC_CAL", elecfile, NULL, &status);
    FITS_open_file(&elecfits, cf_cal_file(elecfile), READONLY, &status);
    FITS_read_key(elecfits, TINT, "ACTIVE_L", &active_l, NULL, &status);
    FITS_read_key(elecfits, TINT, "ACTIVE_R", &active_r, NULL, &status);
    FITS_close_file(elecfits, &status);
    cf_verbose(3, "Limits to the active area: left=%d, right=%d",
	   active_l, active_r);

    /*
     *  Get the number of events in the input file.
     */
    FITS_read_key(infits, TLONG, "NAXIS2", &nevents, NULL, &status);
    cf_verbose(3, "Raw data file contains %ld entries", nevents);

    /* Generate the tform array */
    sprintf(fmt_byte,  "%ldB", nevents);
    sprintf(fmt_float, "%ldE", nevents);
    sprintf(fmt_short, "%ldI", nevents);

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
    FITS_create_tbl(outfits, BINARY_TBL, nrows, tfields, ttype, tform,
			tunit, extname, &status);

    /* Write TSCALE and TZERO entries to header */
    for (i=0; i<4; i++) {
	sprintf(keyword, "TUNIT%ld", i+6);
	if (fits_read_keyword(outfits, keyword, card, NULL, &status))
		cf_if_fits_error(status);
	FITS_insert_key_flt(outfits, tscal[i].keyword, tscal[i].value,
		-2, NULL, &status);
	FITS_insert_key_flt(outfits, tzero[i].keyword, tzero[i].value,
		-4, NULL, &status);
    }

    /*
     *  Allocate space for the input times, coords, and pha.
     */
    intime = (float *) cf_malloc(sizeof(float) * nevents);
    inx    = (short *) cf_malloc(sizeof(short) * nevents);
    iny    = (short *) cf_malloc(sizeof(short) * nevents);
    inpha  = (char *)  cf_malloc(sizeof(char)  * nevents);
    /*
     *  Allocate space for the new output columns.
     */
    dumarrf  = (float *) cf_malloc(sizeof(float) * nevents);
    channel  = (unsigned char *) cf_calloc(nevents, sizeof(unsigned char));
    location = (unsigned char *) cf_calloc(nevents, sizeof(unsigned char));
    temporal = (unsigned char *) cf_calloc(nevents, sizeof(unsigned char));

    /*
     *  Read the data from the input file and write to the output
     *  file in a new format
     *
     *  read the times
     */
    FITS_get_colnum(infits, TRUE, "TIME", &tcol, &status);
    FITS_read_col(infits, TFLOAT, tcol, frow, felem, nevents, &intnull,
		  intime, &anynull, &status);
    FITS_write_col(outfits, TFLOAT, 1, frow, felem, nevents, intime, &status);

    /* process X and Y positions  */
    FITS_get_colnum(infits, TRUE, "X", &xcol, &status);
    FITS_read_col(infits, TSHORT, xcol, frow, felem, nevents, &intnull,
		  inx, &anynull, &status);
    FITS_write_col(outfits, TSHORT, 2, frow, felem, nevents, inx, &status);
  
    FITS_get_colnum(infits, TRUE, "Y", &ycol, &status);
    FITS_read_col(infits, TSHORT, ycol, frow, felem, nevents, &intnull,
		  iny, &anynull, &status);
    FITS_write_col(outfits, TSHORT, 3, frow, felem, nevents, iny, &status);

    /* now do the PHA column */
    FITS_get_colnum(infits, TRUE, "PHA", &phacol, &status);
    FITS_read_col(infits, TBYTE, phacol, frow, felem, nevents, &intnull,
		  inpha, &anynull, &status);
    FITS_write_col(outfits, TBYTE, 4, frow, felem, nevents, inpha, &status);

    /* fill weights with a dummy float - all equal to 1. */ 
    for (i=0; i<nevents; i++)
	dumarrf[i] = 1.;
    FITS_write_col(outfits, TFLOAT, 5, frow, felem, nevents, dumarrf, &status);

    /* fill XFARF, YFARF, X and Y arrays with 0. */
    for (i=0; i<nevents; i++)
	dumarrf[i] = 0.;
    for (colval=6; colval<10; colval++) {
        FITS_write_col(outfits, TFLOAT, colval, frow, felem, nevents,
		       dumarrf, &status);
    }

    /* Default value of temporal array is TEMPORAL_HV */
    for (i=0; i<nevents; i++)
	temporal[i] = TEMPORAL_HV;
    FITS_write_col(outfits, TBYTE, 10, frow, felem, nevents, channel, &status);
    FITS_write_col(outfits, TBYTE, 11, frow, felem, nevents, temporal, &status);

    /* Set flags in the screening array for photons outside of the active
	region of the detector. */
    for (i=0; i<nevents; i++) 
	if (inx[i] < active_l || inx[i] > active_r) 
	    location[i] = location[i] | LOCATION_SHLD;

    /* Write location flags to output file. */
    FITS_movabs_hdu(outfits, 2, &hdutype, &status);
    FITS_write_col(outfits, TBYTE, 12, frow, felem, nevents,
		       location, &status);

    /* Fill lambda and energy arrays with 0. */
    FITS_write_col(outfits, TFLOAT, 13, frow, felem, nevents, dumarrf,
	   &status);
    FITS_write_col(outfits, TFLOAT, 14, frow, felem, nevents, dumarrf,
	   &status);

    free(intime);
    free(inx);
    free(iny);
    free(inpha);
    free(dumarrf);
    free(channel);
    free(location);
    free(temporal);

    cf_verbose(3, "Exiting cf_init_photon_list.");
    return EXIT_SUCCESS;
}


/*****************************************************************************/

int main(int argc, char *argv[])

/*****************************************************************************/
{
    char CF_PRGM_ID[] = "cf_ttag_init"; 
    char CF_VER_NUM[] = "1.42";

    fitsfile *infits, *outfits;
    int   morekeys=0;
    int   modify_hkexists=FALSE;
    int	  status=0, hdutype, optc;
    long  nphot ;

    char opts[] = "ehv:";
    char usage[] =
	"Usage:\n"
	"  cf_ttag_init [-he] [-v level] rawttagfile idf_file\n";
    char option[] = 
	"Options:\n"
	"  -h:  this help message\n"
	"  -e:  set header keyword HKEXISTS to YES in IDF file\n"
	"  -v:  verbosity level (default is 1; 0 is silent)\n";
 
    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
        case 'e':
            modify_hkexists = TRUE;
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

    /*  Exit if data are not in TTAG mode. */
    if ( cf_proc_check(infits, CF_PRGM_ID) ) {
	FITS_close_file(infits, &status);
	return (-1);
    }

    /* Issue warning if there are no events in the data set. */
    FITS_read_key(infits, TLONG, "NEVENTS", &nphot, NULL, &status); 
    cf_verbose(3,"Number of events = %d ",nphot) ;
    if (nphot < 1) {
       cf_verbose(1, "No photons found in the data.") ;
    } 


    /*  Create output file. */
    FITS_create_file(&outfits, argv[optind], &status);
    /*  Copy primary image header from the input file. */
    FITS_copy_hdu(infits, outfits, morekeys, &status);

    /* Update header keywords. */
    FITS_update_key(outfits, TSTRING, "FILENAME", argv[optind], NULL, &status);
    FITS_update_key(outfits, TSTRING, "FILETYPE",
    		             "INTERMEDIATE DATA FILE", NULL, &status);
    if (modify_hkexists)
	FITS_update_key(outfits, TSTRING, "HKEXISTS", "YES", NULL, &status);

    /* Populate the calibration keywords in the primary header. */
    cf_fuv_init(outfits);

    /* Write background limits to file header */
     cf_set_background_limits(outfits); 

    /* Fill in the first extension with the photon information */
     if (nphot<1) {
       cf_init_null_photon_list(outfits);
     }
     else {
       cf_init_photon_list(infits, outfits);
     }

    /* Copy good-time intervals from input to output. */
    cf_verbose(3, "Copying good-time intervals from input to output.");
    FITS_movabs_hdu(infits, 3, &hdutype, &status);
    FITS_create_hdu(outfits, &status);
    FITS_copy_hdu(infits, outfits, 0, &status);

    /*  Append a table containing the orbital timeline */
    cf_timeline(outfits);
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    /*  Update processing flags. */
    cf_proc_update(outfits, CF_PRGM_ID, "COMPLETE");

    FITS_close_file(infits, &status);
    FITS_close_file(outfits, &status);

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
    return EXIT_SUCCESS;
}
