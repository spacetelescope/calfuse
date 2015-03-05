/**************************************************************************
 *           Johns Hopkins University  
 *          Center for Astrophysical Sciences     
 *           FUSE
 *************************************************************************
 *
 * synopsis:    cf_jitter input_file output_file     
 *
 * Description: Reads a FITS file containing housekeeping data and
 *              generates a FITS file containing the jitter information 
 *              for the exposure.
 *
 * Arguments: input_file    Input raw FITS housekeeping file for an exposure
 *                          The header of this file is used as a template
 *                          for the output file
 *            output_file   FITS table containing the jitter information 
 *     
 * 
 * Contents of output file: 
 *   1. Time of sample (in seconds from start)
 *   2. Shifts dx and dy (in arcsec) of the pointing from a reference position
 *   3. Tracking quality flag
 *	5 = dx, dy from known FPDs and q_cmd from telemetry
 *	4 = good dx, dy from ACS q_est and q_cmd from telemetry
 *	3 = maybe good dx, dy from ACS q_est
 *	2 = dx, dy from known FPDs, but q_cmd from FITS header coordinates
 *	1 = dx, dy from ACS q_est, but q_cmd from FITS header coordinates
 *	0 = no pointing information (missing telemetry)
 *     -1 = Pointing is assumed to be bad (never achieved known track)
 *
 * Codes to JIT_STAT keyword:
 *	0	All OK
 *	1	Problems with jitter data
 *	2	AT_CMD_ATT flag never set
 *	3	Commanded ACS quaternions not available
 *	4	No star field information
 *	5	Problems with tabulated start time
 *
 *                                 
 * HISTORY:
 *       01/22/02   v1.0     Robinson  Begin work
 *       08/14/02   v1.1     Robinson  Populate JIT_STATUS keyword
 *       09/20/02   v1.11    Robinson  - update FILENAME, FILETYPE keywords
 *                                     - change EXP_DUR keyword to reflect
 *                                       total time duration of jitter data
 *       11/04/02   v1.2     Robinson  Change method in which gaps are filled.
 *                                     We no longer assume that the data is 
 *                                     tabulated at regular intervals. Also 
 *                                     check that tabulated times start BEFORE
 *                                     the exposure start time
 *       12/17/02   v1.3     Robinson  Slight change in definitions of TRAKFLG
 *                                     values.
 *       01/17/03   v1.4     Robinson  Slight change in calculating summary of
 *                                     jitter properties
 *       05/22/03   v1.5     rdr       Change initialization of arrays
 *       06/10/03   v1.6     rdr       Correct problem in determining ref 
 *                                     quaternion
 *       06/23/03   v1.7     rdr       Adjust calculation on some of the
 *                                     keywords
 *       06/25/03   v1.8     rdr       Correct fault in setting the SLEWFLG
 *                                     keyword
 *       06/27/03   v1.9     rdr       Corrected error in selecting the 
 *					reference quaternion
 *       04/05/05   v2.0     wvd       Modify logic for setting TRAKFLG.
 *					Use only FPD quaternions for dx & dy.
 *       06/24/05   v2.1     wvd       De-bug code for use in Linux boxes.
 *       07/06/05   v2.2     wvd       Add JIT_VERS keyword to jitter file.
 *       11/30/05   v2.3     wvd       Don't add 0.5 to tval and time[i] 
 *					before converting to type long.
 *       06/28/06   v2.4     wvd       Add "JITTER KEYWORDS" to file header.
 *       08/21/06   v2.5     wvd       Always calculate dx and dy.
 *					For data taken in one-wheel mode,
 *					jitter is bad if sfknown < 1.
 *       08/25/06   v2.6     wvd       Three-wheel mode was BEFORE DEC2001,
 *					not after.  Use qvalid flag to check
 *					validity of FPD quaternions.
 *       08/25/06   v2.7     wvd       New tracking flag values:
 *					5 = FPD, 4 = ACS, 0 = no info
 *       11/30/06   v3.0     wvd       New code from Tom Ake:
 *					If telemetry unavailable, use
 *					target coordinates to generate
 *					commanded quaternion.
 *					Modify logic of tracking flags.
 *       12/19/06   v3.1     wvd       Add flag to force use of target
 *					coordinates.
 *       03/29/07   v3.2     wvd       Modify to handle time period when
 *					ACS telemetry rate was cut in half. 
 *					Make code generally more robust.
 *       03/30/07   v3.2.1   wvd       Correct typos in last change.
 *       04/07/07   v3.3     wvd       Add #include <unistd.h>
 *       04/13/07   v3.4     wvd       If program fails, set NEXTEND = 0
 *					in jitter file header.
 *       05/04/07   v3.5     wvd       Force fourth component of q_meas
 *					to be positive.
 *	 08/24/07   v3.6     bot       changed TINT to TLONG for time l.1191.
 *	 05/30/08   v3.7     wvd       If computing coordinates from file
 *					header and aperture = RFPT, do not
 *					add an aperture offset.
 *
 **************************************************************************/

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calfuse.h"

#define HSKPVERS_MIN    3.1

static char CF_PRGM_ID[] = "cf_jitter";    
static char CF_VER_NUM[] = "3.7";  


/* Subroutine to normalize a quaternion */

void qnorm(double *quat) {

  double quat3_sq;

  quat3_sq = 1. - quat[0]*quat[0] - quat[1]*quat[1] - quat[2]*quat[2];
  if (quat3_sq > 0) quat[3] = sqrt(quat3_sq);
  else quat[3] = 0.;

}


/* Subroutine to perform matrix multiplication of two quaternions */

void quatprod(double *tquat, double *quat_ref, double *dquat) {

  int i, j ;
  double multmat[4][4] ; 

  /* fill in the array with values from the tquat vector */
  multmat[0][0] =  tquat[3] ;
  multmat[0][1] = -tquat[2] ;
  multmat[0][2] =  tquat[1] ;
  multmat[0][3] = -tquat[0] ;
  multmat[1][0] =  tquat[2] ;
  multmat[1][1] =  tquat[3] ;
  multmat[1][2] = -tquat[0] ;
  multmat[1][3] = -tquat[1] ;
  multmat[2][0] = -tquat[1] ; 
  multmat[2][1] =  tquat[0] ;
  multmat[2][2] =  tquat[3] ;
  multmat[2][3] = -tquat[2] ;
  multmat[3][0] =  tquat[0] ;
  multmat[3][1] =  tquat[1] ;
  multmat[3][2] =  tquat[2] ;
  multmat[3][3] =  tquat[3] ;
  
  /* do the matrix multiplication */
  for (i=0 ; i<4 ; i++) {
     dquat[i] = 0. ;
     for (j=0; j<4 ; j++) dquat[i] = dquat[i] + multmat[j][i] * quat_ref[j] ; 
  }

  if (dquat[3] < 0) for (i=0; i<4 ; i++) dquat[i] = -dquat[i];

}


int main(int argc, char *argv[])  
{
  fitsfile *infits, *outfits ;
  int hdutype, tfields, dstatus ;
  int nvalid, nvalid5, rtime ;
  int pflg=0 ;
  int tflg, qvalidt ;
  int tcolnum, fpdquat1_colnum,  fpdquat2_colnum,  fpdquat3_colnum ;
  int aquat1_colnum, aquat2_colnum, aquat3_colnum, slew_colnum ;
  int status, tcol, nextend=0;
  int dxcol, dycol, trkcol, anynull, intnull=0 ;
  int qvalid_colnum, cents1_colnum, cents2_colnum, cents3_colnum;
  int cents4_colnum, cents5_colnum, cents6_colnum, acaflg_colnum ;
  int acscmd1_colnum, acscmd2_colnum, acscmd3_colnum, sfknown_colnum ;
  int frow, felem, slewflg, ngs_used, jitlgx, jitlgy ;
  int *cents1, *cents2, *cents3, *cents4, *cents5, *cents6 , *qvalid ;
  int *sfknown, *slew;
  int cents1t=0, cents2t=0, cents3t=0, cents4t=0, cents5t=0, cents6t=0 ;
  int sfknownt=0 , slewt=0, slewtime=0 ;
  long i, numq5, numq4, numq3, numq2, numq1, numq0, numqa ; 
  long  npts, tval, n_invert=0;
  long *time, qtime=0 ;
  float tfine, tcoarse, tnull, knowntrk, acstrk, fjitlgx, fjitlgy ;
  float targtrk ;
  float ngs1, ngs2, ngs3, ngs4, ngs5, ngs6, gs_used, expdur, exptime ;
  float dxrms, dyrms, dxave, dyave, dxrms5, dyrms5 ;
  float *acaflg ;
  float *dx , *dy ;
  double *fq1, *fq2,*fq3, *mjd ;
  double *acscmd1, *acscmd2, *acscmd3, fq1t, fq2t, fq3t  ;
  double tstrt ;
  double tquat[4], dquat[4], quat_ref[4] ;
  double *aq1, *aq2, *aq3, aq1t, aq2t, aq3t ;
  short *trakflg ;
  char buffer[FLEN_CARD];
  char *ttype[] = { "TIME", "DX", "DY", "TRKFLG" } ;
  char *tform[] = { "1J", "1E", "1E", "1I" } ;
  char *tunit[] = { "seconds", "arcsec", "arcsec", " "} ;
  char textname[7]="jitter" ;
  char blank[]="                                                                                ";
  char jitter_keywords[]="              JITTER KEYWORDS                                                   ";

  /* New variables (11/30/06) */
  int maxi, qcmd_prob, trakflg_init, *aqsok, aq_fill ;
  float targ_ra, targ_dec, targ_roll, x_off, y_off, fpalif1x, fpalif2x ;
  double cra, sra, cdc, sdc, crl, srl, rad2deg ;
  double dq_sid[4], aq_old[4], rot[9] ;
  char aper_id[FLEN_VALUE], fescent[FLEN_VALUE], hskpvers[FLEN_VALUE];

  /* New variables (12/19/06) */
  int optc;
  int force_target_coords = FALSE;

  char opts[] = "hfv:";
  char usage[] =
    "Usage:\n"
    "  cf_jitter [-hf] [-v level] housekeeping_file jitter_file\n";
  char option[] =
    "Options:\n"
    "  -h:  this help message\n"
    "  -f:  force use of target coordinates \n"
    "  -v:  verbosity level (=1; 0 is silent)\n";

  verbose_level = 1;

  /* Check number of options and arguments */
  while ((optc = getopt(argc, argv, opts)) != -1) {
    switch(optc) {
    case 'h':
      printf("%s\n%s", usage, option);
      return 0;
    case 'f':
      force_target_coords=TRUE;
      break;
    case 'v':
      verbose_level = atoi(optarg);
      break;
    }
  }

  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

  if (argc < optind+2) {
    printf("%s", usage);
    return 1;
  }


  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");


  /* initialize necessary variables */
     hdutype=0 ;
     status=0 ;
     frow=1 ;
     felem=1 ;
     anynull=0 ;

  /* Open the input FITS file */
     cf_verbose(3, "opening the file %s for input",argv[optind]);
     FITS_open_file(&infits, argv[optind], READONLY, &status); 

  /* Read the starting time and duration from the FITS header */
     FITS_read_key(infits, TDOUBLE, "EXPSTART", &tstrt, buffer, &status) ;
     cf_verbose(2, "start time (MJD) = %15.9f",tstrt) ;
     FITS_read_key(infits, TFLOAT, "EXPTIME", &exptime, buffer, &status) ;
     cf_verbose(2, "expected exposure duration (sec) = %.0f",exptime) ;
 
  /* Create the output file */
     FITS_create_file(&outfits, argv[optind+1], &status); 

  /* Copy the header and empty primary image of the input file to 
     the output */
     FITS_copy_hdu(infits, outfits, 0, &status) ;

  /* Add "JITTER KEYWORDS" string to header. */
     FITS_write_record(outfits, blank, &status);
     FITS_write_record(outfits, jitter_keywords, &status);
     FITS_write_record(outfits, blank, &status);

  /* Populate the STATUS keyword with 1.
     Reset keyword if everything goes well. */
     dstatus = 1 ;
     FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
         "Problems with jitter data", &status) ;

  /* Update FILENAME and FILETYPE keywords */
     FITS_update_key(outfits, TSTRING, "FILENAME",argv[2], NULL, &status) ;
     FITS_update_key(outfits, TSTRING, "FILETYPE","JITTER", NULL, &status) ;

  /* Insert JIT_VERS keyword */
     FITS_update_key(outfits, TSTRING, "JIT_VERS", CF_VER_NUM,
	"cf_jitter version number", &status) ;

  /* Read coordinate information in case we can't determine quat_ref from telemetry */
     /* Don't bother reading these keywords if the HSKPVERS is too low. */
     fits_read_key(infits, TSTRING, "HSKPVERS", hskpvers, NULL, &status);
     if (status != 0 || atof(hskpvers) < HSKPVERS_MIN) {
        status = 0;
        cf_verbose(1, "Housekeeping file format is out of date.");
     }
     else {
	FITS_read_key(infits, TFLOAT, "RA_TARG", &targ_ra, buffer, &status) ;
	FITS_read_key(infits, TFLOAT, "DEC_TARG", &targ_dec, buffer, &status) ;
	FITS_read_key(infits, TFLOAT, "APER_PA", &targ_roll, buffer, &status) ;
	FITS_read_key(infits, TSTRING, "APERTURE", &aper_id, buffer, &status) ;
	FITS_read_key(infits, TFLOAT, "FPALIF1X", &fpalif1x, buffer, &status) ;
	FITS_read_key(infits, TFLOAT, "FPALIF2X", &fpalif2x, buffer, &status) ;
	FITS_read_key(infits, TSTRING, "FESCENT", &fescent, buffer, &status) ;
     }

  /* move to the first extension, which contains the table */
     FITS_movabs_hdu(infits, 2, &hdutype, &status) ;

  /* determine the number of times which have been tabulated */
     FITS_read_key(infits, TLONG, "NAXIS2", &npts, NULL, &status) ;
     cf_verbose(2, "number of samples in the table = %ld", npts) ;

  /*   set up arrays to contain the data to be read from the table  */
           mjd = (double *) cf_calloc(npts, sizeof(double)) ;
           qvalid = (int *) cf_calloc(npts, sizeof(int)) ;
           fq1 = (double *) cf_calloc(npts, sizeof(double) ) ;
           fq2 = (double *) cf_calloc(npts, sizeof(double) ) ;
           fq3 = (double *) cf_calloc(npts, sizeof(double) ) ;
           aq1 = (double *) cf_calloc(npts, sizeof(double) ) ;
           aq2 = (double *) cf_calloc(npts, sizeof(double) ) ;
           aq3 = (double *) cf_calloc(npts, sizeof(double)) ; 
           acaflg = (float *) cf_calloc(npts, sizeof(float)) ;
           acscmd1 = (double *) cf_calloc(npts, sizeof(double)) ; 
           acscmd2 = (double *) cf_calloc(npts, sizeof(double) ) ; 
           acscmd3 = (double *) cf_calloc(npts, sizeof(double) ) ; 
           cents1 = (int *) cf_calloc(npts, sizeof(int)) ;
           cents2 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents3 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents4 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents5 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents6 = (int *) cf_calloc(npts, sizeof(int) ) ;
           sfknown = (int *) cf_calloc(npts, sizeof(int) ) ;
           slew = (int *) cf_calloc(npts, sizeof(int) ) ;

  /* read the data */
	   FITS_get_colnum(infits, CASEINSEN, "MJD", &tcolnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, tcolnum, frow, felem, npts, &intnull,
			   mjd, &anynull, &status) ;

	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_MEAS_Q_VALID",
		&qvalid_colnum, &status) ;
	   FITS_read_col(infits, TINT, qvalid_colnum, frow, felem, npts,
		&intnull, qvalid, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_Q_ECI2BDY_1",
		&fpdquat1_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, fpdquat1_colnum , frow, felem, npts,
		&intnull, fq1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_Q_ECI2BDY_2",
		&fpdquat2_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, fpdquat2_colnum , frow, felem, npts,
		&intnull, fq2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_Q_ECI2BDY_3",
		&fpdquat3_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, fpdquat3_colnum , frow, felem, npts,
		&intnull, fq3, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AATTQECI2BDY_1",
		&aquat1_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, aquat1_colnum , frow, felem, npts,
		&intnull, aq1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AATTQECI2BDY_2",
		&aquat2_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, aquat2_colnum, frow, felem, npts,
		&intnull, aq2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AATTQECI2BDY_3",
		&aquat3_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, aquat3_colnum, frow, felem, npts,
		&intnull, aq3, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_AT_CMD_ATT", &acaflg_colnum,
		&status) ;
	   FITS_read_col(infits, TFLOAT, acaflg_colnum, frow, felem, npts,
		&intnull, acaflg, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AQECI2BDYCMD_1",
		&acscmd1_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, acscmd1_colnum, frow, felem, npts,
		&intnull, acscmd1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AQECI2BDYCMD_2",
		&acscmd2_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, acscmd2_colnum, frow, felem, npts,
		&intnull, acscmd2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AQECI2BDYCMD_3",
		&acscmd3_colnum, &status) ;
	   FITS_read_col(infits, TDOUBLE, acscmd3_colnum, frow, felem, npts,
		&intnull, acscmd3, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID1_STATUS",
		&cents1_colnum, &status) ;
	   FITS_read_col(infits, TINT, cents1_colnum, frow, felem, npts,
		&intnull, cents1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID2_STATUS",
		&cents2_colnum, &status) ;
	   FITS_read_col(infits, TINT, cents2_colnum, frow, felem, npts,
		&intnull, cents2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID3_STATUS",
		&cents3_colnum, &status) ;
	   FITS_read_col(infits, TINT, cents3_colnum, frow, felem, npts,
		&intnull, cents3, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID4_STATUS",
		&cents4_colnum, &status) ;
	   FITS_read_col(infits, TINT, cents4_colnum, frow, felem, npts,
		&intnull, cents4, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID5_STATUS",
		&cents5_colnum, &status) ;
	   FITS_read_col(infits, TINT, cents5_colnum, frow, felem, npts,
		&intnull, cents5, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID6_STATUS",
		&cents6_colnum, &status) ;
	   FITS_read_col(infits, TINT, cents6_colnum, frow, felem, npts,
		&intnull, cents6, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_STARFLD_KNWN",
		&sfknown_colnum, &status) ;
	   FITS_read_col(infits, TINT, sfknown_colnum, frow, felem, npts,
		&intnull, sfknown, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_NEW_CMD_Q", &slew_colnum,
                &status) ;
	   FITS_read_col(infits, TINT, slew_colnum, frow, felem, npts,
		&intnull, slew, &anynull, &status) ;

  /* close the input file */
      FITS_close_file(infits, &status) ;

  /* Check that the tabulated times start BEFORE the start of the
	exposure.  Otherwise there may be missing data or an error
	in the tabulated start time. */
     tval=cf_nlong(floor((mjd[0]-tstrt) * 86400.));
     if (tval > 0) {
        cf_if_warning("Possible missing data or error in the start time.  "
        "The first tabulated time is %ld seconds from obs start.", tval) ;
	/* update JIT_STAT keyword */
        dstatus = 5 ;
        FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
            "Problems with tab start time", &status) ;
     }      

  /* fill in the gaps in the data set 

     - ACS quaternions are generally tabulated every other second
     - FPD quaternions can be every second or every other second
     - ACS flags (acaflg) - tabulated every 16 seconds

     NOTE:  there can be irregularities in the spacing.  Thus, save the last
     tabulated value and use this to fill in missing data for a given time. */

  /* Array to track whether we can trust ACS q_meas for dx,dy computation */
     aqsok = (int *) cf_calloc(npts, sizeof(int) ) ;

  /* ACS quaternions. Fill up to aq_fill missing values with the last good one.
     Set aqsok[i] to -1 if we leave a missing value missing.*/

     /* aq_fill = 2 ; replace up to 2 consecutive missing values */

     /* Fix for 2007:072:18:02 - 085:13:38, when ACS telemetry rate
	was cut in half (MOCR 1902). */
     aq_fill = 4;
     aq1t = 0. ;
     aq2t = 0. ;
     aq3t = 0. ;
     tflg = -1 ;

     /* If the first ACS aq_fill values are missing,
	use the one aq_fill steps later. */
     for (i = aq_fill; i > 0; i--) if (aq1[i] > -1 && aq2[i] > -1) {
	aq1t = aq1[i];
	aq2t = aq2[i];
	aq3t = aq3[i];
	tflg = aq_fill;
     }

     for (i=0 ; i<npts ; i++) {
       aqsok[i] = -1 ; /* assume q_meas is missing */
       if ( aq1[i] > -1 && aq2[i] > -1 ) {
         aqsok[i] = 0 ;
         aq1t = aq1[i] ;
         aq2t = aq2[i] ;
         aq3t = aq3[i] ;
         tflg = aq_fill ; /* fill factor */
       } 
       else if (tflg > 0) {
         aqsok[i] = 0 ;
         aq1[i] = aq1t ;
         aq2[i] = aq2t ;
         aq3[i] = aq3t ;
         tflg-- ;
       }
     }    

      /* FPD quaternions and qvalid flags 
         Only update the one following a good quaternion. If there is
	 another one with no quaternion after that, do not update it. */

     fq1t = 0. ;
     fq2t = 0. ;
     fq3t = 0. ;
     qvalidt = 0 ;
     tflg = -1 ;
      for (i=0 ; i<npts ; i++) {
        if ( (fabs(fq1[i]) > 0.) && (fabs(fq2[i]) > 0.) ){

   /* The ACS has a strict convention that the 4th component of the
   quaternion be positive. But if it is close to zero, the IDS FPD
   component may be negative. The housekeeping file only has components
   1-3. As a best estimate, we will choose the sign of the FPD 4th
   component so that the quaternion is closest in pointing to the ACS
   quaternion. Only known tracking FPDs need to be corrected. */


if (aq1[i] > -1 && aq2[i] > -1 && sfknown[i] == 1) {
  quat_ref[0] = -aq1[i] ;
  quat_ref[1] = -aq2[i] ;
  quat_ref[2] = -aq3[i] ;
  qnorm(quat_ref) ;
  if (quat_ref[3] < 0.03) {
    tquat[0] = fq1[i] ;
    tquat[1] = fq2[i] ;
    tquat[2] = fq3[i] ;
    qnorm(tquat) ;
    quatprod(tquat,quat_ref,dquat) ;
    fq1t = dquat[0]*dquat[0] + dquat[1]*dquat[1] ;
    tquat[3] = -tquat[3] ;
    quatprod(tquat,quat_ref,dquat) ;
    if (dquat[0]*dquat[0]+dquat[1]*dquat[1] < fq1t) {
      fq1[i] = -fq1[i] ;
      fq2[i] = -fq2[i] ;
      fq3[i] = -fq3[i] ; 
      n_invert++;
    }
  }
}

	  fq1t=fq1[i] ;
	  fq2t=fq2[i] ;
	  fq3t=fq3[i] ; 
          qvalidt = qvalid[i] ;
          cents1t = cents1[i] ;
          cents2t = cents2[i] ;
          cents3t = cents3[i] ;
          cents4t = cents4[i] ;
          cents5t = cents5[i] ;
          cents6t = cents6[i] ;
          sfknownt = sfknown[i] ;
          slewt = slew[i] ;
          tflg = 1 ; }
        else if (tflg > 0) {
          fq1[i] = fq1t ;
          fq2[i] = fq2t ;
          fq3[i] = fq3t ;
          qvalid[i] = qvalidt ;
          cents1[i] = cents1t ; 
          cents2[i] = cents2t ; 
          cents3[i] = cents3t ; 
          cents4[i] = cents4t ; 
          cents5[i] = cents5t ; 
          cents6[i] = cents6t ; 
          sfknown[i] = sfknownt ;
          slew[i] = slewt ;
          tflg = -1 ; }
	}
if (n_invert > 0) cf_verbose(1, "Number of known FPD quaternions inverted: %ld", n_invert);


       /* ACAFLG flag */

       tflg = -10 ;          

       /* Assume that the first real ACA value is the one to start with. */
	for (i = 0; i < npts-1; i++) if (acaflg[i] > -1) {
	    tflg = acaflg[i];
	    break;
	}

       for (i=0 ; i<npts-1 ; i++) {
         if (acaflg[i] > -1) tflg = acaflg[i] ;
         else acaflg[i] = tflg ;
        } 
		       

 /********** derive the final quats arrays and the tracking flags *******/

  /*    set up arrays to contain the data */
         time = (long *) cf_calloc(npts, sizeof(long) ) ;
	 trakflg = (short *) cf_calloc(npts, sizeof(short) ) ; 
         dx = (float *) cf_calloc(npts, sizeof(float) ) ;
         dy = (float *) cf_calloc(npts, sizeof(float) ) ;

 /* fill in the time array with times (in seconds) from the start
	of the exposure */
	 for (i=0 ; i<npts ; i++) 
           time[i] = cf_nlong(floor((mjd[i]-tstrt) * 86400.));


   /* Determine the reference quat by looking at the commanded ACS
   quaternion (in arrays acscmd1, acscmd2 and acscmd3) when they have
   converged to the fpd commanded quaternion - i.e. when the AT_CMD_ATT
   flag is 1 (held in the array acaflg).  Make sure that the time of
   convergence is AFTER the start of the exposure and that the star field
   is known (i.e. sfknown = 1).  Also, because of the 16s sample time of
   the ascflg array, require that the condition exists for at least 30
   consecutive seconds before accepting the quaternions.  qtime = time at
   which the quaternions were found. This is the earliest time in the
   exposure for which reliable positions are available. */

    rtime=0 ;
    qtime = 0;
    for (i=0; i<npts; i++) {
	if (acaflg[i] == 1 && time[i] > 0 && sfknown[i] == 1) {
	    rtime++ ;
	    if (rtime == 1) qtime=time[i] ;
	    if (rtime > 30 && acscmd1[i] > -1 &&
	    	acscmd2[i] > -1 && acscmd3[i] > -1 ) {
		quat_ref[0]=acscmd1[i] ;
		quat_ref[1]=acscmd2[i] ;
		quat_ref[2]=acscmd3[i] ;
		qnorm(quat_ref);
		qcmd_prob = 0;
		break;
	    }
	}
	else rtime=0 ; 
    }

     /* If rtime <= 30, then the flag was never set
	and the reference quaternion not determined */
    if (rtime <= 30 || force_target_coords) {

	/* If user has requested target coordinates... */
	if (force_target_coords) {
	    cf_verbose(1, "Computing reference quaternion from target coordinates.");
            dstatus = 2 ;
            FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
               "User requests use of target coordinates", &status) ;
	}
	/* If something is wrong... */
	else {
	 cf_verbose(1, "Problem determining the reference quaternions") ;
         pflg = 0 ;
	 for (i=0 ; i<npts; i++) if (acaflg[i] == 1) pflg=1 ;
	 if (pflg == 0) {
               cf_verbose(1, " AT_CMD_ATT flag was never set.") ;
         	/* update JIT_STAT keyword */
               dstatus = 2 ;
               FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
               "AT_CMD_ATT flag never set", &status) ;
	 }

         pflg = 0 ;
	 for (i=0 ; i<npts; i++) if (acscmd1[i] > -1 &&
                acscmd2[i] > -1 && acscmd3[i] > -1) pflg++ ;
	 if (pflg <= 30) {
               cf_verbose(1, " Commanded ACS quaternions not available.") ;
         	/* update JIT_STAT keyword */
               dstatus = 3 ;
               FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
               "Commanded ACS quaternions not available", &status) ;
	 }

         pflg = 0 ;
	 for (i=0 ; i<npts; i++) if (sfknown[i] == 1) pflg++;
	    if (pflg <= 30)  {
                cf_verbose(1, " NO STAR FIELD INFORMATION") ;
   	        /* update JIT_STAT keyword */
                dstatus = 4 ;
                FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
               "No star field information", &status) ;
	    }
	}

     /* If we lack a reference quaternion and can't calculate one from
	information in the housekeeping file header, set the JIT_STAT
	keyword to 1, close the jitter file, and return. */
     if (atof(hskpvers) < HSKPVERS_MIN) {
	dstatus = 1 ;
	FITS_update_key(outfits, TINT, "NEXTEND",  &nextend, NULL, &status);
        FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, NULL, &status);
	FITS_close_file(outfits, &status);
	cf_verbose(1, "Unable to calculate reference quaternion from target coordinates.");
	return 1;
     }

     /* Compute quat_ref from the target coordinates instead */
     qcmd_prob = -3 ;
     rad2deg = 57.2957951 ;
     cra = cos(targ_ra/rad2deg) ;
     sra = sin(targ_ra/rad2deg) ;
     cdc = cos(targ_dec/rad2deg) ;
     sdc = sin(targ_dec/rad2deg) ;
     targ_roll=270.-targ_roll ;
     crl = cos(targ_roll/rad2deg) ;
     srl = sin(targ_roll/rad2deg) ;
     
     /* The FUSE rotation matrix */
     rot[0] = -srl*sra + crl*sdc*cra ; 
     rot[1] = -crl*sra - srl*sdc*cra ; 
     rot[2] =  cdc*cra ;
     rot[3] =  srl*cra + crl*sdc*sra ;
     rot[4] =  crl*cra - srl*sdc*sra ;
     rot[5] =  cdc*sra ;
     rot[6] = -crl*cdc ;
     rot[7] =  srl*cdc ;
     rot[8] =  sdc ;
     
     /*  quat_ref for the target, before updating for aperture and FPA position */
     
     quat_ref[0] = 0.5 * sqrt(1. + rot[0] + rot[4] + rot[8]) ;
     quat_ref[1] = 0.5 * sqrt(1. + rot[0] - rot[4] - rot[8]) ;
     quat_ref[2] = 0.5 * sqrt(1. + rot[4] - rot[0] - rot[8]) ;
     quat_ref[3] = 0.5 * sqrt(1. + rot[8] - rot[0] - rot[4]) ;
     
     maxi=0 ; /* going to choose the largest divisor */
     for (i=1 ; i<4 ; i++) if (quat_ref[i] > quat_ref[maxi]) maxi=i;
     switch (maxi) {
       case 0:  quat_ref[3] = quat_ref[0] ;
                quat_ref[0] = (rot[7] - rot[5])/quat_ref[3]/4. ;
                quat_ref[1] = (rot[2] - rot[6])/quat_ref[3]/4. ;
                quat_ref[2] = (rot[3] - rot[1])/quat_ref[3]/4. ;
                break ;
       case 1:  quat_ref[0] = quat_ref[1] ;
                quat_ref[3] = (rot[7] - rot[5])/quat_ref[0]/4. ;
                if (quat_ref[3] < 0) {
                  quat_ref[0] = -quat_ref[0] ;
                  quat_ref[3] = -quat_ref[3] ;
                } 
                quat_ref[1] = (rot[1] + rot[3])/quat_ref[0]/4. ;
                quat_ref[2] = (rot[2] + rot[6])/quat_ref[0]/4. ;
                 break ; 
       case 2:  quat_ref[1] = quat_ref[2] ;
                quat_ref[3] = (rot[2] - rot[6])/quat_ref[1]/4. ;
                if (quat_ref[3] < 0) {
                  quat_ref[1] = -quat_ref[1] ; 
                  quat_ref[3] = -quat_ref[3] ;
                }
                quat_ref[0] = (rot[1] + rot[3])/quat_ref[1]/4. ;
                quat_ref[2] = (rot[5] + rot[7])/quat_ref[1]/4. ;
                break ;
       default: quat_ref[2] = quat_ref[3] ;
                quat_ref[3]=(rot[3] - rot[1])/quat_ref[2]/4. ;
                if (quat_ref[3] < 0) {
                  quat_ref[2] = -quat_ref[2] ;
                  quat_ref[3] = -quat_ref[3] ;
                }
                quat_ref[0] = (rot[6] + rot[2])/quat_ref[2]/4. ;
                quat_ref[1] = (rot[5] + rot[7])/quat_ref[2]/4. ;
     } 
     
     /* Compute x_offset due to FPA shift, first in microns, then in arcsec */
     if (!strncmp(fescent, "FES A", 5)) x_off = 117. - fpalif1x;
     else x_off = 175. - fpalif2x ;
     x_off = x_off/10.88*cos(30.5/rad2deg) ;
     y_off = 0. ;
     
     /* Now add in the aperture offset (in arcsec) */
     if (strncmp(aper_id, "RFPT", 4)) {		/* If not RFPT, assume LWRS. */
        x_off = x_off + 55.18 ;
        y_off = 118.07 ;
        if (!strncmp(aper_id, "MDRS", 4)) y_off = -90.18;
        else if (!strncmp(aper_id, "HIRS", 4)) y_off = 10.27;
        if (tstrt < 51708.9167) {
          y_off = 112.07 ;
          if (!strncmp(aper_id, "MDRS", 4)) y_off = -97.18;
          else if (!strncmp(aper_id, "HIRS", 4)) y_off = 2.27;
        } else if (tstrt < 51796.3000) {
          y_off = 116.07 ;
          if (!strncmp(aper_id, "MDRS", 4)) y_off = -93.18;
          else if (!strncmp(aper_id, "HIRS", 4)) y_off = 6.27;
        }
     }
     
     /* Set up a slew for the offsets and update quat_ref to be the expected q_cmd */
     dquat[0] = -y_off/2./206264.81; 
     dquat[1] =  x_off/2./206264.81; 
     dquat[2] = 0.; 
     qnorm(dquat) ;
     quatprod(dquat,quat_ref,tquat); 
     for (i=0 ; i<4 ; i++) quat_ref[i] = tquat[i]; 

    }

    if (qcmd_prob < 0) {
       FITS_update_key(outfits,TSTRING,"QREF_SRC","COORDINATES","Source of reference quaternion", &status) ;
       cf_verbose(1, "Reference quat values derived from FITS header coordinates.");
       cf_verbose(1, "Used x_off = %7.2f, y_off = %7.2f arcseconds.", x_off, y_off);
    }
    else {
       FITS_update_key(outfits,TSTRING,"QREF_SRC","TELEMETRY","Source of reference quaternion", &status) ;
       cf_verbose(1, "Reference quat values represent ACS commanded position.");
    }

    for (i=0; i<4 ; i++) cf_verbose(2, "for i=%ld, quat=%18.15f",i,quat_ref[i]) ;              
    FITS_update_key(outfits,TDOUBLE,"QREF_0",quat_ref,"Reference quaternion", &status) ;
    FITS_update_key(outfits,TDOUBLE,"QREF_1",quat_ref+1, NULL, &status) ;
    FITS_update_key(outfits,TDOUBLE,"QREF_2",quat_ref+2, NULL, &status) ;
    FITS_update_key(outfits,TDOUBLE,"QREF_3",quat_ref+3, NULL, &status) ;

      /* Determine the shift in the quaternions and derive the changes
	 in x and y that they represent:
         DX (from q2 values : q2=dx/2) and DY (from q1 values : q1=-dy/2)
         in arcseconds - the conversion factor used is the number of arcsec
         per radian * 2 */

/* Determine initial trakflg. If we had IDS telemetry but never reached known
track, we want to start with trakflg_init=-1. If telemetry missing, temporarily
set trakflg_init to -2. */

trakflg_init=0 ; tflg=0 ;
for (i=0 ; i<npts ; i++){
  if (sfknown[i] == 1) {
    tflg=1 ;
    break ; 
  }
} 
if (tflg == 0) { /* known track flag never set */
  trakflg_init = -1 ;
  for (i=0 ; i<npts ; i++) {
    if (acaflg[i] > -1) {  /* use ACA flag to say if IDS tlm missing */
      tflg = 1 ;
      break ;
    }
  }  
  if (tflg == 0) trakflg_init = -2;   /* ktrk flag never=1 because IDS tlm missing */
}

/* set up aqsok array to flag whether ACS q_meas is a good estimate of pointing
 -1 = ACS q_meas missing (from filling in the gaps earlier)
  0 = don't trust it (e.g., in unknown track)
  1 = trust it (had been in known track, so good for a while)
  2 = maybe ok (started expo with IDS idled (or missing FPDs) and ACA=1)
  3 = we corrected ACS q_meas during unknown track period from the 1st known q 
*/

/* First, if we got to known track at some point, correct ACS q's during unknown
track that eventually got to known track (aqsok=3) by going backward in time.
Look for transition from known to unknown and use the delta in ACS q's as an
estimate of the dq from starID (dq_sid). Apply this to ACS q's while in unknown
track. Because we filled missing ACS q's with previous ones, we have to skip by
aq_fill values to make sure we have new q's when calculating dq_sid. */

if (trakflg_init > -1) {
  tflg = 0  ; /* used to indicate we have a dq_sid to apply */
  for (i=npts-1 ; i>aq_fill ; i--) {
    if (sfknown[i] == 1 && sfknown[i-aq_fill-1] == 1) {
      tflg = 0 ; /* don't use previous dq_sid anymore */
    }
    else if (sfknown[i] == 1 && 	/* Now in known track */
	(sfknown[i-aq_fill-1] == 0 && qvalid[i-aq_fill-1] == 1) &&  /* Were in unknown track */
	tflg == 0 &&		/* We need to compute a quaternion. */
        aqsok[i] > -1 && aqsok[i-aq_fill-1] > -1) {    /* calculate dq_sid */
      tquat[0] = fq1[i] ;
      tquat[1] = fq2[i] ;
      tquat[2] = fq3[i] ;
      qnorm(tquat) ;
      aq_old[0] = -aq1[i-aq_fill-1] ;
      aq_old[1] = -aq2[i-aq_fill-1] ;
      aq_old[2] = -aq3[i-aq_fill-1] ;
      qnorm(aq_old) ;
      quatprod(tquat,aq_old,dq_sid) ;
      tflg=1 ;
    }
    else if (sfknown[i] == 0 && tflg == 1 && aqsok[i] > -1) { /* correct ACS q_meas */
      aq_old[0]=aq1[i] ;
      aq_old[1]=aq2[i] ;
      aq_old[2]=aq3[i] ;
      qnorm(aq_old) ;
      quatprod(dq_sid,aq_old,tquat) ;
      aq1[i]=tquat[0] ;
      aq2[i]=tquat[1] ;
      aq3[i]=tquat[2] ;
      aqsok[i]=3 ; /* ACS q_meas now believable */
    }
    else if (sfknown[i] == -1) {
        tflg = 0 ; /* IDS idled or missing tlm */
    }
  }
}               /* end trakflg > -1 */

/* If IDS telemetry missing, start with trakflg=0, else don't trust ACS q's yet (trkflg = -1). */
if (trakflg_init == -2) trakflg_init = 0; 
else trakflg_init = -1;

/* If we start in corrected unknown track, initialize trakflg to 4 or 1, else trakflg_init = -1. */
if (aqsok[0] == 3) trakflg_init = 4 + qcmd_prob;

/* If we're already in known track at the start, trakflg_init = 5 or 2. */
if (sfknown[0] == 1) trakflg_init = 5 + qcmd_prob;

/* Now assign aqsok 1-2 cases, trakflg, and compute dqs. First we need to 
make quat_ref its inverse */
quat_ref[3] = -quat_ref[3] ;

/* if start with no FPDs and ACA=1, there may have been a previous exposure
that lost track at the end. Set first aqsok to maybe ok. */
if (sfknown[0] == -1 && acaflg[0] == 1 && aqsok[0] > -1) {
	aqsok[0] = 2;
	trakflg_init = 3;
}

trakflg[0] = trakflg_init ;
for (i=1 ; i<npts; i++) {
  trakflg[i] = trakflg_init ;
  if (sfknown[i] == 1) {
    aqsok[i] = 1 ; /* show we can start trusting the ACS q's if we lose stars later */
    trakflg[i] = 5 + qcmd_prob ; /* will be 5 or 2 */
    tquat[0]=fq1[i] ;
    tquat[1]=fq2[i] ; 
    tquat[2]=fq3[i] ;
    qnorm(tquat) ; 
    trakflg_init = 0 ;
  } 
  else if (aqsok[i-1] > 0 || aqsok[i] > 0) { /* 1,2,3 = already had been trusting the ACS q's */
     /* continue to trust them under the following conditions:
        - we haven't done a slew and we started with q's maybe ok or corrected 
        - IDS is idled or we have missing FPDs
     */
     if ( (aqsok[i-1] == 2 && acaflg[i] == 1) || (qvalid[i] < 1) || (aqsok[i] == 3) ) {
	/* continue to trust them */
      if (aqsok[i] != 3) aqsok[i] = aqsok[i-1]; 
      if (aq1[i] == -1 && aq2[i] == -1) trakflg[i] = 0;
      else if (aqsok[i] == 2) trakflg[i] = 3 ;
      else trakflg[i] = 4 + qcmd_prob ;		/* aqsok=1,3 */

      tquat[0]=aq1[i] ;
      tquat[1]=aq2[i] ;
      tquat[2]=aq3[i] ;
      qnorm(tquat) ;
      trakflg_init = 0 ;
    } 
    else if (qvalid[i] == 1) {
	trakflg[i] = -1 ;
	trakflg_init = -1;
    }
  }

  if (trakflg[i] > 0 && tquat[0] > -1 && tquat[1] > -1) { 
    quatprod(tquat,quat_ref,dquat) ;
    dx[i] = (float) dquat[1]*206264.81*2. ;
    dy[i] = (float) -dquat[0]*206264.81*2. ;
  } 
  else {
    dx[i] = 0. ; dy[i] = 0. ;
  }
}


/*Provide keywords in the top-level header which describe the exposure */

   /* determine the time span in the jitter file (in seconds) and set EXP_DUR to
        this value  */
       expdur = (float) (time[npts-1] - time[0]) + 1. ;
       FITS_update_key(outfits, TFLOAT, "EXP_DUR", &expdur, 
         "[s] duration of jitter data", &status) ;

       /* check for a short exposure and reset exptime if necessary */
       if (time[npts-1] < exptime) {
	 exptime = time[npts-1] ;
         cf_verbose(1, "Short exposure. Reset exptime to %.0f",exptime) ;
       }

	 /* determine the statistics of the guiding */
         numq5=0 ;
         numq4=0 ;
         numq3=0 ;
         numq2=0 ;
         numq1=0 ;
         numq0=0 ;
         numqa=0 ;

	 for (i=0; i<npts ; i++) {
	   cf_verbose(4, "at i=%5d, trakflg=%2d, qvalid=%2d, time=%4ld",
	      i,trakflg[i], qvalid[i], time[i] );
           if(trakflg[i] == 5 && qvalid[i] > 0 && time[i] > 0 &&
		time[i] <= exptime) 
                   numq5++ ;
           if(trakflg[i] == 4 && time[i] > 0 && time[i] <= exptime) numq4++ ;
           if(trakflg[i] == 3 && time[i] > 0 && time[i] <= exptime) numq3++ ;
           if(trakflg[i] == 2 && time[i] > 0 && time[i] <= exptime) numq2++ ;
           if(trakflg[i] == 1 && time[i] > 0 && time[i] <= exptime) numq1++ ;
           if(trakflg[i] == 0 && time[i] > 0 && time[i] <= exptime) numq0++ ;
           if(trakflg[i] ==-1 && time[i] > 0 && time[i] <= exptime) numqa++ ; }

         cf_verbose(2, "duration of jitter data = %.0f, exposure time = %.0f",
		expdur, exptime) ;
	 cf_verbose(2, "number of samples with trkflg =  5: %ld", numq5);
	 cf_verbose(2, "number of samples with trkflg =  4: %ld", numq4);
	 cf_verbose(2, "number of samples with trkflg =  3: %ld", numq3);
	 cf_verbose(2, "number of samples with trkflg =  2: %ld", numq2);
	 cf_verbose(2, "number of samples with trkflg =  1: %ld", numq1);
	 cf_verbose(2, "number of samples with trkflg =  0: %ld", numq0);
	 cf_verbose(2, "number of samples with trkflg = -1: %ld", numqa);


       /* tfine, tcoarse and tnull are the fraction of the total exposure
          when the exposure is in fine track, coarse track and no track */
       tfine = (float) (numq2+numq5) /  exptime ;
       tcoarse = (float) (numq1+numq3+numq4) /  exptime ;
       tnull = (float) (numqa+numq0) /  exptime ;
       FITS_update_key(outfits, TFLOAT, "KNOWNTRK", &tfine, 
           "fraction of exp in known track", &status) ;
       FITS_update_key(outfits, TFLOAT, "ACSTRK", &tcoarse, 
           "fraction of exp in ACS track", &status) ;
       FITS_update_key(outfits, TFLOAT, "NOGDEINF", &tnull, 
           "fraction of exp with no guiding info", &status) ;

       /* determine the properties of the guide stars */
       ngs1=0. ;
       ngs2=0. ;
       ngs3=0. ;
       ngs4=0. ;
       ngs5=0. ;
       ngs6=0. ;
       for (i=0; i<npts ; i++) {
         if (cents1[i] == 2 && sfknown[i] == 1 && time[i] > 0 && time[i] <= exptime) ngs1+=1. ;
         if (cents2[i] == 2 && sfknown[i] == 1 && time[i] > 0 && time[i] <= exptime) ngs2+=1. ;
         if (cents3[i] == 2 && sfknown[i] == 1 && time[i] > 0 && time[i] <= exptime) ngs3+=1. ;
         if (cents4[i] == 2 && sfknown[i] == 1 && time[i] > 0 && time[i] <= exptime) ngs4+=1. ;
         if (cents5[i] == 2 && sfknown[i] == 1 && time[i] > 0 && time[i] <= exptime) ngs5+=1. ;
         if (cents6[i] == 2 && sfknown[i] == 1 && time[i] > 0 && time[i] <= exptime) ngs6+=1. ; }

       /* convert to fraction of observing time */
       ngs1 = ngs1 / exptime ;
       ngs2 = ngs2 / exptime ;
       ngs3 = ngs3 / exptime ;
       ngs4 = ngs4 / exptime ;
       ngs5 = ngs5 / exptime ;
       ngs6 = ngs6 / exptime ;

	/* determine the number of guide stars used during the exposure */
	ngs_used = 0 ;
	if (ngs1 > 0) ngs_used += 1 ; 
	if (ngs2 > 0) ngs_used += 1 ; 
	if (ngs3 > 0) ngs_used += 1 ; 
	if (ngs4 > 0) ngs_used += 1 ; 
	if (ngs5 > 0) ngs_used += 1 ; 
	if (ngs6 > 0) ngs_used += 1 ; 
	cf_verbose(2, "maximum number of guide stars used: %d",ngs_used) ;

       /* Populate header keywords */
        FITS_update_key(outfits, TINT, "NGS_USED", &ngs_used, 
           "number of guide stars used", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS1_USED", &ngs1, 
           "fraction of exp using known guide star 1", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS2_USED", &ngs2, 
           "fraction of exp using known guide star 2", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS3_USED", &ngs3, 
           "fraction of exp using known guide star 3", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS4_USED", &ngs4, 
           "fraction of exp using known guide star 4", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS5_USED", &ngs5, 
           "fraction of exp using known guide star 5", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS6_USED", &ngs6, 
           "fraction of exp using known guide star 6", &status) ;

	/* Determine times of known and unknown tracking and
		tracking on the target */ 
	knowntrk = 0. ;
        acstrk = 0. ;
        targtrk = 0. ;
        gs_used = 0. ;
	for (i=0 ; i< npts ; i++) {
          if (qvalid[i] == 1 && time[i] > 0 && time[i] <= exptime) gs_used+=1.;
          if ((trakflg[i] == 5 || trakflg[i] == 2) &&
		fabs(dx[i]) < 15 && fabs(dy[i]) < 15 && time[i] > 0 && time[i] <= exptime) 
                knowntrk += 1. ;
          if ((trakflg[i] == 4 || trakflg[i] == 3 || trakflg[i] == 1) &&
		fabs(dx[i]) < 15 && fabs(dy[i]) < 15 && time[i] > 0 && time[i] <= exptime) 
                acstrk += 1. ; 
          if (trakflg[i] > 0 &&
                fabs(dx[i]) < 15 && fabs(dy[i]) < 15 && time[i] > 0 && time[i] <= exptime)
	        targtrk += 1. ;
     }
        cf_verbose(2, "seconds when guide stars are used: %.0f",gs_used ) ;
        cf_verbose(2, "seconds tracking on known stars with error < 15 arcsec: %.0f",knowntrk) ;
        cf_verbose(2, "seconds tracking on ACS with error < 15 arcsec: %.0f", acstrk) ;
        cf_verbose(2, "seconds tracking on target with error < 15 arcsec: %.0f",targtrk ) ;
	knowntrk = knowntrk /  exptime ;
	acstrk = acstrk / exptime ;
        targtrk = targtrk / exptime ;
        gs_used = gs_used / exptime ;
         FITS_update_key(outfits, TFLOAT, "KTRKFINE", &knowntrk, 
            "exp fraction, err < 15 arcsec, known stars", &status) ;
         FITS_update_key(outfits, TFLOAT, "ATRKFINE", &acstrk, 
            "exp fraction, err < 15 arcsec, ACS", &status) ;
        FITS_update_key(outfits, TFLOAT, "TOTLFINE", &targtrk, 
            "exp fraction, err < 15 arcsec, total", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS_INUSE", &gs_used, 
            "fraction of exp tracking on guide stars", &status) ;


	/* determine whether a new attitude is commanded during the exposure */
        slewflg = -1 ;
	for (i=0; i<npts ; i++) 
          if(slew[i] > 0 && time[i] > 0 && time[i] <= exptime) {
		slewflg = 1 ;
		break;
	}
        FITS_update_key(outfits, TINT, "SLEWFLG", &slewflg, 
            "Slew commanded during obs (if > 0)", &status) ;
        if (slewflg > 0) cf_verbose(1, "SLEW COMMANDED DURING EXPOSURE!!!") ; 

	/* determine the slewing time - note that there is a 16s time
           resolution on the sampling of acaflg. Thus, slews of less than 16s
           are not real */
        slewtime=0 ;
        FITS_update_key(outfits, TINT, "SLEWTIME", &slewtime, 
            "time spent slewing (sec)", &status) ;

        for (i=0; i<npts; i++) 
	  if (acaflg[i] < 1 && time[i] > 0 && time[i] <= exptime ) slewtime++ ;
        if (slewtime > 16) {
             FITS_update_key(outfits, TINT, "SLEWTIME", &slewtime, 
                "time spent slewing (sec)", &status) ;
             slewflg=1 ;
             FITS_update_key(outfits, TINT, "SLEWFLG", &slewflg, 
                "slew commanded during obs (if > 0)", &status) ;
             cf_verbose(1, "SLEW OCCURRED DURING EXPOSURE!!!") ;
             cf_verbose(2, "Time (sec) spent slewing = %d",slewtime) ;
	} 

	/* determine basic properties of the jitter */
	dxrms = 0. ;
        dxrms5 = 0. ;
        dyrms = 0. ;
        dyrms5 = 0. ;
        dxave = 0. ;
        dyave = 0. ;
        nvalid = 0 ;
        nvalid5 = 0 ;
	for (i=0 ; i< npts ; i++) 
	  if(sfknown[i] == 1 && acaflg[i] == 1 && time[i] > qtime && time[i] <= exptime) {
	    cf_verbose(4, "i = %ld, dx = %f, dy = %f", i, dx[i], dy[i]);
            nvalid += 1 ;
            dxrms += dx[i]*dx[i]  ;
            dxave += dx[i] ;
            dyrms += dy[i] * dy[i] ;
            dyave += dy[i] ;
            if (fabs(dx[i]) < 15 && fabs(dy[i]) < 15) {
              dxrms5 += dx[i] * dx[i] ;
              dyrms5 += dy[i] * dy[i] ;
              nvalid5 += 1 ; }
	  }
        cf_verbose(2, "number of points used to determine rms jitter = %d", nvalid) ;
  
	if (nvalid > 0) dxrms = sqrt( dxrms / (float) nvalid) ;
	if (nvalid > 0) dyrms = sqrt( dyrms / (float) nvalid) ;
	if (nvalid5 > 0) dxrms5 = sqrt( dxrms5 / (float) nvalid5) ;
	if (nvalid5 > 0) dyrms5 = sqrt( dyrms5 / (float) nvalid5) ;
        if (nvalid > 0) dxave = dxave / (float) nvalid ;
        if (nvalid > 0) dyave = dyave / (float) nvalid ;
        cf_verbose(2, "jitter properties:") ;
	cf_verbose(2, "dxave = %f, dxrms=%f, dxrms5 = %f", dxave, dxrms, dxrms5) ;
	cf_verbose(2, "dyave = %f, dyrms=%f, dyrms5 = %f", dyave, dyrms, dyrms5) ;
          FITS_update_key(outfits, TFLOAT, "POSAVG_X", &dxave, 
            "[arcsec] mean DX during known track", &status) ;
          FITS_update_key(outfits, TFLOAT, "POSAVG_Y", &dyave, 
            "[arcsec] mean DY during known track", &status) ;

          FITS_update_key(outfits, TFLOAT, "X_JITTER", &dxrms, 
            "[arcsec] sigma of DX during known track", &status) ;
          FITS_update_key(outfits, TFLOAT, "Y_JITTER", &dyrms, 
            "[arcsec] sigma of DY during known track", &status) ;

          FITS_update_key(outfits, TFLOAT, "XJIT_15", &dxrms5, 
            "[arcsec] sigma of DX for err < 15 arcsec", &status) ;
          FITS_update_key(outfits, TFLOAT, "YJIT_15", &dyrms5, 
            "[arcsec] sigma of DY for err < 15 arcsec", &status) ;

  /* determine fraction of exposure where the dy and dx values exceed 2 
     sigma from the mean */
          jitlgx=0 ;
          jitlgy=0 ;
	  for (i=0 ; i<npts ; i++) {
            if ( fabs(dx[i]) > 2.*dxrms) jitlgx++ ;
            if ( fabs(dy[i]) > 2.*dyrms) jitlgy++ ;
		}

	  fjitlgx= (float) jitlgx / (float) npts ;
	  fjitlgy = (float) jitlgy / (float) npts ;

          FITS_update_key(outfits, TFLOAT, "X_JITLRG", &fjitlgx, 
            "frac of DX more than 2-sigma from POSAVE_X", &status) ;
          FITS_update_key(outfits, TFLOAT, "Y_JITLRG", &fjitlgy, 
            "frac of DY more than 2-sigma from POSAVE_Y", &status) ;



  /******* Generate the table and put in guiding info *****/
  
      tfields = 4 ; 
  	 FITS_create_tbl(outfits, BINARY_TBL, 0, tfields, ttype,
		  tform, tunit, textname, &status) ;

      fits_get_colnum(outfits, TRUE, "TIME", &tcol, &status) ;
      FITS_write_col(outfits, TLONG, tcol, 1, 1, npts, time, &status) ;
      fits_get_colnum(outfits, TRUE, "DX", &dxcol, &status) ;
      FITS_write_col(outfits, TFLOAT, dxcol, 1, 1, npts, dx, &status) ;
      fits_get_colnum(outfits, TRUE, "DY", &dycol, &status) ;
      FITS_write_col(outfits, TFLOAT, dycol, 1, 1, npts, dy, &status) ;
      fits_get_colnum(outfits, TRUE, "TRKFLG", &trkcol, &status) ;
      FITS_write_col(outfits, TSHORT, trkcol, 1, 1, npts, trakflg, &status) ;

 /**********************************************************************/

    /* Set the status flag to 0 -- unless it has already been changed. */
    if (dstatus == 1) {
       FITS_movabs_hdu(outfits, 1, &hdutype, &status) ;
       dstatus = 0 ;
       FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
            "Status of jitter data is good", &status) ;
    }

    /* close the output file */
    FITS_close_file(outfits, &status) ; 

    /* free up storage space */
    free(mjd) ;
    free(qvalid) ;
    free(fq1) ;
    free(fq2) ;
    free(fq3) ;
    free(aq1) ;
    free(aq2) ;
    free(aq3) ;
    free(aqsok) ;
    free(acscmd1) ;
    free(acscmd2) ;
    free(acscmd3) ;
    free(time) ;
    free(trakflg) ;
    free(dx) ;
    free(dy) ;

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");

    return 0 ;
} 
