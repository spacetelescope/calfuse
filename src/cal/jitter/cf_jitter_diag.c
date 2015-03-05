/**************************************************************************
 *           Johns Hopkins University  
 *          Center for Astrophysical Sciences    
 *           FUSE
 *************************************************************************
 *
 * synopsis:    cf_jitter_diag input_file raw_file output_file      
 *
 * Description: Reads a FITS file containing housekeeping data and
 *              generates a FITS file containing the jitter information 
 *              for the observation. This program is used to generate 
 *              diagnostics of the jitter performance during an observation
 *
 * Arguments: input_file    Input raw FITS housekeeping file for an observation
 *                            The header of this is used as a template
 *                            for the output file
 *            raw_file      Raw data file - used to get the start time and the 
 *                            exposure time
 *            output_file   FITS table containing the jitter information 
 *     
 * 
 *Contents of output file: 
 *   1. Time of sample (in seconds from start)
 *   2. shifts dx and dy (in arcsec) of the pointing from a reference position
 *   3. tracking quality flag (best=5, worst=0)
 *
 *                                 
 * HISTORY:       01/22/02   v1.0     Robinson  Begin work
 *                08/14/02   v1.1     Robinson  populate JIT_STATUS keyword
 *                09/20/02   v1.11    Robinson  - update FILENAME, FILETYPE keywords
 *                                              - change EXP_DUR keyword to reflect
 *                                                total time duration of jitter data
 *                11/04/02   v1.2     Robinson  change method in which gaps are filled
 *                                              we no longer assume that the data is 
 *                                              tabulated at regular intervals. Also 
 *                                              checks that tabulated times start BEFORE
 *                                              the observation start time
 *                12/17/02   v1.3     Robinson  Slight change in definitions of TRAKFLG
 *                                              values.
 *                01/17/03   v1.4     Robinson  slight change in calculating summary of
 *                                              jitter properties
 *                05/22/03   v1.5     rdr       changed initialization of arrays
 *                06/10/03   v1.6     rdr       Correct problem in determining ref 
 *                                              quaternian
 *                06/23/03   v1.7     rdr       Adjust calculation on some of the 
 *                                               keywords
 *                06/25/03   v1.8     rdr       Added call to raw data file to determine
 *                                               actual start and exposure times
 *                06/27/03   v1.9     rdr       Corrected procedure for determining 
 *                                               reference quaternian
 *                03/03/04   v1.10    bjg       include unistd.h
 *		  06/24/05   v1.11    wvd       Delete unused variables.
 *		  08/24/07   v1.12    bot       changed TINT to TLONG for time l.826.
 *
 **************************************************************************/
#include <unistd.h>
#include <string.h>   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "cf_jitter_diag";   
static char CF_VER_NUM[] = "1.12";  

int quatprod(double *, double *, double *) ;    

int main(int argc, char *argv[])  
{
  FILE *output, *ftest ;
  fitsfile *infits, *outfits ;
  int hdutype, tfields, dstatus ;
  int nvalid, nvalid5, optc, rtime ;
  int aattmode_colnum, pflg=0 ;
  int *aattflg, tflg, qvalidt ;
  int tcolnum, fpdquat1_colnum,  fpdquat2_colnum,  fpdquat3_colnum ;
  int aquat1_colnum, aquat2_colnum, aquat3_colnum, slew_colnum ;
  int status, tcol;
  int dxcol, dycol, trkcol, anynull, intnull=0 ;
  int fpdexp_colnum, qvalid_colnum, cents1_colnum, cents2_colnum, cents3_colnum ;
  int cents4_colnum, cents5_colnum, cents6_colnum, acsflg_colnum ;
  int acscmd1_colnum, acscmd2_colnum, acscmd3_colnum, sfknwn_colnum ;
  int frow, felem, slewflg, ngs_used, jitlgx, jitlgy ;
  int *cents1, *cents2, *cents3, *cents4, *cents5, *cents6 , *qvalid ;
  int *sfknwn, *slew ;
  int cents1t=0, cents2t=0, cents3t=0, cents4t=0, cents5t=0, cents6t=0 ;
  int sfknwnt=0 , slewt=0, slewtime=0 ;
  long i, j, numq4, numq5, numq3, numq2, numqa, numq0 ; 
  long  npts, tval ;
  long *time, qtime=0;
  float tfine, tcoarse, tnull, knowntrk, unknowntrk, fjitlgx, fjitlgy ;
  float targtrk ;
  float ngs1, ngs2, ngs3, ngs4, ngs5, ngs6, gs_used, expdur, exptime ;
  float dxrms, dyrms, dxave, dyave, dxrms5, dyrms5 ;
  float *fpdexp, *acsflg ;
  float *dx , *dy ;
  double *fq1, *fq2,*fq3, *dq1, *dq2, *dq3, *mjd ;
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

char opts[] = "hv:";
    char usage[] =
	"Usage:\n"
	"  cf_jitr_diag [-h] [-v level] hskpfile rawfile outfile\n";
    char option[] = 
	"Options:\n"
	"  -h:  this help message\n"
	"  -v:  verbosity level (=1; 0 is silent)\n";

    verbose_level = 1;

    while ((optc = getopt(argc, argv, opts)) != -1) {
	switch(optc) {

	case 'h':
	    printf("%s\n%s", usage, option);
	    return 0;
	case 'v':
	    verbose_level = atoi(optarg);
	    break;
	}
    }

  cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

    if (argc != optind+3) {
        printf("%s", usage);
        cf_if_error("Incorrect number of arguments");
    }

  cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin execution");

  /* initialize necessary variables */
  hdutype=0 ;
  status=0 ;
  frow=1 ;
  felem=1 ;
  anynull=0 ;


  /* check that the housekeeping and raw data files are available */
  cf_verbose(3,"checking for housekeeping file : %s ",argv[optind]) ;
  ftest=fopen(argv[optind],"r") ;
   if (ftest == NULL) {
      cf_verbose(1,"No input housekeeping file available - returning ") ;
      return 0; }
   fclose( ftest ) ;
  ftest=fopen(argv[optind+1],"r") ;
   if (ftest == NULL) {
      cf_verbose(1,"No input raw data file available - returning ") ;
      return 0; }
   fclose( ftest ) ;

 /* Open the input raw data FITS file and read the start and exposure times */
  cf_verbose(3,"opening the file %s for input ",argv[optind+1]);
  FITS_open_file(&infits, argv[optind+1], READONLY, &status); 
 
  /* read the starting time and duration  from the FITS header */
  FITS_read_key(infits, TDOUBLE, "EXPSTART", &tstrt, buffer, &status) ;
  cf_verbose(3,"start time (MJD) = %15.9f ",tstrt) ;
  FITS_read_key(infits, TFLOAT, "EXPTIME", &exptime, buffer, &status) ;
  exptime=(int) exptime ;
   cf_verbose(3,"exposure duration (sec) = %f ",exptime) ;

   /* close the raw data file and open the housekeeping file */
      FITS_close_file(infits, &status) ;
      FITS_open_file(&infits, argv[optind], READONLY, &status); 
 
      /* check for the existance of the output file */      
      output = fopen(argv[optind+2], "r") ; 
      if (output != NULL) {
        cf_verbose(1,"output file already exists - returning ") ;
        return 1 ; } 
      else fclose(output) ;  

  /* Create the output file */
     FITS_create_file(&outfits, argv[optind+2], &status); 

  /* Copy the header and empty primary image of the input file to 
     the output */
     FITS_copy_hdu(infits, outfits, 0, &status) ;

     /* populate the STATUS keyword with 1 - 
        will reset keyword if everyhing goes well */
     dstatus = 1 ;
       FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
            "problems with jitter data", &status) ;

   /* update FILENAME, FILETYPE, EXPTIME and EXPSTART keywords */
     FITS_update_key(outfits, TSTRING, "FILENAME",argv[optind+2], NULL, &status) ;
     FITS_update_key(outfits, TSTRING, "FILETYPE","JITTER", NULL, &status) ;
     FITS_update_key(outfits, TFLOAT, "EXPTIME",&exptime, NULL, &status) ;
     FITS_update_key(outfits, TDOUBLE, "EXPSTART",&tstrt, NULL, &status) ;

  /* move to the first extension, which contains the table */
     FITS_movabs_hdu(infits, 2, &hdutype, &status) ;

  /* determine the number of times which have been tabulated */
     FITS_read_key(infits, TLONG, "NAXIS2", &npts, NULL, &status) ;
     cf_verbose(3,"number of samples in the table = %d ",npts) ;

  /*   set up arrays to contain the data to be read from the table  */
           mjd = (double *) cf_calloc(npts, sizeof(double)) ;
           fpdexp = (float *) cf_calloc(npts, sizeof(float) ) ;
           qvalid = (int *) cf_calloc(npts, sizeof(int)) ;
           fq1 = (double *) cf_calloc(npts, sizeof(double) ) ;
           fq2 = (double *) cf_calloc(npts, sizeof(double) ) ;
           fq3 = (double *) cf_calloc(npts, sizeof(double) ) ;
           aq1 = (double *) cf_calloc(npts, sizeof(double) ) ;
           aq2 = (double *) cf_calloc(npts, sizeof(double) ) ;
           aq3 = (double *) cf_calloc(npts, sizeof(double)) ; 
	   aattflg = (int *) cf_calloc(npts, sizeof(int)) ;
           acsflg = (float *) cf_calloc(npts, sizeof(float)) ;
           acscmd1 = (double *) cf_calloc(npts, sizeof(double)) ; 
           acscmd2 = (double *) cf_calloc(npts, sizeof(double) ) ; 
           acscmd3 = (double *) cf_calloc(npts, sizeof(double) ) ; 
           cents1 = (int *) cf_calloc(npts, sizeof(int)) ;
           cents2 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents3 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents4 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents5 = (int *) cf_calloc(npts, sizeof(int) ) ;
           cents6 = (int *) cf_calloc(npts, sizeof(int) ) ;
           sfknwn = (int *) cf_calloc(npts, sizeof(int) ) ;
           slew = (int *) cf_calloc(npts, sizeof(int) ) ;

  /* read the data */
	   FITS_get_colnum(infits, CASEINSEN, "MJD", &tcolnum, &status) ;
	   cf_verbose(5,"column number for time = %d ",tcolnum) ; 
	     FITS_read_col(infits, TDOUBLE, tcolnum, frow, felem, npts, &intnull,
			   mjd, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_EXP_DURATION", &fpdexp_colnum, 
                   &status) ;
	   /*printf("column number for fdpexp = %d \n",fpdexp_colnum) ;*/
	     FITS_read_col(infits, TFLOAT, fpdexp_colnum, frow, felem, npts, &intnull,
		     fpdexp, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_MEAS_Q_VALID", &qvalid_colnum, 
                  &status) ;
	   /*printf("column number for qvalid = %d \n",qvalid_colnum) ;*/
	     FITS_read_col(infits, TINT, qvalid_colnum, frow, felem, npts, &intnull,
			   qvalid, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_Q_ECI2BDY_1", &fpdquat1_colnum, 
                  &status) ;
	   /*printf("column number for fq1 = %d \n",fpdquat1_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, fpdquat1_colnum , frow, felem, npts, &intnull,
			  fq1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_Q_ECI2BDY_2", &fpdquat2_colnum, 
                  &status) ;
	   /*printf("column number for fq2 = %d \n",fpdquat2_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, fpdquat2_colnum , frow, felem, npts, &intnull,
			   fq2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_Q_ECI2BDY_3", &fpdquat3_colnum, 
                  &status) ;
	   /*printf("column number for fq3 = %d \n",fpdquat3_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, fpdquat3_colnum , frow, felem, npts, &intnull,
			   fq3, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AATTQECI2BDY_1", &aquat1_colnum, &status) ;
	   /*printf("column number for aq1 = %d \n",aquat1_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, aquat1_colnum , frow, felem, npts, &intnull,
			   aq1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AATTQECI2BDY_2", &aquat2_colnum, &status) ;
	   /*printf("column number for aq2 = %d \n",aquat2_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, aquat2_colnum, frow, felem, npts, &intnull,
			   aq2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AATTQECI2BDY_3", &aquat3_colnum, &status) ;
	   /*printf("column number for aq3 = %d \n",aquat3_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, aquat3_colnum, frow, felem, npts, &intnull,
			   aq3, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AATTMODE", &aattmode_colnum, &status) ;
	   /*printf("column number for aattmode = %d \n",aattmode_colnum) ;*/
	     FITS_read_col(infits, TINT, aattmode_colnum, frow, felem, npts, &intnull,
			   aattflg, &anynull, &status) ;

	   FITS_get_colnum(infits, CASEINSEN, "I_AT_CMD_ATT", &acsflg_colnum, &status) ;
	   /*printf("column number for acsflg = %d \n",acsflg_colnum) ;*/
	     FITS_read_col(infits, TFLOAT, acsflg_colnum, frow, felem, npts, &intnull,
			   acsflg, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AQECI2BDYCMD_1", &acscmd1_colnum, &status) ;
	   /*printf("column number for acscmd1 = %d \n",acscmd1_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, acscmd1_colnum, frow, felem, npts, &intnull,
			   acscmd1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AQECI2BDYCMD_2", &acscmd2_colnum, &status) ;
	   /*printf("column number for acscmd2 = %d \n",acscmd2_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, acscmd2_colnum, frow, felem, npts, &intnull,
			   acscmd2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "AQECI2BDYCMD_3", &acscmd3_colnum, &status) ;
	   /* printf("column number for acscmd3 = %d \n",acscmd3_colnum) ;*/
	     FITS_read_col(infits, TDOUBLE, acscmd3_colnum, frow, felem, npts, &intnull,
			   acscmd3, &anynull, &status) ;

	   FITS_get_colnum(infits, CASEINSEN, "CENTROID1_STATUS", &cents1_colnum,
                &status) ;
	   /* printf("column number for guide star 1 status = %d \n",cents1_colnum) ; */
	     FITS_read_col(infits, TINT, cents1_colnum, frow, felem, npts, &intnull,
			   cents1, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID2_STATUS", &cents2_colnum,
                &status) ;
	   /* printf("column number for guide star 2 status = %d \n",cents2_colnum) ; */
	     FITS_read_col(infits, TINT, cents2_colnum, frow, felem, npts, &intnull,
			   cents2, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID3_STATUS", &cents3_colnum,
                &status) ;
	   /* printf("column number for guide star 3 status = %d \n",cents3_colnum) ; */
	     FITS_read_col(infits, TINT, cents3_colnum, frow, felem, npts, &intnull,
			   cents3, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID4_STATUS", &cents4_colnum,
                &status) ;
	   /* printf("column number for guide star 4 status = %d \n",cents4_colnum) ; */
	     FITS_read_col(infits, TINT, cents4_colnum, frow, felem, npts, &intnull,
			   cents4, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID5_STATUS", &cents5_colnum,
                &status) ;
	   /* printf("column number for guide star 5 status = %d \n",cents5_colnum) ; */
	     FITS_read_col(infits, TINT, cents5_colnum, frow, felem, npts, &intnull,
			   cents5, &anynull, &status) ;
	   FITS_get_colnum(infits, CASEINSEN, "CENTROID6_STATUS", &cents6_colnum,
                &status) ;
	   /* printf("column number for guide star 6 status = %d \n",cents6_colnum) ; */
	     FITS_read_col(infits, TINT, cents6_colnum, frow, felem, npts, &intnull,
			   cents6, &anynull, &status) ;

	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_STARFLD_KNWN", &sfknwn_colnum,
                &status) ;
	   /* printf("column number for starfield known status = %d \n",sfknwn_colnum) ; */
	     FITS_read_col(infits, TINT, sfknwn_colnum, frow, felem, npts, &intnull,
			   sfknwn, &anynull, &status) ;

	   FITS_get_colnum(infits, CASEINSEN, "I_FPD_NEW_CMD_Q", &slew_colnum,
                &status) ;
	   /* printf("column number for slew status = %d \n",slew_colnum) ; */
	     FITS_read_col(infits, TINT, slew_colnum, frow, felem, npts, &intnull,
			   slew, &anynull, &status) ;


  /* close the input file */
      FITS_close_file(infits, &status) ;

  /* check for problems with the tabulated start time */
      if (fabs(mjd[0]-tstrt) > 1) {
	cf_verbose(1,"Problems with tabulated start time ") ;
        cf_verbose(1,"  Value of tabulated start time = %15.7f ",tstrt) ;
        cf_verbose(1,"  First tabulated MJD = %15.7f ",mjd[0]) ;
	/* update JIT_STAT keyword */
          dstatus = 1 ;
          FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
            "problems with tab start time", &status) ;
        FITS_close_file(outfits, &status) ;
	return 1; }

  /* make sure that the tabulates times start BEFORE the start of the observation,
     Otherwise there may be missing data or an error in the tabulated start time */
       tval=(long) (0.5 + (mjd[0]-tstrt) * 86400.) ;
      if (tval > 0) {
       cf_verbose(1,"WARNING: There may be missing data or an error in the start time \n") ;
       cf_verbose(1,"         The first tabulated time is %d seconds from obs start ", tval) ;
	/* update JIT_STAT keyword */
          dstatus = 1 ;
          FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
            "problems with tab start time", &status) ;
        FITS_close_file(outfits, &status) ;
        return 1 ;
 }      


  /* fill in the gaps in the data set 

     - ACS quaternians are generally tabulated every other second
     - FPD quaternians can be every second or every other second
     - AATTMODE flags tabulates every 5 seconds
     - ACS flags (acsflg) - tabulated every 16 seconds

NOTE:  there can be irregularities in the spacing. Thus, save the last
tabulated value and use this in filling in missing data for a given time */

      /* ACS quaternians */
      aq1t = 0. ;
      aq2t = 0. ;
      aq3t = 0. ;
      for (i=0 ; i<npts ; i++) {
        if ( aq1[i] > -1 && aq2[i] > -1){
	  aq1t=aq1[i] ;
	  aq2t=aq2[i] ;
	  aq3t=aq3[i] ; }
        else {
          aq1[i] = aq1t ;
          aq2[i] = aq2t ;
          aq3[i] = aq3t ; }
	}

      /* FPD quaternians and qvalid flags 
         only update the one following a good quaternian. if there is another one
         with no quaternian after that, do not update it */

     fq1t = 0. ;
     fq2t = 0. ;
     fq3t = 0. ;
     qvalidt = 0 ;
     tflg = -1 ;
      for (i=0 ; i<npts ; i++) {
        if ( (fabs(fq1[i]) > 0.) && (fabs(fq2[i]) > 0.) ){
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
          sfknwnt = sfknwn[i] ;
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
          sfknwn[i] = sfknwnt ;
          slew[i] = slewt ;
          tflg = -1 ; }
	}


       /* AATTMODE flag - update as long as quaternians are
          present */
 
       tflg = -1 ;          
       for (i=0 ; i<npts-1 ; i++) {
         if (aattflg[i] > -1) tflg = aattflg[i] ;
         else {
	     aattflg[i] = -10 ;
             if (aq1[i] > -1) aattflg[i]=3 ;
      	     if ((fabs(fq1[i]) > 0.) && ( fabs(fq2[i]) > 0. ) ) aattflg[i]=tflg ; }
        } 


       /* ACSFLG flag */

       tflg = -10 ;          
       for (i=0 ; i<npts-1 ; i++) {
         if (acsflg[i] > -1) tflg = acsflg[i] ;
         else acsflg[i] = tflg ;
        } 
		       

 /********** derive the final quats arrays and the tracking flags *******/

  /*    set up arrays to contain the data */
         time = (long *) cf_calloc(npts, sizeof(long) ) ;
	 trakflg = (short *) cf_calloc(npts, sizeof(short) ) ; 
         dq1 = (double *) cf_calloc(npts, sizeof(double) ) ;
         dq2 = (double *) cf_calloc(npts, sizeof(double)) ;
         dq3 = (double *) cf_calloc(npts, sizeof(double) ) ;
         dx = (float *) cf_calloc(npts, sizeof(float) ) ;
         dy = (float *) cf_calloc(npts, sizeof(float) ) ;

 /* fill in the time array with times (in seconds) from the start of the observation */
	 for (i=0 ; i<npts ; i++) {
           tval=(long) (0.5 + (mjd[i]-tstrt) * 86400.) ;
           time[i] = tval ; }


 /* determine the reference quat by looking at the commanded ACS quaternian (in arrays
     acscmd1, acscmd2 and acscmd3) when they have converged to the fpd commanded
     quaternian - i.e. when the AT_CMD_ATT flag is 1 (held in the array acsflg) -
     Make sure that the time for convergence is AFTER the start of the observation and that
     the starfield is known (i.e. sfknwn = 1). Also, because of the 16s sample time of 
     the ascflg array, require that the condition exists for at least 30 consequtive seconds 
     before accepting the quaternians.
       qtime = time at which the quaternians were found. This is the earliest time in the exposure
               for which reliable positions are available*/

	 rtime=0 ;
	 for (i=0 ; i<npts; i++) {
           if (acsflg[i] == 1 &&  time[i] > 0 && sfknwn[i] == 1) {
             rtime += 1 ;
             cf_verbose(4,"quat ref cond exists at time=%d, rtime=%d ",time[i],rtime) ;
	     if(rtime == 1) qtime=time[i] ;
             if (rtime > 30 && acscmd1[i] > -1 && acscmd2[i] > -1 && acscmd3[i] > -1 ) {
  	       quat_ref[0]=acscmd1[i] ;
               quat_ref[1]=acscmd2[i] ;
               quat_ref[2]=acscmd3[i] ;
               quat_ref[3]=sqrt( 1. - quat_ref[0]*quat_ref[0] - quat_ref[1]*quat_ref[1] - quat_ref[2]*quat_ref[2] ) ;
	       goto skip ; }
	   }
           else rtime=0 ; 
	 }

     /* if we are here, then the flag was never set and the reference quat not determined */
	 cf_verbose(1," Problem determining the reference quaternians ") ;
	 /* determine cause of failure */
         pflg = 0 ;
	 for (i=0 ; i<npts; i++) if (acsflg[i] == 1) pflg=1 ;
	   if (pflg == 0) {
               cf_verbose(1," AT_CMD_ATT flag was never set  ") ;
         	/* update JIT_STAT keyword */
                   dstatus = 2 ;
                   FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
                   "AT_CMD_ATT flg never set", &status) ;
	   }
         pflg = 0 ;
	 for (i=0 ; i<npts; i++) if (sfknwn[i] == 1) pflg=1 ;
	   if (pflg == 0)  {
               printf(" UNKNOWN STAR FIELD \n") ;
   	         /* update JIT_STAT keyword */
                    dstatus = 3 ;
                    FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
                      "unknown star field", &status) ;
	   }
         FITS_close_file(outfits, &status) ;
         return 1 ;
 skip: ;
	cf_verbose(3," ") ;
        cf_verbose(3,"reference quat values - from ACS commanded pos ") ;
        cf_verbose(3,"found at time %d \n",qtime) ;
      for (i=0; i<4 ; i++) cf_verbose(3,"for i=%d, quat=%18.15f ",i,quat_ref[i]) ;              

      /* determine the value of the tracking flags and identify which quaternians to use */

	 for (i=0; i<npts ; i++) {
           if (fq1[i] > -1 && fq2[i] > -1 && qvalid[i] > 0) {
             dq1[i] = fq1[i] ; 
             dq2[i] = fq2[i] ;
             dq3[i] = fq3[i] ;
             if (sfknwn[i] == 1 && aattflg[i] == 5) trakflg[i] = (short) 5 ;
             if (sfknwn[i] == 0 && aattflg[i] == 5) trakflg[i] = (short) 4 ; 
             if (aattflg[i] == 4) trakflg[i] = (short) 3 ; }
	   else if (aq1[i] > -1 && aq2[i] > -1 ) {
             dq1[i] = aq1[i] ;
             dq2[i] = aq2[i] ;
             dq3[i] = aq3[i] ;
             if(aattflg[i] == 4) trakflg[i] = (short) 2 ;
             if(aattflg[i] == 3) trakflg[i] = (short) 1 ;}
           else {
             dq1[i] = -10. ;
             dq2[i] = -10. ;
             dq3[i] = -10. ; 
             trakflg[i]= (short) 0 ; }
       	   }


      /* determine the shift in the quaturnians and derive the changes
	 in x and y which they represent :
           DX (from q2 values : dq2=dx/2) and DY (from q1 values : dq1=-dy/2)
           in arcseconds - the conversion factor used is the number of arcsec
            per radian * 2 */
       for (i=0l ; i<npts ; i++) {
	 if (trakflg[i] > 0) {
            tquat[0] = dq1[i] ;
            tquat[1] = dq2[i] ;
            tquat[2] = dq3[i] ;
            tquat[3] = sqrt(1. - dq1[i]*dq1[i] - dq2[i]*dq2[i] - dq3[i]*dq3[i] ) ;

	  /* initialize dquat */
	    for (j=0; j<4 ; j++) dquat[j]=0. ;

 	  /* now determine the change in quat */
	     quatprod(tquat, quat_ref, dquat) ;

 	  /* convert of arcsec in x and y */
 	    dy[i] = (float) -dquat[0] * 412541. ;
	    dx[i] = (float) dquat[1] * 412541. ; } 
	 else {
	    dx[i] = -999. ;
            dy[i] = -999. ; }
    }

/****** provide keywords in the top level header which describe the observations ****/

    /* determine the time span in the jitter file (in seconds) and set EXP_DUR to
        this value  */
       expdur = (float) (time[npts-1] - time[0]) + 1. ;
       FITS_update_key(outfits, TFLOAT, "EXP_DUR", &expdur, 
         "duration of jitter data (sec)", &status) ;

       /* check for a short exposure and reset exptime if necessary */
       if (time[npts-1] < exptime) {
	 exptime = time[npts-1] ;
         cf_verbose(2,"short exposure- reset exptime to %f ",exptime) ;
       }

	 /* determine the statistics of the guiding */
         numq4=0 ;
         numq5=0 ;
         numq3=0 ;
         numq2=0 ;
         numqa=0 ;
         numq0=0 ;
	 cf_verbose(3,"exposure duration = %f ",expdur) ;

	 for (i=0; i<npts ; i++) {
	   cf_verbose(6,"at i=%d, trakflg=%d, qvalid=%d, time=%f ",
	      i,trakflg[i], qvalid[i], time[i] ) ;
           if(trakflg[i] == 5 && qvalid[i] > 0 && time[i] > 0 && time[i] <= exptime) 
                   numq5+=1 ;
           if(trakflg[i] == 4 && time[i] > 0 && time[i] <= exptime) numq4+=1 ;
           if(trakflg[i] == 3 && time[i] > 0 && time[i] <= exptime) numq3+=1 ;
           if(trakflg[i] == 2 && time[i] > 0 && time[i] <= exptime) numq2+=1 ;
           if(trakflg[i] == 1 && time[i] > 0 && time[i] <= exptime) numqa+=1 ;
           if(trakflg[i] == 0 && time[i] > 0 && time[i] <= exptime) numq0+=1 ; }

         cf_verbose(3," ") ;
         cf_verbose(3,"duration of jitter data =%f, exposure time = %f ",expdur,exptime) ;
	 cf_verbose(3,"number of samples in fine track with known GS = %d ",numq5) ;
	 cf_verbose(3,"number of samples in fine track with unknown GS = %d ",numq4) ;
         cf_verbose(3,"number of samples in coarse track = %d ",numq3) ;
         cf_verbose(3,"number of samples using ACS quaternians = %d ",numq2) ;
         cf_verbose(3,"number of samples with no quaternians = %d ",numq0) ;


       /* tfine, tcoarse and tnull are the fraction of the total exposure
          when the observation is in fine track, coarse track and no track */
       tfine = (float) (numq5+numq4) /  exptime ;
       tcoarse = (float) (numq3+numq2+numqa) /  exptime ;
       tnull = (float) numq0 /  exptime ;
       FITS_update_key(outfits, TFLOAT, "FINE_GDE", &tfine, 
           "fraction of exp in fine guide", &status) ;
       FITS_update_key(outfits, TFLOAT, "COARSGDE", &tcoarse, 
           "fraction of exp in coarse guide", &status) ;
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
         if (cents1[i] == 2 && time[i] > 0 && time[i] <= exptime) ngs1+=1. ;
         if (cents2[i] == 2 && time[i] > 0 && time[i] <= exptime) ngs2+=1. ;
         if (cents3[i] == 2 && time[i] > 0 && time[i] <= exptime) ngs3+=1. ;
         if (cents4[i] == 2 && time[i] > 0 && time[i] <= exptime) ngs4+=1. ;
         if (cents5[i] == 2 && time[i] > 0 && time[i] <= exptime) ngs5+=1. ;
         if (cents6[i] == 2 && time[i] > 0 && time[i] <= exptime) ngs6+=1. ; }
       /* convet to fraction of observing time and put into keywords */
       ngs1 = ngs1 / exptime ;
        FITS_update_key(outfits, TFLOAT, "GS1_USED", &ngs1, 
           "fraction of exp using guide star 1", &status) ;
       ngs2 = ngs2 / exptime ;
        FITS_update_key(outfits, TFLOAT, "GS2_USED", &ngs2, 
           "fraction of exp using guide star 2", &status) ;
       ngs3 = ngs3 / exptime ;
        FITS_update_key(outfits, TFLOAT, "GS3_USED", &ngs3, 
           "fraction of exp using guide star 3", &status) ;
       ngs4 = ngs4 / exptime ;
        FITS_update_key(outfits, TFLOAT, "GS4_USED", &ngs4, 
           "fraction of exp using guide star 4", &status) ;
       ngs5 = ngs5 / exptime ;
        FITS_update_key(outfits, TFLOAT, "GS5_USED", &ngs5, 
           "fraction of exp using guide star 5", &status) ;
       ngs6 = ngs6 / exptime ;
        FITS_update_key(outfits, TFLOAT, "GS6_USED", &ngs6, 
           "fraction of exp using guide star 6", &status) ;

	/* determine the number of guide stars used during the observation */
	ngs_used = 0 ;
	if (ngs1 > 0) ngs_used += 1 ; 
	if (ngs2 > 0) ngs_used += 1 ; 
	if (ngs3 > 0) ngs_used += 1 ; 
	if (ngs4 > 0) ngs_used += 1 ; 
	if (ngs5 > 0) ngs_used += 1 ; 
	if (ngs6 > 0) ngs_used += 1 ; 
        FITS_update_key(outfits, TINT, "NGS_USED", &ngs_used, 
           "number of guide stars used", &status) ;
        cf_verbose(3," ") ;
	cf_verbose(3,"number of guide stars used = %d \n",ngs_used) ;


	/* determine times of known and unknown tracking and tracking on the target*/ 
	knowntrk = 0. ;
        unknowntrk = 0. ;
        targtrk = 0. ;
        gs_used = 0. ;
	for (i=0 ; i< npts ; i++) {
          if (qvalid[i] == 1 && time[i] > 0 && time[i] <= exptime) gs_used+=1. ;
          if (sfknwn[i] == 1 && qvalid[i] ==1 && time[i] > 0 && time[i] <= exptime) 
                knowntrk+=1. ;
          if (sfknwn[i] == 0 && qvalid[i] ==1 && time[i] > 0 && time[i] <= exptime) 
                unknowntrk+=1. ; 
          if (sfknwn[i] == 1 && acsflg[i] == 1 && time[i] > 0 && time[i] <= exptime)
	    targtrk += 1. ;
     }

	/* check for digitization errors in targtrk (caused by the 16s cadence of acsflg) */
	if (targtrk+16 > exptime) targtrk = exptime ;

        cf_verbose(3,"seconds when guide stars are used = %f",gs_used ) ;
        cf_verbose(3,"seconds tracking on known stars = %f ",knowntrk) ;
        cf_verbose(3,"seconds tracking on unknown stars = %f ",unknowntrk) ;
        cf_verbose(3,"seconds tracking on the target = %f ",targtrk ) ;
	knowntrk = knowntrk /  exptime ;
	unknowntrk = unknowntrk / exptime ;
        targtrk = targtrk / exptime ;
        gs_used = gs_used / exptime ;
         FITS_update_key(outfits, TFLOAT, "KNOWNTRK", &knowntrk, 
            "fraction of exp tracking on known stars", &status) ;
         FITS_update_key(outfits, TFLOAT, "UNKWNTRK", &unknowntrk, 
            "fraction of exp tracking on unknown stars", &status) ;
        FITS_update_key(outfits, TFLOAT, "TARGTRK", &targtrk, 
            "fraction of exp tracking on target", &status) ;
        FITS_update_key(outfits, TFLOAT, "GS_INUSE", &gs_used, 
            "fraction of exp tracking on guide stars", &status) ;


	/* determine whether a new attitude is commanded during the exposure */
        slewflg = -1 ;
	for (i=0; i<npts ; i++) 
          if(time[i] > 0 && slew[i] > 0 && time[i] <= exptime) slewflg = 1 ;
        FITS_update_key(outfits, TINT, "SLEWFLG", &slewflg, 
            "slew commanded during obs (if > 0)", &status) ;
        if (slewflg > 0) {
	  cf_verbose(1," ") ;
          cf_verbose(1,"SLEW COMMANDED DURING OBSERVATION!!! ") ; }

	/* determine the slewing time - note that there is a 16s time
           resolution on the sampling of acsflg. Thus, slews of less than 16s
           are not real */
        slewtime=0 ;
        FITS_update_key(outfits, TINT, "SLEWTIME", &slewtime, 
            "time spent slewing (sec)", &status) ;

        for (i=0; i<npts; i++) 
	  if (acsflg[i] < 1 && time[i] > 0 && time[i] <= exptime ) slewtime+=1 ;
        if (slewtime > 16) {
             FITS_update_key(outfits, TINT, "SLEWTIME", &slewtime, 
                "time spent slewing (sec)", &status) ;
             slewflg=1 ;
             FITS_update_key(outfits, TINT, "SLEWFLG", &slewflg, 
                "slew commanded during obs (if > 0)", &status) ;
             cf_verbose(3,"SLEW OCCURRED DURING OBSERVATION!!! ") ;
             cf_verbose(3,"Time (sec) spent slewing = %d \n",slewtime) ;
	} 

	/* determine basic properties of the jitter 
            - these properties are only determined during the time of the exposure where the 
              spacecraft is guiding on the target */
	dxrms = 0. ;
        dxrms5 = 0. ;
        dyrms = 0. ;
        dyrms5 = 0. ;
        dxave = 0. ;
        dyave = 0. ;
        nvalid = 0 ;
        nvalid5 = 0 ;
	for (i=0 ; i< npts ; i++) 
	  if(sfknwn[i] == 1 && acsflg[i] == 1 && time[i] > qtime  && time[i] <= exptime) {
            nvalid += 1 ;
            dxrms = dxrms + dx[i]*dx[i]  ;
            dxave = dxave + dx[i] ;
            dyrms = dyrms + dy[i] * dy[i] ;
            dyave = dyave + dy[i] ;
            if (time[i] > exptime -  300) {
              dxrms5 = dxrms5 + dx[i] * dx[i] ;
              dyrms5 =dyrms5 + dy[i] * dy[i] ;
              nvalid5 += 1 ; }
	  }
        cf_verbose(3,"number of points used to determine rms jitter = %d ",nvalid) ;
  
	if (nvalid > 0) dxrms = sqrt( dxrms / (float) nvalid) ;
	if (nvalid > 0) dyrms = sqrt( dyrms / (float) nvalid) ;
	if (nvalid5 > 0) dxrms5 = sqrt( dxrms5 / (float) nvalid5) ;
	if (nvalid5 > 0) dyrms5 = sqrt( dyrms5 / (float) nvalid5) ;
        if (nvalid > 0) dxave = dxave / (float) nvalid ;
        if (nvalid > 0) dyave = dyave / (float) nvalid ;
        cf_verbose(3," ") ;
        cf_verbose(3,"jitter properties ") ;
	cf_verbose(3,"dxave = %f, dxrms=%f, dxrms5 = %f ", dxave, dxrms, dxrms5) ;
	cf_verbose(3,"dyave = %f, dyrms=%f, dyrms5 = %f ", dyave, dyrms, dyrms5) ;
          FITS_update_key(outfits, TFLOAT, "POSAVG_X", &dxave, 
            "[arcsec] mean DX during exposure", &status) ;
          FITS_update_key(outfits, TFLOAT, "POSAVG_Y", &dyave, 
            "[arcsec] mean DY during exposure", &status) ;

          FITS_update_key(outfits, TFLOAT, "X_JITTER", &dxrms, 
            "[arcsec] sigma of DX during exposure", &status) ;
          FITS_update_key(outfits, TFLOAT, "Y_JITTER", &dyrms, 
            "[arcsec] sigma of DY during exposure", &status) ;

          FITS_update_key(outfits, TFLOAT, "X_JIT_5M", &dxrms5, 
            "[arcsec] sigma of DX during last 5 min of exp", &status) ;
          FITS_update_key(outfits, TFLOAT, "Y_JIT_5M", &dyrms5, 
            "[arcsec] sigma of DY during last 5 min of exp", &status) ;

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

      /* set the status flag as good */
       FITS_movabs_hdu(outfits, 1, &hdutype, &status) ;
       dstatus = 0 ;
       FITS_update_key(outfits, TINT, "JIT_STAT", &dstatus, 
            "status of jitter data is good", &status) ;

      /* close the output file */
    FITS_close_file(outfits, &status) ; 

    /* free up storage space */
    free(mjd) ;
    free(fpdexp) ;
    free(qvalid) ;
    free(fq1) ;
    free(fq2) ;
    free(fq3) ;
    free(aq1) ;
    free(aq2) ;
    free(aq3) ;
    free(aattflg) ;
    free(acscmd1) ;
    free(acscmd2) ;
    free(acscmd3) ;
    free(time) ;
    free(trakflg) ;
    free(dq1) ;
    free(dq2) ;
    free(dq3) ;
    free(dx) ;
    free(dy) ;

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");

    return 0 ;

} 


/* procedure to do a matrix multiplication of two quaturnians to determine
   the displacement of a measured value from a reference */
 
int quatprod(double *tquat, double *quat_ref, double *dquat) {

  int i, j ;
  double multmat[4][4], qinv[4] ; 

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
  

  /* specify the inverse of quat_ref */
  qinv[0] = -quat_ref[0] ;
  qinv[1] = -quat_ref[1] ;
  qinv[2] = -quat_ref[2] ;
  qinv[3] =  quat_ref[3] ;
 
  /* do the matrix multiplication */
  for (i=0 ; i<4 ; i++) {
    dquat[i] = 0. ;
     for (j=0; j<4 ; j++) dquat[i] = dquat[i] + multmat[j][i]*qinv[j] ; 
}

  return 0 ;

}
