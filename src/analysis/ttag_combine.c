/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	ttag_combine outfile file1 file2 file3 file4 file5 ...
 *
 * Description: Reads multiple contiguous ttag files and concatenates them 
 *              into a single file.  For now, the input files should be
 *              listed in time order.
 *
 * History:	08/16/99	emurphy	 Begin work
 *              08/27/99        emurphy  Added comments.  Good time intervals
 *                                       now match observing times.
 *              10/09/99        emurphy  exptime now sum of individual exptimes
 *              08/20/01 v1.2	wvd	 Complain if FPA moves between exposures.
 *					 Set MAXGTI = 2000.
 *              10/01/01 v1.3   RDR      Update all segment and active image counters
 *                                       DATE keyword now contains the date 
 *                                       NEVENTS reflects total number of events in
 *                                         combined file
 *                                       EXP_ID defaults to 999 to signify combined
 *                                         data
 *                                       Print warning if FPA x position changes
 *                                         by more than 30 microns
 *                                       PLANTIME now reflects total planned exposure
 *                                       Will transfer DX, DY and DNFLG columns
 *                                         if available 
 *              10/09/01 v1.4	wvd	 Delete redundant check for FPA splits.
 *					 Program complains if input files not in
 *					 time order.
 *              01/15/02 v1.5	wvd	 Update keyword EXPNIGHT.
 *              02/13/02 v1.6	wvd	 Respond gracefully if EXPNIGHT does not exist.
 *              02/20/02     	wvd	 Update keyword RAWTIME, if it exists.
 *              04/26/02 v1.7	wvd	 Change vtime to type time_t
 *              12/11/02 v1.8	wvd	 Change EXPNIGHT to type long
 *              02/20/03 v1.9   wvd      Correct errors in pointer use.
 *					 Change abs() to fabs().
 *              03/22/04 v1.10  bjg      Change cf_fpa_pos() to cf_read_fpa_pos()
 *              04/05/04 v1.11  bjg      Include math.h
 *                                       Remove unused variables
 *                                       Change formats to match arg type
 *                                       in printf
 *
 *                                       
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "calfuse.h"

static char CF_PRGM_ID[] = "ttag_combine";
static char CF_VER_NUM[] = "1.11";

#define MAXGTI          2000      /*  Maximum number of good time intervals */

static int 
cf_dtcor_read_counters(fitsfile *infits, long *de, long *fe, long *ai1a,
   long *ai1b, long *ai2a, long *ai2b, long *sic, long *lif, long *as,
   double *ctm )
{
    char cntstr[FLEN_CARD], detstr[FLEN_CARD];
    long de_b, de_e, fe_b, fe_e;
    long ai1a_b, ai1a_e, ai1b_b, ai1b_e, ai2a_b, ai2a_e, ai2b_b, ai2b_e ;
    long sic_b, sic_e, lif_b, lif_e, as_b, as_e ;
    double ctm_b, ctm_e ;
    int status=0;

    FITS_read_key(infits, TSTRING, "DETECTOR", detstr, NULL, &status);

    /* Segment 1A Counters */

    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"DE_B");
    FITS_read_key(infits, TLONG, cntstr, &de_b, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"DE_E");
    FITS_read_key(infits, TLONG, cntstr, &de_e, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"FE_B");
    FITS_read_key(infits, TLONG, cntstr, &fe_b, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"FE_E");
    FITS_read_key(infits, TLONG, cntstr, &fe_e, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"SIC_B");
    FITS_read_key(infits, TLONG, cntstr, &sic_b, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"SIC_E");
    FITS_read_key(infits, TLONG, cntstr, &sic_e, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"LIF_B");
    FITS_read_key(infits, TLONG, cntstr, &lif_b, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"LIF_E");
    FITS_read_key(infits, TLONG, cntstr, &lif_e, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"AS_B");
    FITS_read_key(infits, TLONG, cntstr, &as_b, NULL, &status);

    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"AS_E");
    FITS_read_key(infits, TLONG, cntstr, &as_e, NULL, &status);

    /* Active Image Counters */
    FITS_read_key(infits, TLONG, "C1AAI_B", &ai1a_b, NULL, &status);
    FITS_read_key(infits, TLONG, "C1AAI_E", &ai1a_e, NULL, &status);
    FITS_read_key(infits, TLONG, "C1BAI_B", &ai1b_b, NULL, &status);
    FITS_read_key(infits, TLONG, "C1BAI_E", &ai1b_e, NULL, &status);
    FITS_read_key(infits, TLONG, "C2AAI_B", &ai2a_b, NULL, &status);
    FITS_read_key(infits, TLONG, "C2AAI_E", &ai2a_e, NULL, &status);
    FITS_read_key(infits, TLONG, "C2BAI_B", &ai2b_b, NULL, &status);
    FITS_read_key(infits, TLONG, "C2BAI_E", &ai2b_e, NULL, &status);
    FITS_read_key(infits, TDOUBLE, "CTIME_B", &ctm_b, NULL, &status);
    FITS_read_key(infits, TDOUBLE, "CTIME_E", &ctm_e, NULL, &status);


    *de = de_e - de_b;
    *fe= fe_e - fe_b;
    *ai1a = ai1a_e - ai1a_b ;
    *ai1b = ai1b_e - ai1b_b ;
    *ai2a = ai2a_e - ai2a_b ;
    *ai2b = ai2b_e - ai2b_b ;
    *sic = sic_e - sic_b ;
    *lif = lif_e - lif_b ;
    *as = as_e - as_b ; 
    *ctm = ctm_e - ctm_b ;

    return 0;
}

int main(int argc, char *argv[])
{
    int      i, j, k;
    char     buffer[FLEN_CARD], comment[FLEN_CARD], stime[FLEN_CARD];
    char     cntstr[FLEN_CARD], detstr[FLEN_CARD];
    char     expid[FLEN_CARD], procstep[FLEN_CARD] ;
    long     de, fe, de_tot, fe_tot, fe_b, de_b, fe_e, de_e;
    long     ai1a, ai1b, ai2a, ai2b, ai1a_tot, ai1b_tot, ai2a_tot;
    long     ai2b_tot,ai1a_b, ai1a_e, ai1b_b, ai1b_e, ai2a_b ;
    long     ai2a_e, ai2b_b, ai2b_e, sic_b, sic_e, sic_tot, sic ;
    long     lif_b, lif_e, lif_tot, lif, as_b, as_e, as_tot, as ;
    long     expnight, indnight;
    int      status=0, nfiles=0, hdutype=0, nkeys, keynum, anynull;
    fitsfile *infits, *outfits;
    long     nevents, jrow, ngti, ngti_in;
    char     *pha_in, *dnflg ;
    short    *x_in, *y_in;
    float    *time_in, *gti_start, *gti_stop;
    float    *gti_start_in, *gti_stop_in, *dx, *dy ;
    float    fpasx0, fpasx, fpalx0, fpalx, dfpasx, dfpalx ;
    double   mjd_init, mjd_start, mjd_end, exptime, mjd_last=0, indtime;
    double   rawtime, trawtime;
    double   ctm, ctm_b, ctm_e, ctm_tot ;
    int      plantime, procflg ;
    long     tplantime ;
    int      tcol_in,  xcol_in,  ycol_in,  phacol_in, startcol ,stopcol;
    int      tcol_out, xcol_out, ycol_out, phacol_out;
    int      dxcol_in, dxcol_out, dycol_in, dycol_out,dnflgcol_in;
    int      dnflgcol_out ;
    int      fpa_split = FALSE;
    time_t   vtime;


    if (argc < 3) {
        printf("Incorrect number of arguments.\n");
	printf("Calling sequence:ttag_combine output_file input_file1 input_file2 input_file3 ...\n");
        exit(1);
    }

    /* Initialize error checking. */
    cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);
    
    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Started execution.");

    /* get and display time */
    vtime = time(NULL) ;
    strcpy(stime,ctime(&vtime)) ;

    jrow=1;
    de_tot=0 ;
    fe_tot=0 ;
    ai1a_tot=0 ;
    ai1b_tot=0 ;
    ai2a_tot=0 ;
    ai2b_tot=0 ;
    sic_tot=0 ;
    lif_tot=0 ;
    as_tot=0 ;
    ctm_tot=0. ;
    procflg=-1 ;
    gti_start  = (float *) cf_malloc(sizeof(float) * 1000);
    gti_stop   = (float *) cf_malloc(sizeof(float) * 1000);

    gti_start[0] = 0.0;
    gti_stop[0] = 0.0;

    nfiles = argc-1;

    exptime=0.0;
    expnight=0;
    trawtime=0.0;
    tplantime=0;

    ngti=0; /* is the number of new good time intervals */

    for (i=2; i<=nfiles; i++) {
             printf("Processing %-60.60s\n",argv[i]);
             FITS_open_file(&infits, argv[i], READONLY, &status);

             if (i==2) {
                 FITS_read_key(infits, TDOUBLE, "EXPSTART", &mjd_init, 
                     NULL, &status);
                 mjd_last=mjd_start=mjd_init;
                 FITS_read_key(infits, TDOUBLE, "EXPEND", &mjd_end, 
                     NULL, &status);
                 FITS_read_key(infits, TSTRING, "INIT_COR", procstep,
                     NULL, &status) ;
		 cf_read_fpa_pos (infits, &fpasx0, &fpalx0);
                 FITS_create_file(&outfits, argv[1], &status);
                 FITS_copy_hdu(infits, outfits, 0, &status);
                 FITS_create_hdu(outfits, &status);
                 FITS_movabs_hdu(infits,2,&hdutype,&status);
                 FITS_get_hdrpos(infits, &nkeys, &keynum, &status);
                 for (j=1; j<=nkeys; j++) {
                     FITS_read_record(infits, j, buffer, &status);
                     FITS_write_record(outfits, buffer, &status);
                 }
	      /* adjust some initial keywords */
                 nevents=0;
                 FITS_update_key(outfits, TLONG, "NAXIS2", &nevents, 
                                 NULL, &status);
	      /* set a flag (procflg) if the file is processed in any way */
		 if ( strncmp(procstep, "COMP", 4) == 0) {
                    printf("file has been processed - setting flag \n") ;
                    procflg = 1 ; }
                  else printf("These are raw data files \n") ;
     	     } else {
                 FITS_read_key(infits, TDOUBLE, "EXPSTART", &mjd_start, 
                     NULL, &status);
                 if (mjd_start < mjd_last) 
                     cf_if_error("Files must be in time order, first to last.");
		 mjd_last=mjd_start;
                 FITS_read_key(infits, TDOUBLE, "EXPEND", &mjd_end, 
                     NULL, &status);
		 cf_read_fpa_pos (infits, &fpasx, &fpalx);
                 dfpasx = fpasx - fpasx0 ;
                 dfpalx = fpalx - fpalx0 ;
		 if (fabs((double)dfpasx) > 30.) {
                    printf("\nLarge change in SIC FPA X position: %f \n",dfpasx) ;
		    fpa_split = TRUE;
		 }
		 if (fabs((double)dfpalx) > 30.) {
                    printf("\nLarge change in LIF FPA X position: %f \n",dfpalx) ;
		    fpa_split = TRUE;
		 }
	     }

             FITS_movabs_hdu(infits,1,&hdutype,&status);

             FITS_read_key(infits, TDOUBLE, "EXPTIME", &indtime, 
                           NULL, &status);
             exptime+=indtime;

	     fits_read_key(infits, TDOUBLE, "RAWTIME", &rawtime,
                           NULL, &status);
	     if (status) 
		status = 0;
	     else
                trawtime+=rawtime;

             FITS_read_key(infits, TINT, "PLANTIME", &plantime, 
                           NULL, &status);
             tplantime+=plantime ;

	     fits_read_key(infits, TLONG, "EXPNIGHT", &indnight,
                           NULL, &status);
	     if (status) 
		status = 0;
	     else
                expnight+=indnight;

       cf_dtcor_read_counters(infits,&de,&fe,&ai1a,&ai1b,&ai2a,&ai2b,
          &sic,&lif,&as,&ctm) ;

             de_tot += de;
             fe_tot += fe;
             ai1a_tot += ai1a ;
             ai1b_tot += ai1b ;
             ai2a_tot += ai2a ;
             ai2b_tot += ai2b ;
             sic_tot += sic ;
             lif_tot += lif ;
             as_tot += as ;
             ctm_tot += ctm ;

             FITS_movabs_hdu(infits,2,&hdutype,&status);
             /* Get the number of events in the input file. */
             FITS_read_key(infits, TLONG, "NAXIS2", &nevents, NULL, &status);

             printf("     EXPTIME=%10.3f  NEVENTS=%10ld\n",indtime,nevents);

             /* Allocate space for the input and output times, coords, 
                pha */
             time_in  = (float *) cf_malloc(sizeof(float) * nevents);
             x_in     = (short *) cf_malloc(sizeof(short) * nevents);
             y_in     = (short *) cf_malloc(sizeof(short) * nevents);
             pha_in   = (char *)  cf_malloc(sizeof(char) * nevents);
             dx       = (float *) cf_malloc(sizeof(float) * nevents);
             dy       = (float *) cf_malloc(sizeof(float) * nevents);
             dnflg    = (char *) cf_malloc(sizeof(char) * nevents);

             FITS_get_colnum(infits, TRUE, "TIME", &tcol_in, &status);
             FITS_get_colnum(infits, TRUE, "X",    &xcol_in, &status);
             FITS_get_colnum(infits, TRUE, "Y",    &ycol_in, &status);
             FITS_get_colnum(infits, TRUE, "PHA",  &phacol_in, &status);
	   if(procflg > 0) {
            FITS_get_colnum(infits, TRUE, "DX",   &dxcol_in, &status);
            FITS_get_colnum(infits, TRUE, "DY",   &dycol_in, &status);
            FITS_get_colnum(infits, TRUE, "DNFLG",&dnflgcol_in, &status); }

             FITS_get_colnum(outfits, TRUE, "TIME", &tcol_out, &status);
             FITS_get_colnum(outfits, TRUE, "X",    &xcol_out, &status);
             FITS_get_colnum(outfits, TRUE, "Y",    &ycol_out, &status);
             FITS_get_colnum(outfits, TRUE, "PHA",  &phacol_out, &status);
	   if (procflg > 0 ) {
            FITS_get_colnum(outfits, TRUE, "DX",   &dxcol_out, &status);
            FITS_get_colnum(outfits, TRUE, "DY",   &dycol_out, &status);
            FITS_get_colnum(outfits, TRUE, "DNFLG",&dnflgcol_out, &status);}

  	     FITS_read_col(infits, TFLOAT, tcol_in, 1, 1, nevents, NULL,
                           time_in, &anynull, &status);
	     FITS_read_col(infits, TSHORT, xcol_in, 1, 1, nevents, NULL,
		           x_in, &anynull, &status);
  	     FITS_read_col(infits, TSHORT, ycol_in, 1, 1, nevents, NULL,
		           y_in, &anynull, &status);
	     FITS_read_col(infits,  TBYTE, phacol_in, 1, 1, nevents, NULL,
		           pha_in, &anynull, &status);
	   if(procflg > 0) {
	    FITS_read_col(infits,  TFLOAT, dxcol_in, 1, 1, nevents, NULL,
		           dx, &anynull, &status);
	    FITS_read_col(infits,  TFLOAT, dycol_in, 1, 1, nevents, NULL,
		           dy, &anynull, &status);
	    FITS_read_col(infits,  TBYTE, dnflgcol_in, 1, 1, nevents, NULL,
		           dnflg, &anynull, &status);}

             for (j=0; j<nevents; j++){
                  time_in[j]+=(float) ((mjd_start-mjd_init)*86400.0);
	     }

             FITS_write_col(outfits, TFLOAT, tcol_out, jrow, 1,
			       nevents, time_in, &status);
	     FITS_write_col(outfits, TSHORT, xcol_out, jrow, 1,
			       nevents, x_in, &status);
	     FITS_write_col(outfits, TSHORT, ycol_out, jrow, 1,
			       nevents, y_in, &status);
	     FITS_write_col(outfits, TBYTE, phacol_out, jrow, 1,
			       nevents, pha_in, &status);
	   if (procflg > 0) {
	    FITS_write_col(outfits, TFLOAT, dxcol_out, jrow, 1,
			       nevents, dx, &status);
	    FITS_write_col(outfits, TFLOAT, dycol_out, jrow, 1,
			       nevents, dy, &status);
	    FITS_write_col(outfits, TBYTE, dnflgcol_out, jrow, 1,
			   nevents, dnflg, &status); }
  	     jrow+=nevents;

             free(time_in);
             free(x_in);
             free(y_in);
             free(pha_in);
             free(dx) ;
             free(dy) ;
             free(dnflg) ;


             /* Copy GTI from input to output */
             FITS_movabs_hdu(infits, 3, &hdutype, &status);
             FITS_read_key(infits,   TLONG, "NAXIS2", &ngti_in, NULL, &status);
             gti_start_in  = (float *) cf_malloc(sizeof(float) * ngti_in);
             gti_stop_in   = (float *) cf_malloc(sizeof(float) * ngti_in);
             FITS_get_colnum(infits, TRUE,  "START",  &startcol, &status);
             FITS_get_colnum(infits, TRUE,  "STOP",   &stopcol, &status);
             FITS_read_col(infits,   TFLOAT, startcol, 1, 1, ngti_in, 0,
  		  gti_start_in, &anynull, &status);
             FITS_read_col(infits,   TFLOAT, stopcol, 1, 1, ngti_in, 0,
		  gti_stop_in, &anynull, &status);
             for (k=0; k<ngti_in; k++) {
                  gti_start[ngti]=gti_start_in[k]+
                                  (float) ((mjd_start-mjd_init)*86400.0);;
                  gti_stop[ngti]=gti_stop_in[k]+
                                  (float) ((mjd_start-mjd_init)*86400.0);;
                  ngti++;
                  if (ngti >= MAXGTI) {
                       sprintf(buffer, "More than %d good time intervals.",MAXGTI) ;
                       cf_if_error(buffer);
		  }
	     }
             free(gti_start_in);
             free(gti_stop_in);

             FITS_close_file(infits, &status);
    }

    FITS_flush_file(outfits, &status);
    
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);
     jrow=jrow-1 ;
    FITS_update_key(outfits, TLONG, "NEVENTS", &jrow, 0, &status) ;
    FITS_update_key(outfits, TDOUBLE, "EXPTIME", &exptime, 0, &status);
    FITS_update_key(outfits, TDOUBLE, "EXPEND", &mjd_end, 0, &status);
    FITS_update_key(outfits, TLONG, "PLANTIME", &tplantime, 0, &status) ;

    if (expnight < 0)
	expnight = 0;
    FITS_update_key(outfits, TLONG, "EXPNIGHT", &expnight, 0, &status);

    if (trawtime - exptime > -1.)
        FITS_update_key(outfits, TDOUBLE, "RAWTIME", &trawtime, 0, &status);

    FITS_read_key(outfits, TSTRING, "DETECTOR", detstr, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"FE_B");
    FITS_read_key(outfits, TLONG, cntstr, &fe_b, NULL, &status);
     sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"DE_B");
    FITS_read_key(outfits, TLONG, cntstr, &de_b, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"SIC_B");
    FITS_read_key(outfits, TLONG, cntstr, &sic_b, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"LIF_B");
    FITS_read_key(outfits, TLONG, cntstr, &lif_b, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"AS_B");
    FITS_read_key(outfits, TLONG, cntstr, &as_b, NULL, &status);

    FITS_read_key(outfits, TLONG, "C1AAI_B", &ai1a_b, NULL, &status);
    FITS_read_key(outfits, TLONG, "C1BAI_B", &ai1b_b, NULL, &status);
    FITS_read_key(outfits, TLONG, "C2AAI_B", &ai2a_b, NULL, &status);
    FITS_read_key(outfits, TLONG, "C2BAI_B", &ai2b_b, NULL, &status);
    FITS_read_key(outfits, TDOUBLE, "CTIME_B", &ctm_b, NULL, &status);
 
    sprintf(comment, "New ID for combined data") ;
    strncpy(expid,"999               ",18) ;
     expid[18]='\0' ;
    FITS_update_key(outfits, TSTRING, "EXP_ID", expid, comment, &status);

    sprintf(comment, "Date and time file was created") ;
    stime[24]='\0' ;
    FITS_update_key(outfits, TSTRING, "DATE", stime, comment, &status) ;

    de_e = de_tot + de_b;
    fe_e = fe_tot + fe_b;
    lif_e = lif_tot + lif_b ;
    sic_e = sic_tot + sic_b ;
    as_e = as_tot + as_b ;
    ctm_e = ctm_tot + ctm_b ;
    ai1a_e = ai1a_tot + ai1a_b ;
    ai1b_e = ai1b_tot + ai1b_b ;
    ai2a_e = ai2a_tot + ai2a_b ;
    ai2b_e = ai2b_tot + ai2b_b ;

    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"DE_E");
    FITS_update_key(outfits, TLONG, cntstr, &de_e, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"FE_E");
    FITS_update_key(outfits, TLONG, cntstr, &fe_e, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"SIC_E");
    FITS_update_key(outfits, TLONG, cntstr, &sic_e, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%5.5s","C",detstr,"LIF_E");
    FITS_update_key(outfits, TLONG, cntstr, &lif_e, NULL, &status);
    sprintf(cntstr, "%1.1s%2.2s%4.4s","C",detstr,"AS_E");
    FITS_update_key(outfits, TLONG, cntstr, &as_e, NULL, &status);

    FITS_update_key(outfits, TLONG, "C1AAI_E", &ai1a_e, NULL, &status);
    FITS_update_key(outfits, TLONG, "C1BAI_E", &ai1b_e, NULL, &status);
    FITS_update_key(outfits, TLONG, "C2AAI_E", &ai2a_e, NULL, &status);
    FITS_update_key(outfits, TLONG, "C2BAI_E", &ai2b_e, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "CTIME_E", &ctm_e, NULL, &status);


    FITS_write_comment(outfits, "  ", &status);
    FITS_write_comment(outfits, "CF_TTAG_COMBINE", &status);
    FITS_write_comment(outfits, "  ", &status);
    FITS_write_comment(outfits, "This file is a concatenation of:", &status);
    for (i=2; i<=nfiles; i++) {
             FITS_write_comment(outfits, argv[i], &status);
    }
    

    FITS_open_file(&infits, argv[2], READONLY, &status);

    /* Move to the second extension which contains the goof time intervals. */
    FITS_movabs_hdu(infits, 3, &hdutype, &status);

    /* Create the second extension in the output file. */
    FITS_create_hdu(outfits, &status);

    /* Copy the header of the input file to the output file. */
    FITS_get_hdrpos(infits, &nkeys, &keynum, &status);
    for (i=1; i<=nkeys; i++) {
	FITS_read_record(infits, i, buffer, &status);
	FITS_write_record(outfits, buffer, &status);
    }

    FITS_get_colnum(outfits,TRUE,"START",&startcol,&status);
    FITS_get_colnum(outfits,TRUE,"STOP",&stopcol,&status);

    FITS_write_col(outfits, TFLOAT, startcol, 1, 1, ngti, gti_start, &status);
    FITS_write_col(outfits, TFLOAT, stopcol, 1, 1, ngti, gti_stop, &status);

    free(gti_start);
    free(gti_stop);

    FITS_close_file(infits, &status);
    FITS_close_file(outfits, &status);

    if (fpa_split)
       cf_if_warning("FPA motion detected. Resolution may be compromised.\n");

    cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Finished execution.");
    return 0;
}
