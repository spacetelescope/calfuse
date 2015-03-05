
/*******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *******************************************************************************
 *
 * Synopsis:    cf_fes_init(fitsfile *fptr)
 *
 * Description: cf_fes_init performs two functions.  First, populates the 
 *              data processing step keywords based upon the type of 
 *              data (FESA or FESB).  Second, it populates the calibration
 *              file keywords based on the time and date of observation,
 *              the detector segment, and the interpolation status of the
 *              calibration file.
 *
 * Arguments:   fitsfile    *fptr	Pointer to input file
 *
 * History:     05/21/98        emm     Begin work.
 *              06/17/98        emm     Modified keywords to reflect updated
 *                                      calfuse design
 *              07/15/98        mlr     copied from emm
 *		08/23/04	wvd	Change from fits_modify_key_str to
 *					FITS_update_key throughout.
 ******************************************************************************/

#include <stdio.h>
#include <string.h>
#include "calfuse.h"

#define  MAXCHARS 120
#define  MASTER_CAL_FILE "master_fes_calib_file.dat"
#define  CF_PRGM_ID      "cf_fes_init"
#define  CF_VER_NUM      "1.4"

int cf_fes_init(fitsfile *fptr) 
{
int   i,j,k,l,flg;
char  comment[FLEN_CARD], instmode[FLEN_CARD], detector[FLEN_CARD];
char  fes[2];
char  keyword_in[9], segment_in[3];
char  filename_in[19]="                  \0";
char  keyword_out1[9], keyword_out2[9];
char  complete[19], error[19], blank[19];
char  linin[120];
int   status, interp_in;
float expstart, aftermjd_in;

struct fes_keyword_tab keytab[NUM_FES_PROC_STEPS]=FES_CALIBRATION_STEP_KEYS;
struct cal_file_tab calkey[NUM_FES_CAL_KEYS]=FES_CALIBRATION_FILE_KEYS;
FILE  *fpin;

cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Begin initializing header");
cf_error_init(CF_PRGM_ID, CF_VER_NUM, stderr);

status=0;

   strncpy(complete,"COMPLETE          \0",19);
   strncpy(error,   "ERROR             \0",19);
   strncpy(blank,   "                  \0",19);

/*  Populate data processing keywords. */

   for (i=0; i<NUM_FES_PROC_STEPS; i++)
        FITS_update_key(fptr,TSTRING,keytab[i].name,
                   keytab[i].value,NULL,&status);


/*  Populate calibration file keywords based on the date of the 
 *  observation and the segment/channel number. */

   FITS_read_key(fptr, TFLOAT, "EXPSTART", &expstart, comment, &status);
   FITS_read_key(fptr, TSTRING, "FES_ID", instmode, comment, &status);
   sscanf(instmode,"%*3c %1c",fes);

/*  Open the Master calibration database file. */

   fpin=NULL;
   fpin=fopen(cf_cal_file(MASTER_CAL_FILE),"r");
   if (fpin == NULL)
             cf_if_warning("Master calibration database file not found");

   /*  Cycle through the master cal file until we run out of lines. */

   /*  Here's how this section works.  We need to keep track of 3 calibration
    *  files for each keyword: the cal file which immediately preceedes the 
    *  date of observation, the cal file which immediately preceedes the cal
    *  file which immediately preceedes the observation, and the cal file 
    *  immediately after the observation.  These are stored in strut calkey
    *  in:
    *     calkey[i].filenames[0] => second cal file preceeding observation
    *     calkey[i].filenames[1] => cal file immediately preceeding observation
    *     calkey[i].filenames[2] => cal file immediately after observation
    *
    *  Here is the plan.  Cycle through the master cal file, reading each
    *  line.  Determine the keyword for each line, and use that to determine
    *  the appropriate index [i] for calkey.  Then, if the detector segment 
    *  is correct, determine where the new file fits into the structure above.
    *  If it is more recent than one or both of the preceeding cal files, 
    *  but is older than the observation (expstart) replace either 
    *  filenames[0] or filenames[1].  If filenames[1] is replaced, shift
    *  the old value down into filenames[0].  If the file is more recent 
    *  than the observation, and closer to the observation than the file 
    *  in filenames[2], replace filenames[2].  The values 
    *  calkey[i].aftermjd[] and calkey[i].interp[] store the date and 
    *  interpolation status of the best files.  This organization allows
    *  the master calibration database file to be written in any order.
    */

   while (fgets(linin, MAXCHARS, fpin) != NULL) {

        /* Check for comment lines */

       if ((linin[0] != '#') && (linin[0] != '\n')) {
              sscanf(linin, "%4c %1c %12c %9f %1d",keyword_in,
                           segment_in,filename_in,&aftermjd_in,&interp_in);

	      /* Determine which keyword we just read in. */
              i=0;
              while ((strncmp(calkey[i].name,keyword_in,4) != 0) && 
                                      (i < NUM_FES_CAL_KEYS)) i++;

              if (i >= NUM_FES_CAL_KEYS) {
                   cf_if_warning("Unrecognized keyword in master calibration file");
              } else {
                   if ((strncmp(segment_in,"  ",1) == 0) || 
                       (strncmp(segment_in,fes,1) == 0)) {
                        j=0;
                        if ((aftermjd_in < expstart) && 
                            (aftermjd_in > calkey[i].aftermjd[0]) &&
                            (aftermjd_in > calkey[i].aftermjd[1])) {
                                  strncpy(calkey[i].filenames[0],
                                          calkey[i].filenames[1],18);
                                  calkey[i].aftermjd[0]=calkey[i].aftermjd[1];
                                  calkey[i].interp[0]=calkey[i].interp[1];
                                  strncpy(calkey[i].filenames[1],
                                          filename_in,18);
                                  calkey[i].aftermjd[1]=aftermjd_in;
                                  calkey[i].interp[1]=interp_in;
                         } else if ((aftermjd_in < expstart) && 
                            (aftermjd_in > calkey[i].aftermjd[0])) {
                                  strncpy(calkey[i].filenames[0],
                                          filename_in,18);
                                  calkey[i].aftermjd[0]=aftermjd_in;
                                  calkey[i].interp[0]=interp_in;
			 }
     
                        if ((aftermjd_in > expstart) && 
                            (aftermjd_in < calkey[i].aftermjd[2])) {
                                  strncpy(calkey[i].filenames[2],
                                          filename_in,18);
                                  calkey[i].aftermjd[2]=aftermjd_in;
                                  calkey[i].interp[2]=interp_in;
			}
		   } /* Endif segment match */
	      } /* end else unrecognized keyword */		   
	}  /* end while linin not comment */
   } /* end while(linin != NULL) */

   fclose(fpin);

   /*  OK, now that we have identified the closest two files before
    *  the observation, and the file after the observation, we can
    *  decide what to do with them.  First, the easy part:  if 
    *  calkey[i].numfiles=1, no interpolation is possible, so we
    *  write out calkey[i].filenames[1], which is the file immediately
    *  preceeding the observation.  If calkey[i].numfiles=2, interpolation
    *  is possible, and we will have to write two keywords into the header.
    *  The decision as to which two keywords to write follows these rules:
    *  Rule 1: If an observation occurs between two files which can be 
    *          interpolated, then an interpolation must be done. 
    *  Rule 2: If the calibration file immediately preceeding an 
    *          observation is to be used in a stepwise manner, then said
    *          file shall be used in stepwise manner. 
    *  Rule 3: If the two calibration files immediately preceeding an 
    *          observation can be interpolated, but the following file 
    *          cannot, then the two previous files will be extrapolated 
    *          forward in time. 
    *  Rule 4: If the closest preceeding calibration file has an interpolation
    *          flag, but the file preceeding it has a stepwise flag, the the 
    *          closest preceeding file will be used in a sterpwise manner.
    *  Rule 5: Never interpolate backward in time.
    *
    *  In fact, there are only 8 possible cases:
    *
    *  Remember, the files were sorted such that the observation 
    *  always occurs between [1] and [2].  Also remember that 
    *  interp[i]=0 => stepwise):
    *
    *    Case               Action
    *
    *  interp[0] = 0
    *  interp[1] = 0    filenames[1] will be used in a stepwise manner
    *  interp[2] = 0
    *
    *  interp[0] = 0
    *  interp[1] = 0    filenames[1] will be used in a stepwise manner
    *  interp[2] = 1
    *
    *  interp[0] = 0
    *  interp[1] = 1    filenames[1] will be used in a stepwise manner
    *  interp[2] = 0
    *
    *  interp[0] = 1
    *  interp[1] = 0    filenames[1] will be used in a stepwise manner
    *  interp[2] = 0
    *
    *  interp[0] = 1
    *  interp[1] = 0    filenames[1] will be used in a stepwise manner
    *  interp[2] = 1
    *
    *  interp[0] = 0
    *  interp[1] = 1    filenames[1] and filenames[2] will be interpolated
    *  interp[2] = 1
    *
    *  interp[0] = 1
    *  interp[1] = 1    filenames[0] and filenames[1] will be extrapolated
    *  interp[2] = 0
    *
    *  interp[0] = 1
    *  interp[1] = 1    filenames[1] and filenames[2] will be interpolated
    *  interp[2] = 1
    *
    */

   /*  Write the calibration files into the header keywords. */
   /*  Cycle through the keywords. */
   for (i=0; i<NUM_FES_CAL_KEYS; i++) {

     if (calkey[i].numfiles == 1) {

          /*  This is the easy case, no interpolation is possible. */

          sprintf(keyword_out1,"%4.4s%1.1s%3.3s",calkey[i].name,"_",
                    calkey[i].extension);
          FITS_update_key(fptr,TSTRING,keyword_out1,calkey[i].filenames[1],
                    NULL,&status); 
     } else {

          /*  OK, interpolation is possible, now decide which two files
           *  to write out. */

          /* Create the two output keywords */
          sprintf(keyword_out1,"%4.4s%1.1s%3.3s",calkey[i].name,"1",
                    calkey[i].extension);
          sprintf(keyword_out2,"%4.4s%1.1s%3.3s",calkey[i].name,"2",
                    calkey[i].extension);

       if (calkey[i].interp[1] == 0) {

	 /*  Easy, preceeding file is a 0, so write filenames[1] to 
          *  both keywords.  Takes care of cases 000, 001, 100, 101.
          */
          FITS_update_key(fptr,TSTRING,keyword_out1,calkey[i].filenames[1],
                    NULL,&status); 
          FITS_update_key(fptr,TSTRING,keyword_out2,calkey[i].filenames[1],
                    NULL,&status); 
       } else if ((calkey[i].interp[1] == 1) && (calkey[i].interp[2] == 1)) {
	 /* Case 011 and 111 */
          FITS_update_key(fptr,TSTRING,keyword_out1,calkey[i].filenames[1],
                    NULL,&status); 
          FITS_update_key(fptr,TSTRING,keyword_out2,calkey[i].filenames[2],
                    NULL,&status); 
       } else if ((calkey[i].interp[0] == 1) && (calkey[i].interp[1] == 1) 
		  && (calkey[i].interp[2] == 0)){
	 /* Case 110 */
          FITS_update_key(fptr,TSTRING,keyword_out1,calkey[i].filenames[0],
                    NULL,&status); 
          FITS_update_key(fptr,TSTRING,keyword_out2,calkey[i].filenames[1],
                    NULL,&status); 

       } else if ((calkey[i].interp[0] == 0) && (calkey[i].interp[1] == 1) 
		  && (calkey[i].interp[2] == 0)){
	  /* Case 010 */
          FITS_update_key(fptr,TSTRING,keyword_out1,calkey[i].filenames[1],
                    NULL,&status); 
          FITS_update_key(fptr,TSTRING,keyword_out2,calkey[i].filenames[1],
                    NULL,&status); 

       } else cf_if_warning("The interpolation/stepwise logic has failed.  I am confused.");

     } /* end else calkey.numfiles == 1 */

   } /* end for i=0; i<NUMCALKEYS */


FITS_update_key(fptr,TSTRING,"INIT_FES",complete,NULL,&status);

cf_timestamp(CF_PRGM_ID, CF_VER_NUM, "Done initializing header");    

return status;
}


