/******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 ******************************************************************************
 *
 * Synopsis:    cf_fes_proc_check(fitsfile *fptr, char *prog_id)
 *
 * Description: cf_fes_proc_check will determine if a given calibration
 *              step is to be performed on the data.  It will also
 *              determine whether all previous steps have been
 *              completed successfully.
 *
 * Arguments:   fitsfile    *fptr	Pointer to input file
 *              char        *prog_id 	Procedure name
 *
 * History:     06/21/98        emm     Begin work.
 *		08/19/2004 1.1	wvd	Move to v3.0, add CF_VER_NUM,
 *					change cf_errmsg to cf_if_error.
 *		09/07/2007 1.2	bot	Removed unused CF_VER_NUM and CF_PRGM_ID
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "calfuse.h"
#define  MAXCHARS 120

int cf_fes_proc_check(fitsfile *fptr, char *prog_id) 
{
    int   i,j, status=0;
    char  comment[FLEN_CARD], key_value[FLEN_CARD];
    char  complete[19], skipped[19], perform[19];

    /*
     *  The calfuse.h file contains the definitions of fes_keyword_tab
     *  NUM_PROC_STEPS, and CALIBRATION_STEP_KEYS.
     */
    struct fes_keyword_tab keytab[NUM_FES_PROC_STEPS]=FES_CALIBRATION_STEP_KEYS;
    status=0;

    strncpy(complete,"COMPLETE          ",19);
    strncpy(skipped, "SKIPPED           ",19);
    strncpy(perform, "PERFORM           ",19);

    j=0;
    /*
     *  First, determine if this procedure is even supposed to be
     *  run on this data.
     */
    while ((strncmp(keytab[j].proc,prog_id, strlen(keytab[j].proc)) != 0) &&
	   (j < NUM_FES_PROC_STEPS)) j++;
    if ((j >= NUM_FES_PROC_STEPS) || 
	(strncmp(keytab[j].value,perform,7))) {
	cf_if_error("Processing step does not "
		  "need to be run on this type of data.");
	fprintf(stderr,"     Step %18.18s does not need "
		"to be run.\n",prog_id);
	return 1;
    }

    status=0;
    fits_read_key_str(fptr, keytab[j].name, key_value, comment, 
		      &status);

    /* Now check to see if the step has already been completed. */
    if (strncmp(key_value,complete,7)==0) {
	cf_if_error("Processing step has already been completed.\n");
	fprintf(stderr,"     Step %18.18s does not need to be run.\n",
		prog_id);
	return 1;
    }

    /*  Now determine if the previous programs are all complete. */
    for (i=0; i<j; i++) {
	fits_read_key_str(fptr, keytab[i].name, key_value, comment, 
			  &status);
	if (strncmp(keytab[i].value,perform,7)==0) {
	    if (strncmp(key_value,complete,8) &&
		strncmp(key_value,skipped,7)) {
		cf_if_error("Processing step not completed.");
		fprintf(stderr,"     Step %8.8s has value "
			"%8.8s.\n",keytab[i].name,key_value);
		fprintf(stderr,"     It should be either "
			"COMPLETE or SKIPPED before running "
			"%18.18s\n",prog_id);
		return 1;
	    } /* End if check on key_value */
	} /* End check on keytab[i].value */
    } /* End for */
   
    /* If we got to here, it must be OK */
    return 0;
}
