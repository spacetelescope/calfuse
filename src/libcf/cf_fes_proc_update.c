/******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 ******************************************************************************
 *
 * Synopsis:    cf_fes_proc_update(fitsfile *fptr, char *prgm_id, 
 *                                 char *key_value)
 *
 * Description: cf_fes_proc_update will update the FITS header keyword 
 *              which corresponds to the given prgm_id in the file 
 *              fptr to the value given in key_value.
 *
 * Arguments:   fitsfile    *fptr	Pointer to input file
 *              char        *prgm_id 	Procedure name
 *              char        *key_value 	Updated procedure status
 *
 * History:     07/21/98        emm     Begin work.
 *              08/19/2004 1.1  wvd     Move to v3.0, add CF_VER_NUM,
 *                                      change cf_errmsg to cf_if_error.
 *					Delete fits_print_err.
 *		09/07/2007 1.2	bot	Removed unused CF_VER_NUM and CF_PRGM_ID
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "calfuse.h"

int cf_fes_proc_update(fitsfile *fptr, char *prgm_id, char *key_value) 
{
    int i=0, status=0;
    /*
     *  The calfuse.h file contains the definitions of fes_keyword_tab
     *  NUM_FES_PROC_STEPS, and FES_CALIBRATION_STEP_KEYS.
     */
    struct fes_keyword_tab keytab[NUM_FES_PROC_STEPS]=FES_CALIBRATION_STEP_KEYS;
    /*
     *  Find the keyword associated with prgm_id by looping 
     *  through keytab[i].proc
     */
    while ((strncmp(keytab[i].proc,prgm_id, strlen(keytab[i].proc)) != 0) &&
	   (i < NUM_FES_PROC_STEPS)) i++;
    if (i < NUM_FES_PROC_STEPS) {
	/*
	 *  We found a match to prgm_id, so change the associated
	 *  keyword in the header.
	 */
	fits_modify_key_str(fptr,keytab[i].name,key_value,
			    "&",&status);
	if (status) {
	    cf_if_error("Error updating keyword");
	    return status;
	}
    }
    else {
	/*
	 *  The given prgm_id did not match any of the known
	 *  keytab[i].hist_proc, so return 1.
	 */
	cf_if_error("Program ID does not apply to this"
		  "type of file, could not update");
	fprintf(stderr,"Failed to update program id:%20.20s\n",
		    prgm_id);
	return 1;
    }
    return status;
}
