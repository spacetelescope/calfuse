/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_cal_file(char *calfile)
 *
 * Description: Returns a pointer to a string with the full calibration 
 *              path (from the CF_CALDIR environment variable) prepended 
 *              to the given file name.
 *
 * Arguments:   char    *calfile        Name of the calibration file.
 *
 * Returns:     A pointer to a string which contains the full path filename.
 *
 * History:     02/26/99       emurphy  Begin and finish work.
 *              03/01/99       emurphy  Added some error checking
 *              03/02/99       emurphy  Added use of cf_if_warning
 *              04/04/99       emurphy  removed cf_error_init
 *              04/08/99       barrett  replaced fitsio.h with calfuse.h
 *              05/18/99       emurphy  Added checks on strlen and ending "/"
 *              06/04/99          peb   Removed 160 character string limit
 *                                      Added cf_malloc function.
 *              08/05/99       emurphy  Created cf_populate_file, to support
 *                                      cf_cal_file and cf_hist_file.c
 *              08/25/99       emurphy  Added cf_parm_file for parameter files.
 *              12/18/03       bjg      Change calfusettag.h to calfuse.h
 *
 ****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "calfuse.h"

static char
*cf_populate_file(char *calfile, char *enviro_name)
{
    char *enviro, *filen;
    int strl, strt;

    enviro = getenv(enviro_name);

    if (enviro) {
        strl = strlen(enviro);
	strt = strl + strlen(calfile) + 16;

	filen = cf_malloc(strt);
        strcpy(filen, enviro);

        /* Check to see if last character is a "/" */
        if (strncmp(filen+strl-1, "/", 1)!=0)
             strcat(filen, "/");

        strcat(filen, calfile);
        return filen;
    }
    else {
        printf("Environment variable %-20.20s undefined.",enviro_name);
        cf_if_warning("Environment variable undefined");
        return calfile;
    }
}

char *cf_cal_file(char *calfile)
{

     return cf_populate_file(calfile, "CF_CALDIR");

}

char *cf_hist_file(char *calfile)
{

     return cf_populate_file(calfile, "CF_HISTDIR");

}

char *cf_parm_file(char *calfile)
{

     return cf_populate_file(calfile, "CF_PARMDIR");

}
