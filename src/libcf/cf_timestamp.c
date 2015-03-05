/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_timestamp(const char *prgmid, const char *vernum, char *msg)
 *
 * Description: Writes a timestamp and the calling program name to stdout.
 *
 * Arguments:   char    *prgmid         Name of the calling program
 *              char	*msg		Optional additional message string
 *
 * Returns:     none
 *
 * History:     01/04/91        mlr     Original "warn_str.c" for libhut
 *		03/03/98	gak	Adapted for placing timestamps in
 *					FUSE pipeline programs.
 *              06/07/99        peb     Added reporting of version number.
 *              10/27/99        emm     Added fflush(stdout)
 *              04/23/01        wvd     Change declaration of tim from
 *                                      long to time_t
 *              03/19/03        peb     Made this function a wrapper for
 *                                      cf_verbose.
 *              03/25/03        peb     Changed output to print time, program,
 *                                      version number and to use verbose_level
 *                                      when printing output
 *                                      level>=1 prints "Finished processing"
 *                                      level>=2 prints "Begin processing"
 *              09/15/03  v1.5  wvd     To reduce size of the trailer file,
 *					if verbose level = 1, print "Finished"
 *					only for top-level routines.
 *              10/17/03  v1.6  wvd     Use strcmp, rather than strncmp, to
 *					set this_is_a_top_level_routine.
 *              10/30/03  1.7   peb     Replaced cftime function with strftime
 *                                      for UNIX compatibility.
 *
 ****************************************************************************/

#include <string.h>
#include <time.h>
#include "calfuse.h"

#define  N_STYM  80

void cf_timestamp(const char *prgmid, const char *vernum, char *msg)
{
    char stym[N_STYM]={'\0'};
    char *top_level_routines[]=TOP_LEVEL_ROUTINES;
    int  i, this_is_a_top_level_routine = FALSE;
    time_t tym;

    time(&tym);
    strftime(stym, N_STYM, "%Y %b %e %T", localtime(&tym));

    /* Compare program ID with list of top-level routines. */
    for (i = 0; i < NTOP_LEVEL_ROUTINES; i++) {
	if (!strcmp(prgmid, top_level_routines[i])) {
	    this_is_a_top_level_routine = TRUE;
	    break;
	}
    }

    /*
     *  The following test should get most "Finished" or "Done" type of
     *  timestamps.  For those not caught if might be best to change them to
     *  one of these two.
     */
    if (((verbose_level == 1 && this_is_a_top_level_routine) ||
	(verbose_level > 1)) &&
	(strncmp(msg, "Finish", 6) == 0 || strncmp(msg, "Done", 4) == 0))
	printf("%s %s-%s: %s\n", stym, prgmid, vernum, "Finished processing");
    else if (verbose_level >= 2 &&
	     (strncmp(msg, "Start", 5) == 0 || strncmp(msg, "Begin", 5) == 0))
	printf("%s %s-%s: %s\n", stym, prgmid, vernum, "Begin processing");

    fflush(stdout);
}
