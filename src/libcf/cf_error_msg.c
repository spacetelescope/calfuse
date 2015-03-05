/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_error_init(char *prgmid, char *vernum, FILE *errfile)
 *
 * Description: Writes a timestamp, the calling program name, and an
 *		error message to stdout for a cfitsio routine.
 *
 * Arguments:   char    *progid         Name of the calling program
 *              char    *vernum         Version of the calling program
 *                                      (should be version control number)
 *              FILE    *errfile        Output file name (default is stderr)
 *
 * Returns:     none
 *
 * History:     02/22/99        peb     Some simple functions to make error
 *                                      reporting easier.  Deficiencies in
 *                                      cfitsio make printing comprehensive
 *                                      messages difficult.
 *              03/15/99        emm     Added cf_if_memory_error.
 *              04/30/99        peb     Added cf_malloc and cf_calloc.
 *              06/07/99        peb     Added reporting of version number 
 *              09/01/99        emm     Added fflush to see if it helps print
 *                                      out error messages on failure.
 *              11/08/99        emm     Added !!!!!!!! to error reporting to
 *                                      better delimit the error message and
 *                                      to make them more obvious.
 *              04/23/01  1.10  wvd     Make cf_malloc & cf_calloc print error
 *					message returned by malloc & calloc
 *              03/07/03  1.2   peb     Added cf_verbose and external variable
 *                                      verbose_level
 *              03/13/03  1.3   peb     Enabled cf_verbose, cf_if_warning, and
 *                                      cf_if_error to accept a variable
 *                                      number of parameters.  Changed the
 *                                      output format to print month, day, time
 *                                      hostname, function name and message.
 *              03/25/03  1.4   peb     Changed output to print time, program,
 *                                      and version number for warning and
 *                                      error messages.
 *                                      Changed output to use verbose_level
 *                                      to modify output format:
 *                                      level==1 prints program, version, and
 *                                      message.
 *                                      level>=2 prints only message.
 *              10/30/03  1.5   peb     Replaced cftime function with strftime
 *                                      for UNIX compatibility.
 *                                      Removed unnecessary header files,
 *                                      malloc.h and unistd.h
 *              04/07/04  1.6   bjg     fflush stdout and error_file
 *
 ****************************************************************************/

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fitsio.h"

#define  N_STYM  40

int   verbose_level=0;

static char *program_name=NULL;
static char *version_no=NULL;
static FILE *error_file=NULL;

void cf_error_init(const char *progid, const char *vernum, FILE *errfile)
{
    static char unknown[] = "unknown file";
    static char version[] = "1.6";
    if (progid) {
        program_name = malloc(strlen(progid)+1);
        strcpy(program_name, progid);
    }
    else {
        program_name = malloc(strlen(unknown)+1);
        strcpy(program_name, unknown);
    }
    if (vernum) {
        version_no = malloc(strlen(vernum)+1);
	strcpy(version_no, vernum);
    }
    else {
        version_no = malloc(strlen(version)+1);
	strcpy(version_no, version);
    }
    if (errfile)
        error_file = errfile;
    else
        error_file = stderr;
}

void cf_verbose(int level, const char *format, ...)
{
    va_list args;
    va_start(args, format);

    if (verbose_level >= level) {
	if (level == 1)
	    printf("   %s-%s: ", program_name, version_no);
	else if (level >= 2)
	    printf("\t");
	vprintf(format, args);
	printf("\n");
	fflush(stdout);
    }

    va_end(args);
}

void cf_if_warning(char *format, ...)
{
    char stym[N_STYM]={'\0'};
    time_t tym;
    va_list args;

    va_start(args, format);
    time(&tym);
    strftime(stym, N_STYM, "%Y %b %e %T", localtime(&tym));

    fprintf(error_file, "%s %s-%s: WARNING - ",
	    stym, program_name, version_no);
    vfprintf(error_file, format, args);
    fprintf(error_file, "\n");
    fflush(error_file);
    va_end(args);
}

void cf_if_error(char *format, ...)
{
    char stym[N_STYM]={'\0'};
    time_t tym;
    va_list args;

    va_start(args, format);
    time(&tym);
    strftime(stym, N_STYM, "%Y %b %e %T", localtime(&tym));

    fprintf(error_file, "%s %s-%s: ERROR - ",
	    stym, program_name, version_no);
    vfprintf(error_file, format, args);
    fprintf(error_file, "\n");

    fflush(error_file);
    va_end(args);
    exit(1);
}

void cf_if_fits_warning(int status)
{
    if (status) {
	char stym[N_STYM]={'\0'};
	time_t tym;

	time(&tym);
	strftime(stym, N_STYM, "%Y %b %e %T", localtime(&tym));

	fprintf(error_file, "%s %s-%s: FITS WARNING - ",
		stym, program_name, version_no);
        fits_report_error(error_file, status);
        fflush(error_file);
    }
}

void cf_if_fits_error(int status)
{
    if (status) {
	char stym[N_STYM]={'\0'};
	time_t tym;

	time(&tym);
	strftime(stym, N_STYM, "%Y %b %e %T", localtime(&tym));

	fprintf(error_file, "%s %s-%s: FITS ERROR - ",
		stym, program_name, version_no);
        fits_report_error(error_file, status);
        fflush(error_file);
        exit(1);
    }
}

void *cf_malloc(size_t size)
{
    void *ptr;

    if (!(ptr = malloc(size))) {
	char stym[N_STYM]={'\0'};
	time_t tym;

	time(&tym);
	strftime(stym, N_STYM, "%Y %b %e %T", localtime(&tym));

	fprintf(error_file, "%s %s-%s: ERROR - %s\n",
		stym, "cf_malloc", "1.0", strerror(errno));

        fflush(error_file);
	exit(1);
	/*
	fprintf(error_file,
		"     \n     malloc: Attemping to allocate %d bytes", size);
	*/
    }
    return ptr;
}

void *cf_calloc(size_t nelem, size_t elsize)
{
    void *ptr;

    if (!(ptr = calloc(nelem, elsize))) {
	char stym[N_STYM]={'\0'};
	time_t tym;

	time(&tym);
	strftime(stym, N_STYM, "%Y %b %e %T", localtime(&tym));

	fprintf(error_file, "%s %s-%s: ERROR - %s\n",
		stym, "cf_calloc", "1.0", strerror(errno));

        fflush(error_file);
	exit(1);
	/*
	fprintf(error_file,
		"     \n     calloc: Attemping to allocate %d bytes",
				nelem * elsize);
	*/
    }
    return ptr;
}

void cf_if_memory_error(int status)
{
    time_t tym;

    if (!status) {
        time(&tym);
        fprintf(error_file, "Memory allocation error in %s-%s: %s",
                program_name, version_no, ctime(&tym));

        /* Print out cfitsio error messages and continue */

        fits_report_error(error_file, status);
	fflush(error_file);
    }
}

/*
int main()
{
    int i;
    cf_error_init("test program", "1.0", stdout);
    cf_if_warning("Non-fatal error occurred");

    for (i=0; i < 10; i++)
        cf_if_fits_warning(i);
    exit(0);
}
*/
