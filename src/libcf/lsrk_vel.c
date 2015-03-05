/****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 ****************************************************************************
 *
 * Synopsis:    double lsrk_vel(double ra, double dec)
 *
 * Description: Computes the Sun's velocity in the direction 
 *              of ra,dec to convert to a LSR.  In this case we
 *              are using standard solar motion of 
 *              20 km/s towards RA 18h Dec +30d (1900).
 *
 * Arguments:   double ra             (deg) J2000.0 right ascension
 *              double dec            (deg) J2000.0 declination 
 *
 * Returns:     double                (km/s) velocity in direction of source
 *                    
 * Calls:       slaRvlsrk             rvlsrk.f
 *
 * History:     03/05/98  E. Murphy   Begin work.
 *              03/05/98  E. Murphy   Initial version working
 *              04/13/99  E. Murphy   Clean up a bit
 *              05/31/00        peb   Implemented cfortran.h calls for slalib
 *                                    functions.
 *              06/15/01  V. Dixon    Multiply return value by -1.0 to agree
 *                                    with astronomical convention.
 *              12/18/03    bjg       Change calfusettag.h to calfuse.h
 *
 * References:  This routine makes use of the Starlink set of astronomical
 *              subroutines (SLALIB).  More information can be found at
 *              http://star-www.rl.ac.uk.
 ****************************************************************************/
 
#include <stdio.h>

#ifdef CFORTRAN
#include "cfortran.h"
PROTOCCALLSFFUN2(FLOAT, SLA_RVLSRK, sla_rvlsrk, FLOAT, FLOAT)
#define slaRvlsrk(R2000, D2000) \
     CCALLSFFUN2(SLA_RVLSRK, sla_rvlsrK, FLOAT, FLOAT, R2000, D2000)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"

double lsrk_vel(double ra, double dec)
{
    return -1.0 * slaRvlsrk((float) ra*RADIAN, (float) dec*RADIAN);
}
