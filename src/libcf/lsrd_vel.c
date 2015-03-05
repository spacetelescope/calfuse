/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    double lsrd_vel(double ra, double dec)
 *
 * Description: Computes the Sun's velocity in the direction 
 *              of ra,dec to convert to a LSR.  In this case we
 *              are using the dynamical solar motion of 
 *              16.6 km/s towards l=53, b=+25.
 *
 * Arguments:   double ra             (deg) J2000.0 right ascension
 *              double dec            (deg) J2000.0 declination 
 *
 * Returns:     double                (km/s) velocity in direction of source
 *                    
 * Calls:       slaRvlsrd             rvslrd.f
 *
 * History:     03/05/98  E. Murphy   Begin work.
 *              03/05/98  E. Murphy   Initial version working
 *              05/31/00        peb   Implemented cfortran.h calls for slalib
 *                                    functions.
 *              06/15/01  V. Dixon    Multiply return value by -1.0 to agree
 *				      with astronomical convention.
 *              12/18/03    bjg       Change calfusettag.h to calfuse.h
 *
 * References:  This routine makes use of the Starlink set of astronomical
 *              subroutines (SLALIB).  More information can be found at
 *              http://star-www.rl.ac.uk.
 ****************************************************************************/
 
#include <stdio.h>

#ifdef CFORTRAN
#include "cfortran.h"
PROTOCCALLSFFUN2(FLOAT, SLA_RVLSRD, sla_rvlsrd, FLOAT, FLOAT)
#define slaRvlsrd(R2000, D2000) \
     CCALLSFFUN2(SLA_RVLSRD, sla_rvlsrd, FLOAT, FLOAT, R2000, D2000)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"
 
double lsrd_vel(double ra, double dec)
{
    return -1.0 * slaRvlsrd((float) ra*RADIAN, (float) dec*RADIAN);
}
