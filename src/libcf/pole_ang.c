/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    double pole_ang(double *pos, double *vel, double ra, double dec)
 *
 * Description: Computes the position of the orbit pole from the given state 
 *		6-vector.  Returns the angle between the pole and the given
 *		RA and DEC.
 *
 * Arguments:   double *pos           (km) 3-vector position of satellite
 *              double *vel           (km/s) 3-vector velocity of satellite
 *              double ra,dec         (deg) position of source
 *
 * Returns:     double                (deg) angle from orbit pole to source
 *
 * History:     07/06/06  wvd   1.1   Adapted from space_vel and sun_ang
 *		04/08/07  wvd	1.2   Cast SGP4_getPole as a void.
 *
 * References:  This routine makes use of the Starlink set of astronomical
 *              subroutines (SLALIB).  More information can be found at
 *              http://star-www.rl.ac.uk.
 *****************************************************************************/
 
#include <stdio.h>

#ifdef CFORTRAN
#include "cfortran.h"
PROTOCCALLSFSUB3(SLA_DCC2S, sla_dcc2s, DOUBLEV, PDOUBLE, PDOUBLE)
#define slaDcc2s(V, A, B) \
     CCALLSFSUB3(SLA_DCC2S, sla_dcc2s, DOUBLEV, PDOUBLE, PDOUBLE, V, A, B)
PROTOCCALLSFFUN1(DOUBLE, SLA_DRANRM, sla_dranrm, DOUBLE)
#define slaDranrm(ANGLE) \
     CCALLSFFUN1(SLA_DRANRM, sla_dranrm, DOUBLE, ANGLE)
PROTOCCALLSFFUN4(DOUBLE, SLA_DSEP, sla_dsep, DOUBLE, DOUBLE, DOUBLE, DOUBLE)
#define slaDsep(A1, B1, A2, B2) \
     CCALLSFFUN4(SLA_DSEP, sla_dsep, DOUBLE, DOUBLE, DOUBLE, DOUBLE, \
                 A1, B1, A2, B2)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"
 
double pole_ang(double pos[3], double vel[3], double ra, double dec)
{
    double pole[3], ra_pole, dec_pole;

    /* Get 3-vector containing the orbit pole. */
    (void) SGP4_getPole(pos, vel, pole);

    /* Convert 3-vector into J2000.0 RA and DEC */
#ifdef CFORTRAN
    slaDcc2s(pole, ra_pole, dec_pole);
#else
    slaDcc2s(pole, &ra_pole, &dec_pole);
#endif

    /*  Compute the angular separation of RA, DEC and RA_POLE, DEC_POLE */
    return slaDsep(ra*RADIAN, dec*RADIAN, slaDranrm(ra_pole), dec_pole)/RADIAN;
}
