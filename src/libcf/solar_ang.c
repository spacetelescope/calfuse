/******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 ******************************************************************************
 *
 * Synopsis:    double solar_ang(double mjdate, double ra, double dec)
 *
 * Description: Computes the angle between the Sun and J2000.0 RA and DEC
 *              for a given date.
 *
 * Arguments:   double mjdate         Modified Julian date of observation
 *              double ra             (deg) J2000.0 right ascension
 *              double dec            (deg) J2000.0 declination 
 *
 * Returns:     double                (deg) angle between Sun and RA,DEC
 *                    
 * Calls:       slaDcc2s              dcc2s.f
 *              slaDsep               dsep.f
 *              slaDranrm             dranrm.f
 *              slaEvp                evp.f
 *
 * History:     03/10/98  E. Murphy   Begin work.
 *              03/10/98  E. Murphy   Initial version working
 *              05/31/00        peb   Implemented cfortran.h calls for slalib
 *                                    functions.
 *              12/18/03    bjg       Change calfusettag.h to calfuse.h
 *		07/21/04  1.4   wvd   Comment out variable vect; unused.
 *
 * References:  This routine makes use of the Starlink set of astronomical
 *              subroutines (SLALIB).  More information can be found at
 *              http://star-www.rl.ac.uk.
 *****************************************************************************/
 
#include <stdio.h>
#include <math.h>

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
PROTOCCALLSFSUB6(SLA_EVP, sla_evp, DOUBLE, DOUBLE, DOUBLEV, DOUBLEV, DOUBLEV, \
		 DOUBLEV)
#define slaEvp(DATE, DEQX, DVB, DPB, DVH, DPH) \
     CCALLSFSUB6(SLA_EVP, sla_evp, DOUBLE, DOUBLE, DOUBLEV, DOUBLEV, DOUBLEV, \
		 DOUBLEV, DATE, DEQX, DVB, DPB, DVH, DPH)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"
 
double solar_ang(double mjdate, double ra, double dec)
{
    /* Define variables. */
    int i;
    /* double dvb[3], dpb[3], dvh[3], dph[3], vect[3], ra_sun, dec_sun; */
    double dvb[3], dpb[3], dvh[3], dph[3], ra_sun, dec_sun;
    /*
     *  Evp returns four 3-vectors containing the barycentric velocity and
     *  position (dvb,dpv) and the heliocentric velocity and position 
     *  (dvh,dph) of the Earth on the date mjdate.  It requires a modified 
     *  Julian date as input.  Output has units of AU for positions and 
     *  AU/s for velocities.
     */
    slaEvp(mjdate, 2000.0, dvb, dpb, dvh, dph);

    /*  Convert the 3-vector position of the Sun into a J2000.0 RA and DEC */

    for (i=0; i<3; i++) dph[i]*=-1.0;

#ifdef CFORTRAN
    slaDcc2s(dph, ra_sun, dec_sun);
#else
    slaDcc2s(dph, &ra_sun, &dec_sun);
#endif

    /*  Compute the angular separation of RA,DEC and RA_SUN and DEC_SUN */

    return slaDsep(ra*RADIAN, dec*RADIAN, slaDranrm(ra_sun), dec_sun)/RADIAN;
}
