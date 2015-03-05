/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    double helio_vel(double mjdate, double ra, double dec)
 *
 * Description: Computes the Earth's orbital velocity in the direction 
 *              of ra,dec to convert to heliocentric velocity.
 *
 * Arguments:   double mjdate         Mod. Julian date of observation
 *              double ra             (deg) J2000.0 right ascension
 *              double dec            (deg) J2000.0 declination 
 *
 * Returns:     double                (km/s) velocity in direction of source
 *                    
 * Calls:       slaEvp                evp.f
 *              slaDcs2c              dcs2c.f
 *              slaDvdv               dvdv.f
 *
 * History:     03/05/98  E. Murphy   Begin work.
 *              03/05/98  E. Murphy   Initial version working
 *              05/31/00        peb   Implemented cfortran.h calls for slalib
 *                                    functions.
 *              06/15/01  V. Dixon    Correct sign of return value to agree
 *				      with astronomical convention.
 *              12/18/03  bjg         Change calfusettag.h to calfuse.h
 *
 * References:  This routine makes use of the Starlink set of astronomical
 *              subroutines (SLALIB).  More information can be found at
 *              http://star-www.rl.ac.uk.
 ****************************************************************************/
 
#include <stdio.h>

#ifdef CFORTRAN
#include "cfortran.h"

PROTOCCALLSFSUB3(SLA_DCS2C, sla_dcs2c, DOUBLE, DOUBLE, DOUBLEV)
#define slaDcs2c(A, B, V) \
     CCALLSFSUB3(SLA_DCS2C, sla_dcs2c, DOUBLE, DOUBLE, DOUBLEV, A, B, V)

PROTOCCALLSFFUN2(DOUBLE, SLA_DVDV, sla_dvdv, DOUBLEV, DOUBLEV)
#define slaDvdv(VA, VB) \
     CCALLSFFUN2(SLA_DVDV, sla_dvdv, DOUBLEV, DOUBLEV, VA, VB)

PROTOCCALLSFSUB6(SLA_EVP, sla_evp, DOUBLE, DOUBLE, DOUBLEV, DOUBLEV, \
		 DOUBLEV, DOUBLEV)
#define slaEvp(DATE, DEQX, DVB, DPB, DVH, DPH) \
     CCALLSFSUB6(SLA_EVP, sla_evp, DOUBLE, DOUBLE, DOUBLEV, DOUBLEV, \
		 DOUBLEV, DOUBLEV, DATE, DEQX, DVB, DPB, DVH, DPH)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"
 
double helio_vel(double mjdate, double ra, double dec)
{
    /* Define variables. */
    double dvb[3], dpb[3], dvh[3], dph[3], vect[3];
    /*
     *  Evp returns four 3-vectors containing the barycentric velocity and
     *  position (dvb,dpb) and the heliocentric velocity and position 
     *  (dvh,dph) of the Earth on the date mjdate.  It requires a modified 
     *  Julian date as input.  Output has units of AU for positions and 
     *  AU/s for velocities.
     */
    slaEvp(mjdate, 2000.0, dvb, dpb, dvh, dph);

    /*  Convert the J2000.0 RA and DEC into a 3-vector position */

    slaDcs2c(RADIAN*ra, RADIAN*dec, vect);
    /*
     *  Compute the dot product of the RA DEC position vector and the 
     *  heliocentric velocity vector.  Mutliply by km/AU to get velocity 
     *  in km/s.
     */
    return slaDvdv(vect, dvh)*149.5978707E6;
}
