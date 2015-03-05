/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    int eclipse(double pos[3], double mjdate)
 *
 * Description: Determines if FUSE is in the day or night portion of the orbit.
 *              
 *
 * Arguments:   double pos            Cartesian vector xyz coord of satellite.
 *              double mjdate          Julian date of observation
 *
 * Returns:     int                   0 if day portion of orbit
 *                                    1 if night portion of orbit
 *                    
 * Calls:       slaDcc2s              dcc2s.f
 *              slaDranrm             dranrm.f
 *              slaDsep               dsep.f
 *              slaEvp                evp.f
 *              slaPreces             preces.f
 *
 * History:     03/10/98  E. Murphy   Begin work.
 *              03/10/98  E. Murphy   Initial version working
 *              07/07/99  E. Murphy   Replaced x,y,z in call with pos
 *              07/20/99  jm          Removed duplicate define statements
 *              08/26/99  E. Murphy   ang_sep now returned.
 *              05/31/00        peb   Implemented cfortran.h calls for slalib
 *                                    functions.
 *              12/18/03  bjg       Change calfusettag.h to calfuse.h
 *              06/17/04  bjg  1.4  Corrected cfortran call to sla_preces   
 *
 * References:  This routine makes use of the Starlink set of astronomical
 *              subroutines (SLALIB).  More information can be found at
 *              http://star-www.rl.ac.uk.
 ****************************************************************************/
 
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
PROTOCCALLSFFUN4(DOUBLE, SLA_DSEP, sla_dsep, DOUBLE, DOUBLE, DOUBLE, \
		 DOUBLE)
#define slaDsep(A1, B1, A2, B2) \
     CCALLSFFUN4(SLA_DSEP, sla_dsep, DOUBLE, DOUBLE, DOUBLE, DOUBLE, \
		 A1, B1, A2, B2)
PROTOCCALLSFSUB6(SLA_EVP, sla_evp, DOUBLE, DOUBLE, DOUBLEV, DOUBLEV, \
		 DOUBLEV, DOUBLEV)
#define slaEvp(DATE, DEQX, DVB, DPB, DVH, DPH) \
     CCALLSFSUB6(SLA_EVP, sla_evp, DOUBLE, DOUBLE, DOUBLEV, DOUBLEV, \
		 DOUBLEV, DOUBLEV, DATE, DEQX, DVB, DPB, DVH, DPH)
PROTOCCALLSFSUB5(SLA_PRECES, sla_preces, STRING, DOUBLE, DOUBLE, PDOUBLE, \
		 PDOUBLE)
#define slaPreces(SYSTEM, EP0, EP1, RA, DC) \
     CCALLSFSUB5(SLA_PRECES, sla_preces, STRING, DOUBLE, DOUBLE, PDOUBLE, \
		 PDOUBLE, SYSTEM, EP0, EP1, RA, DC)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"

int eclipse(double *pos, double mjdate, double *ang_sep)
{
    /* Define variables. */
    int    i, retcode;
    double dvb[3], dpb[3], dvh[3], dph[3], epoch2;
    double r, ra_earth, dec_earth, ra_sun, dec_sun, ang_earth, ang_sun;
    char fk5_st[10]="FK5";

    /* Calculate the J2000.0 apparent RA and DEC of the Earth */
    r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    ra_earth=atan2(-pos[1],-pos[0]);
    dec_earth=asin(-pos[2]/r);

    /*  We must precess the J2000 RA and DEC to date of observation. This
     *  is a very small correction, and could probably be ignored.
     */

    epoch2=2000.0-((51544.0-mjdate)/365.25);

#ifdef CFORTRAN
    slaPreces(fk5_st, 2000.0, epoch2, ra_earth, dec_earth);
#else
    slaPreces(fk5_st, 2000.0, epoch2, &ra_earth, &dec_earth);
#endif

    /*  Evp returns four 3-vectors containing the barycentric velocity and
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

    ra_sun=slaDranrm(ra_sun);

    /*  We must precess the J2000 RA and DEC to date of observation. This
     *  is a very small correction, and could probably be ignored.
     */

    epoch2=2000.0-((51544.0-mjdate)/365.25);

#ifdef CFORTRAN
    slaPreces(fk5_st, 2000.0, epoch2, ra_sun, dec_sun);
#else
    slaPreces(fk5_st, 2000.0, epoch2, &ra_sun, &dec_sun);
#endif

    /* Compute the angular sizes of the Earth and the Sun. */

    ang_earth=asin(RE/r);

    ang_sun=asin(RS/(sqrt(dph[0]*dph[0]+dph[1]*dph[1]+dph[2]*dph[2])*AU));

    /* Compute the angular separation of the Earth and the Sun. */

    *ang_sep=slaDsep(ra_earth, dec_earth, ra_sun, dec_sun);

    retcode=0;
    if (*ang_sep < ang_earth+ang_sun) retcode=1;

    /* Convert ang_sep into degrees. */
    *ang_sep /= RADIAN;

    return retcode;
}
