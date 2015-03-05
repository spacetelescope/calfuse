/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    double state_limb(double pos[3], double mjdate,
 *              double ra, double dec, double *zdis,tint day_limb)
 *
 * Description: Computes the earth limb angle from J2000 RA and DEC given
 *              the 3-vector position of FUSE.
 *
 * Arguments:   double pos            (km) 3-vector position of satellite
 *              double mjdate         (days) Modified Julian date 
 *              double ra,dec         (deg) J2000.0 RA and DEC
 *              double zdist          (deg) zenith angle
 *              integer day_limb      output : 1 if bright limb, 0 if dark limb
 *
 * Returns:     double                (deg) Earth limb angle
 *              integer               day_limb will have the value 1 or 0 when
 *                                    returning in the main program.
 *                    
 * Calls:       slaPreces             preces.f
 *		slaDcc2s
 *		slaDranrm
 *		slaEvp
 *
 * History:     03/10/98  E. Murphy   Begin work.
 *              03/10/98  E. Murphy   Initial version working
 *              07/07/99  E. Murphy   Changed call to use pos instead of x,y,z
 *              08/26/99  E. Murphy   Added zenith distance.
 *              05/31/00        peb   Implemented cfortran.h calls for slalib
 *                                    functions.
 *              Ake, T. 1998 in The Scientific Impact of the Goddard
 *                   High Resolution Spectrograph, ed. J. C. Brandt et al.,
 *                   ASP Conference Series, in preparation.
 *
 *              09/23/00   v1.5  ma   Now this function determine bright or
 *                                    dark limb
 *                                    To do so, the procedure has a new
 *                                    input argument: day_limb.
 *                                    day_limb = 1 if bright, 0 if dark.
 *              09/24/00   v1.6  ma   Fixed a bug in the input argument
 *              09/27/00   v1.7  ma   Fixed bug in limbvec calc.
 *		10/02/00   v1.8  jwk  Add fortran wrappers for slaDcc2s,
 *                                    slaEvp, and slaDranrm
 *              12/18/03  bjg  1.3    Change calfusettag.h to calfuse.h
 *              06/17/04  bjg         Corrected cfortran call to sla_preces
 *              07/22/04  bjg  1.4    Remove unused variables  
 ****************************************************************************/
 
#include <math.h>

#ifdef CFORTRAN
#include "cfortran.h"
PROTOCCALLSFSUB5(SLA_PRECES, sla_preces, STRING, DOUBLE, DOUBLE, PDOUBLE, \
		 PDOUBLE)
#define slaPreces(SYSTEM, EP0, EP1, RA, DC) \
     CCALLSFSUB5(SLA_PRECES, sla_preces, STRING, DOUBLE, DOUBLE, PDOUBLE, \
		 PDOUBLE, SYSTEM, EP0, EP1, RA, DC)
PROTOCCALLSFSUB3(SLA_DCC2S, sla_dcc2s, DOUBLEV, PDOUBLE, PDOUBLE)
#define slaDcc2s(V, A, B) \
     CCALLSFSUB3(SLA_DCC2S, sla_dcc2s, DOUBLEV, PDOUBLE, PDOUBLE, V, A, B)
PROTOCCALLSFFUN1(DOUBLE, SLA_DRANRM, sla_dranrm, DOUBLE)
#define slaDranrm(ANGLE) \
     CCALLSFFUN1(SLA_DRANRM, sla_dranrm, DOUBLE, ANGLE)
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

double state_limb(double pos[3], double mjdate, 
                  double ra, double dec, double *zdist,int *day_limb)
{
    double dvb[3], dpb[3], dvh[3], dph[3], limbvec[3], epoch2;
    double ra_earth, dec_earth, r, ndist; 
    double ra_sun, dec_sun;
    int    i;    
    char fk5_st[10]="FK5";

    ra*=RADIAN;
    dec*=RADIAN;

    r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    ra_earth=atan2(-pos[1],-pos[0]);
    dec_earth=asin(-pos[2]/r);
    /*  
     *  We must precess the J2000 RA and DEC to date of observation. This
     *  is a very small correction, and could probably be ignored.
     */
    epoch2=2000.0-((51544.0-mjdate)/365.25);

#ifdef CFORTRAN
    slaPreces(fk5_st, 2000.0, epoch2, ra_earth, dec_earth);
#else
    slaPreces(fk5_st, 2000.0, epoch2, &ra_earth, &dec_earth);
#endif

    /* Now we need the position of the sun (same code as in eclipse.c) */

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

    /* ndist is the angular distance (RADIAN) from the nadir to the target */
    ndist = (acos(sin(dec)*sin(dec_earth)+
                  cos(dec)*cos(dec_earth)*cos(ra-ra_earth)));

    /* Now we want to calculate the limb vector (originating at the center
       of Earth and crossing the target-FUSE line perpendicularily)*/

    limbvec[2]=pos[2]+r*cos(ndist)*sin(dec);
    limbvec[1]=pos[1]+r*cos(ndist)*cos(dec)*sin(ra);
    limbvec[0]=pos[0]+r*cos(ndist)*cos(dec)*cos(ra);

    /* Now we do the scalar product with the sun vector; then
       if the result is <0 the limb is in the day zone (less
       than 90 degrees between sun vector and limb vector)*/

    if ((limbvec[0]*cos(dec_sun)*cos(ra_sun) + 
	 limbvec[1]*cos(dec_sun)*sin(ra_sun)+limbvec[2]*sin(dec_sun)) >0) {
	*day_limb=1;
    } else {
	*day_limb=0;
    }
    *zdist = 180.0 - ndist/RADIAN;

    return  ndist/RADIAN-(asin(RE/r)/RADIAN);
}
