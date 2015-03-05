/******************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 ******************************************************************************
 *
 * Synopsis:    void state_geod(double pos,
 *              double mjdate, double lon, double lat)
 *
 * Description: Computes the geocentric longitude and latitude of FUSE
 *              from the given state 6-vector.
 *
 * Arguments:   double pos            (km) 3-vector xyz position of satellite
 *              double mjdate         Modified Julian date of observation
 *
 * Returns:     double lon, lat       (deg) Geocentric longitude and latitude
 *
 * Calls:       slaPreces             preces.f
 *
 * History:     03/04/98  E. Murphy   Begin work.
 *              03/04/98  E. Murphy   Initial version working
 *              03/15/99  E. murphy   Changed function definition to void 
 *                                    (was double)
 *              07/07/99  E. Murphy   Changed call to use pos instead of x,y,z
 *              05/31/00        peb   Implemented cfortran.h calls for slalib
 *                                    functions.
 *              12/18/03    bjg       Change calfusettag.h to calfuse.h
 *              06/17/04  bjg  1.4  Corrected cfortran call to sla_preces  
 *
 *              Ake, T. 1998 in The Scientific Impact of the Goddard
 *                   High Resolution Spectrograph, ed. J. C. Brandt et al.,
 *                   ASP Conference Series, in preparation.
 *****************************************************************************/
 
#include <stdio.h>
#include <math.h>

#ifdef CFORTRAN
#include "cfortran.h"
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

void state_geod(double pos[3], double mjdate, 
		double *lon, double *lat, double *gmst)
{
    double ra, dec, epoch2, c1, r, intmjd, frcmjd;
    char fk5_st[10]="FK5";

    /* pos[0]=x
       pos[1]=y
       pos[2]=z
       */

    r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    ra=atan2(pos[1],pos[0]);
    dec=asin(pos[2]/r);

    /*  We must precess the J2000 RA and DEC to date of observation. This
     *  is a very small correction, and could probably be ignored.
     */

    epoch2=2000.0-((51544.0-mjdate)/365.25);

#ifdef CFORTRAN
    slaPreces(fk5_st, 2000.0, epoch2, ra, dec);
#else
    slaPreces(fk5_st, 2000.0, epoch2, &ra, &dec);
#endif

    /*  Calculate the Greenwich Mean Sidereal Time in seconds 
     *  Astronomical Almanac, 1998, p. B6. 
     */

    frcmjd=modf(mjdate, &intmjd);

    c1=(intmjd-51544.5)/36525.0;

    /* The following line gives the GMST at 0 hr UT */
    *gmst=24110.54841 + 8640184.812866*c1+0.093104*c1*c1-6.2E-6*c1*c1*c1;
    /* The following adds on the number of seconds since the UT above */
    *gmst+=frcmjd*1.00273791*86400.0;

    *gmst*=(2*M_PI)/86400.0;

    while (*gmst < 0.0) *gmst+=2*M_PI;

    while (*gmst > 2*M_PI) *gmst-=2*M_PI;

    *lon=(ra-*gmst)/RADIAN;

    while (*lon > 180.0) *lon-=360.0;

    while (*lon < -180.0) *lon+=360.0;

    *lat=dec/RADIAN;

    *gmst*=24.0/(2*M_PI);
}
