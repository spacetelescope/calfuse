/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    double geod_mag(double lon, double lat)
 *
 * Description: Computes the magnetic latitude of FUSE from the
 *              given geocentric longitude and latitude.
 *
 * Arguments:   double lon,lat        (deg) Geocentric longitude and latitude
 *
 * Returns:     double                Geomagnetic latitude.
 *                    
 * History:     03/11/98  E. Murphy   Begin work.
 *              03/11/98  E. Murphy   Initial version working
 *              04/13/99  E. Murphy   Moved PI and RADIAN to calfuse.h
 *              12/18/03  bjg         Change calfusettag.h to calfuse.h
 *
 *              Ake, T. 1998 in The Scientific Impact of the Goddard
 *                   High Resolution Spectrograph, ed. J. C. Brandt et al.,
 *                   ASP Conference Series, in preparation.
 ****************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "calfuse.h"
 
double geod_mag(double lon, double lat)
{
    double lat_rad, c1;

    lat_rad=lat*RADIAN;

    c1=sin(lat_rad)*cos(11.4*RADIAN)-
        cos(lat_rad)*cos((lon+69.8)*RADIAN)*sin(11.4*RADIAN);

    return asin(c1)/RADIAN;
}
