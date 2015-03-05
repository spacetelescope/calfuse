/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    double space_vel(double *vel, double ra, double dec)
 *
 * Description: Computes the projection of the spacecraft orbital velocity 
 *		vector onto a unit vector in the direction ra, dec.
 *
 * Arguments:   double *vel           (km/s) 3-vector velocity of satellite
 *              double ra,dec         (deg) position of target
 *
 * Returns:     double                (km/s) velocity in direction of target
 *                    
 * History:     03/04/98  E. Murphy   Begin work.
 *              03/04/98  E. Murphy   Initial version working
 *              07/07/99  E. Murphy   Converted x,y,z and vx,vy,vz to pos,vel
 *              06/15/01  V. Dixon    Comment out calculation of r, not used.
 *              12/18/03  bjg         Change calfusettag.h to calfuse.h
 *		07/21/04  1.4  wvd    Delete defn of r, as it is not used.
 *		03/04/07  1.5  wvd    Rewrite program to compute dot product
 *				      of target and velocity vectors.
 *
 ****************************************************************************/
 
#include <math.h>
#include "calfuse.h"
 
double space_vel(double vel[3], double ra, double dec)
{

    double x, y, z, phi;

    ra  *= RADIAN;
    dec *= RADIAN;

    z = sin(dec);
    phi = cos(dec);
    y = phi * sin(ra);
    x = phi * cos(ra);

    return x*vel[0] + y*vel[1] + z*vel[2];

}
