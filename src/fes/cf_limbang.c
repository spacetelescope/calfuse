/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_limbang_calc(fitsfile *outfits, double mjd)
 *              cf_min_limbang(fitsfile *outfits, 
 *                                           double mjd_start, double mjd_end)
 * Description: cf_limbang_calc - Calculates the limb angle at a given mjd
 *              cf_min_limbang  - determines the minimum limb angle between
 *                                 the start and end times. The limb angle is
 *                                 calculated every 8.6 seconds( .0001 MJD)
 *
 * Arguments:   fitsfile  *outfits       Pointer to FITS file containing the
 *                                      input data and orbital elements
 *              double    mjd           The limb angle calculation is done for
 *                                       this time.  The time is given as 
 *                                       a Modified Julian Date.
 *              double mjd_start        The range of times for which the minimum
 *              double mjd_end          limb angle should be determined.
 *                                      Both are given as Modified Julian Dates.
 *
 * History:     08/09/99        mlr     copied cf_velang.c to start
 *              08/09/99        mlr     removed excess calls then made the
 *                                      limb_ang calls - same as cf_check_point
 *
 ****************************************************************************/

#include <stdio.h>
#include "calfuse.h"
#include "sgp4.h"
static char CF_PRGM_ID[] = "cf_limbang";
 
SGP4   set_orbit_parms(fitsfile *);

double cf_limbang_calc(fitsfile *outfits, double mjd)
{
    int    status=0, isday_dummy;
    char   comment[FLEN_CARD];
    double ra, dec, pos[3], vel[3];
    double lim_ang, zdist_dummy;
    SGP4   sgp4;

    sgp4 = set_orbit_parms(outfits);

    /* get the state vector at time mjd */
    SGP4_getStateVector(sgp4, mjd, pos, vel);
    SGP4_precess(pos, mjd, MJD2000); SGP4_precess(vel, mjd, MJD2000);

    FITS_read_key(outfits, TDOUBLE, "RA_TARG", &ra, comment, &status);
    FITS_read_key(outfits, TDOUBLE, "DEC_TARG", &dec, comment, &status);

    lim_ang = state_limb(pos, mjd, ra, dec, &zdist_dummy, &isday_dummy);

    return(lim_ang);
}

void cf_min_limbang(fitsfile *outfits, double mjd_start, double mjd_end)
{
     
     int status = 0;
     double lim_ang, min_lim;
     double mjd;     

    mjd = mjd_start;
    min_lim = cf_limbang_calc(outfits, mjd);

#ifdef DEBUG
    printf("mjd: %lf mjd_end = %lf minlim = %lf\n", mjd, mjd_end, min_lim);
#endif
    /* .0001 = 8.6 seconds */
    while ( (mjd_end - mjd) > .0001 )
    {
          mjd = mjd + .0001;
          lim_ang = cf_limbang_calc(outfits, mjd);
#ifdef DEBUG
          printf("mjd: %lf mjd_end = %lf minlim = %lf, limang = %lf\n", 
                    mjd, mjd_end, min_lim, lim_ang);
#endif
          if (lim_ang < min_lim) min_lim = lim_ang;
    }

    mjd = mjd_end;
    lim_ang = cf_limbang_calc(outfits, mjd);

#ifdef DEBUG
    printf("mjd: %lf mjd_end = %lf minlim = %lf, limang = %lf\n", 
                 mjd, mjd_end, min_lim, lim_ang);
#endif

    if (lim_ang < min_lim)  min_lim = lim_ang;

#ifdef DEBUG
    printf("mjd: %lf mjd_end = %lf minlim = %lf\n", mjd, mjd_end, min_lim);    
#endif

    FITS_update_key(outfits, TDOUBLE, "MIN_LIMB", &min_lim, NULL, &status);

}
