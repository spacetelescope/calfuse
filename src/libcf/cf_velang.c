/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    cf_velang(fitsfile *outfits, double mjd)
 *
 * Description: Calculates the sun angle, moon angle, and conversions to
 *              heliocentric and LSR velocities at given mjd, which should be
 *              the middle of the exposure.  For FUV images, the needed
 *              data is in the primary header, so we have to go back
 *              up to the top.  For FES images, the data should be in
 *              each HDU, so you should call cf_velang_calc instead.
 *
 * Arguments:   fitsfile  *outfits       Pointer to FITS file containing the
 *                                      input data and orbital elements
 *              double    mjd           The Modified Julian Date of the 
 *                                      photon arrival time.
 *
 * Calls:       slaDsep         dsep.f
 *              slaRdplan       rdplan.f
 *
 * History:     07/13/98        emm     Begin work
 *              06/10/99        emm     Added keyword header update, added
 *                                      calculation of MOONANG.
 *              06/11/99        emm     Fixed some bugs that Jayant found.
 *              06/21/99         jm     Corrected name
 *              07/16/99        emm     Changed GEO_LON to GEO_LONG
 *              07/21/99        peb     RA and Dec now use RA_TARG and DEC_TARG
 *              07/22/99	emm     Added V_GEOCEN update.
 *              07/23/99	emm     Added cf_velang_calc for FES images.
 *              05/31/00        peb	Implemented cfortran.h calls for slalib
 *                                 	functions.
 *		04/01/03  v1.3	wvd	Replaced printf with cf_verbose
 *		04/04/03  v1.4	wvd	Modifed calls to cf_verbose.
 *              12/18/03  v1.5  bjg     Change calfusettag.h to calfuse.h
 *              07/06/06  v1.6  wvd     Add call to pole_ang and populate 
 *					header keyword POLEANGL.
 *              03/07/07  v1.7  wvd     Remove pos from call to space_vel. 
 *
 ****************************************************************************/

#include <stdio.h>

#ifdef CFORTRAN
#include "cfortran.h"
PROTOCCALLSFFUN4(DOUBLE, SLA_DSEP, sla_dsep, DOUBLE, DOUBLE, DOUBLE, \
		 DOUBLE)
#define slaDsep(A1, B1, A2, B2) \
     CCALLSFFUN4(SLA_DSEP, sla_dsep, DOUBLE, DOUBLE, DOUBLE, DOUBLE, \
		 A1, B1, A2, B2)
PROTOCCALLSFSUB7(SLA_RDPLAN, sla_rdplan, DOUBLE, INT, DOUBLE, DOUBLE, \
		 PDOUBLE, PDOUBLE, PDOUBLE)
#define slaRdplan(DATE, NP, ELONG, PHI, RA, DEC, DIAM) \
     CCALLSFSUB7(SLA_RDPLAN, sla_rdplan, DOUBLE, INT, DOUBLE, DOUBLE, \
		 PDOUBLE, PDOUBLE, PDOUBLE, \
		 DATE, NP, ELONG, PHI, RA, DEC, DIAM)
#else
#include "slalib.h"
#include "slamac.h"
#endif

#include "calfuse.h"
#include "sgp4.h"

SGP4 set_orbit_parms(fitsfile *);

static void
cf_velang_calc(fitsfile *outfits, double mjd)
{
    int    status=0;
    char   comment[FLEN_CARD];
    double ra, dec, geo_lon, geo_lat, mag_lat, pos[3], vel[3], gmst, helvel;
    double ra_moon, dec_moon, moon_ang, dia_moon, vlsr1, vlsr2, sun_ang;
    double geo_vel, pole_angle;
    SGP4   sgp4;

    /* get the state vector at time mjd */
    sgp4 = set_orbit_parms(outfits);
    SGP4_getStateVector(sgp4, mjd, pos, vel);
    SGP4_precess(pos, mjd, MJD2000); SGP4_precess(vel, mjd, MJD2000);

    cf_verbose(2,"State vector at MJD=%14.6f days", mjd);
    cf_verbose(2,"x=%14.6f y=%14.6f z=%14.6f km", pos[0], pos[1], pos[2]);
    cf_verbose(2,"vx=%14.6f vy=%14.6f vz=%14.6f km/s",vel[0], vel[1], vel[2]);

    state_geod(pos, mjd, &geo_lon, &geo_lat, &gmst);

    mag_lat = geod_mag(geo_lon, geo_lat);
    /*
     *  Note: this function has been changed to use RA_TARG and DEC_TARG
     *  keywords, because the RA_APER and DEC_APER keywords are currently
     *  not being set.
     */
    FITS_read_key(outfits, TDOUBLE, "RA_TARG", &ra, comment, &status);
    FITS_read_key(outfits, TDOUBLE, "DEC_TARG", &dec, comment, &status);

    geo_vel= space_vel(vel, ra, dec);
    helvel = helio_vel(mjd, ra, dec);
    vlsr1  = lsrk_vel(ra, dec);
    vlsr2  = lsrd_vel(ra, dec);

    /* Calculate sun angle */
    sun_ang = solar_ang(mjd, ra, dec);

    /* Calculate moon angle */
#ifdef CFORTRAN
    slaRdplan(mjd, 3, geo_lon, geo_lat, ra_moon, dec_moon, dia_moon);
#else
    slaRdplan(mjd, 3, geo_lon, geo_lat, &ra_moon, &dec_moon, &dia_moon);
#endif
    moon_ang = slaDsep(ra*RADIAN, dec*RADIAN, ra_moon, dec_moon)/RADIAN;

    /* Calculate sun angle */
    pole_angle = pole_ang(pos, vel, ra, dec);

    FITS_update_key(outfits, TDOUBLE, "SUNANGLE", &sun_ang, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "MOONANGL", &moon_ang, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "POLEANGL", &pole_angle, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "GEO_LONG", &geo_lon, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "GEO_LAT",  &geo_lat, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "MAG_LAT",  &mag_lat, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "V_GEOCEN", &geo_vel, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "V_HELIO",  &helvel, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "V_LSRDYN", &vlsr2, NULL, &status);
    FITS_update_key(outfits, TDOUBLE, "V_LSRSTD", &vlsr1, NULL, &status);
}


void cf_velang(fitsfile *outfits, double mjd)
{
    int    status=0, hdunum=0, hdutype;
    /*
     *  The outfits file pointer may be down in one of the extensions
     *  whereas the orbital info is in the primary hdu.  Therefore,
     *  record the current position, move to the primary header, 
     *  read the orbital data, and return to the previous hdu.
     */
    FITS_get_hdu_num(outfits, &hdunum);
    FITS_movabs_hdu(outfits, 1, &hdutype, &status);

    cf_velang_calc(outfits, mjd);

    FITS_movabs_hdu(outfits, hdunum, &hdutype, &status);

}
