/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * Synopsis:	SGP4 set_orbit_parms(fitsfile *infits)
 *
 * Description: Initialized the sgp4 structure and reads the orbital 
 *              paramaters from the infits header and places the data 
 *              in the sgp4 structure.  It also calculates a number of
 *              time invariant quantities in the sgp4 structure.
 *
 * Arguments:	fitsfile   *infits        Input FITS file.
 *
 * Return:      SGP4       sgp4	          Structure with orb parms.
 *
 * History:	07/07/99	emurphy	 Begin and finished work
 *              07/14/99        emurphy  Fixed minor bugs
 *              07/21/99        peb      Changed function to return SGP4
 *                                       pointer. (This fixes a bug.)
 *              08/05/99        emurphy  Converted to use set_orbit_parms_calc
 *                                       so that FES pipeline can use these 
 *                                       routines.
 *              12/18/03        bjg      Change calfusettag.h to calfuse.h
 *              04/07/07  1.4	wvd	 Delete CF_PRGM_ID and CF_VER_NUM, 
 *					 as they are not used.
 *
 ****************************************************************************/

#include <stdio.h>
#include "calfuse.h"
#include "sgp4.h"

SGP4 set_orbit_parms_calc(fitsfile *infits)
{
    int    status=0;
    char   comment[FLEN_CARD];
    double epchtime, td, tf, inclinat, eccentry, meananom, argperig, rascascn;
    double n0dt, n0dt2, bstar, mean_motion;
    SGP4   sgp4;
             
    /*  Get the orbital data from the header of the infits file. */
    FITS_read_key(infits, TDOUBLE, "EPCHTIMD", &td, comment, &status);
    FITS_read_key(infits, TDOUBLE, "EPCHTIMF", &tf, comment, &status);
    epchtime=td+tf;    
    FITS_read_key(infits, TDOUBLE, "INCLINAT", &inclinat, comment, &status);
    FITS_read_key(infits, TDOUBLE, "ECCENTRY", &eccentry, comment, &status);
    FITS_read_key(infits, TDOUBLE, "MEANANOM", &meananom, comment, &status);
    FITS_read_key(infits, TDOUBLE, "ARGPERIG", &argperig, comment, &status);
    FITS_read_key(infits, TDOUBLE, "RASCASCN", &rascascn, comment, &status);
    FITS_read_key(infits, TDOUBLE, "FDM2COEF", &n0dt, comment, &status);
    FITS_read_key(infits, TDOUBLE, "SDM6COEF", &n0dt2, comment, &status);
    FITS_read_key(infits, TDOUBLE, "DRAGCOEF", &bstar, comment, &status);
    FITS_read_key(infits, TDOUBLE, "MEANMOTN", &mean_motion, comment, &status);

    sgp4 = SGP4_create();

    SGP4_set(sgp4, epchtime, n0dt, n0dt2,  bstar, inclinat,
	     rascascn, eccentry, argperig, meananom, mean_motion);

    /* do time invariant initializations */
    SGP4_init(sgp4);
    return sgp4;
}

SGP4 set_orbit_parms(fitsfile *infits)
{
    int    status=0, hdunum=0, hdutype;
    SGP4   sgp4;
    /*
     *  The infits file pointer may be down in one of the extensions
     *  whereas the orbital info is in the primary hdu.  Therefore,
     *  record the current position, move to the primary header, 
     *  set the orbital data, and return to the previous hdu.
     */
    FITS_get_hdu_num(infits, &hdunum);
    FITS_movabs_hdu(infits, 1, &hdutype, &status);

    sgp4=set_orbit_parms_calc(infits);

    FITS_movabs_hdu(infits, hdunum, &hdutype, &status);

    return sgp4;

}
