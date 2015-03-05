
/*****************************************************************************
 *	      Johns Hopkins University
 *	      Center For Astrophysical Sciences
 *	      FUSE
 *****************************************************************************
 *
 * 
 * From Modified Julian Day and RA,DEC of target, 
 * computes Heliocentic Modified Julian Day 
 *
 * Uses slalib to compute the light time delay along the target direction.
 *
 *
 * History:	10/06/04	bjg           
 *
 ****************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#ifdef CFORTRAN
#include "cfortran.h"

PROTOCCALLSFSUB6(SLA_DJCL, sla_djcl, DOUBLE, PINT, PINT, PINT, PDOUBLE, \
                 PINT)
#define slaDjcl(DJM, IY, IM, ID, FD, J) \
     CCALLSFSUB6(SLA_DJCL, sla_djcl, DOUBLE, PINT, PINT, PINT, PDOUBLE, \
                 PINT, DJM, IY, IM, ID, FD, J)


PROTOCCALLSFSUB6(SLA_CALYD, sla_calyd, INT, INT, INT, PINT, PINT, \
                 PINT)
#define slaCalyd(IY, IM, ID, NY, ND, J) \
     CCALLSFSUB6(SLA_CALYD, sla_calyd, INT, INT, INT, PINT, PINT, \
                 PINT, IY, IM, ID, NY, ND, J)

PROTOCCALLSFSUB7(SLA_ECOR, sla_ecor, FLOAT, FLOAT, INT, INT, FLOAT, PFLOAT, \
                 PFLOAT)
#define slaEcor(RM, DM, IY, ID, FD, RV, TL) \
     CCALLSFSUB7(SLA_ECOR, sla_ecor,  FLOAT, FLOAT, INT, INT, FLOAT, PFLOAT, \
                 PFLOAT, RM, DM, IY, ID, FD, RV, TL)

#else
#include "slalib.h"
#include "slamac.h"
#endif


#include "calfuse.h"


#define DEG2RAD (M_PI/180.0f)

double gethmjd(double mjd, int ra_h, int ra_m, float ra_s, int dec_d, int dec_m, float dec_s){


  int status;

  int year, month, day;
  int jyear, jday; 
  double frac;
  
  
  float ra, dec;
  float lt, v;
    
  ra  = ( ra_h  +  ra_m/60.0 +  ra_s/3600.0 ) * 15.0 * DEG2RAD;
  dec = ( dec_d + dec_m/60.0 + dec_s/3600.0 ) * DEG2RAD;

  
#ifdef CFORTRAN
  slaDjcl(mjd, year, month, day, frac, status);
#else
  slaDjcl(mjd, &year, &month, &day, &frac, &status);
#endif  
  
#ifdef CFORTRAN
  slaCalyd(year, month, day, jyear, jday, status);
#else
  slaCalyd(year, month, day, &jyear, &jday, &status);
#endif

#ifdef CFORTRAN
  slaEcor(ra, dec, jyear, jday, (float) frac, v, lt);
#else
  slaEcor(ra, dec, jyear, jday, (float) frac, &v, &lt);
#endif

  return (double) (mjd+lt/(3600.0*24.0));

}
