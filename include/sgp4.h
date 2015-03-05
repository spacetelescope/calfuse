/* H+
 *      Title     : sgp4.h
 *      Author    : Bryce A. Roberts
 *      Date      : 23 February 1999
 *      Synopsis  : header file for SGP4 propagation routines
 *      SCCS      : @(#)sgp4.h	1.1 03/04/03
 *      Revisions :
 *      mm/dd/yy name   description
 * H-
 */

#ifndef _SGP4_H
#define _SGP4_H

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>

#define MJD2000      (51544.0)

/* the data structure itself */
struct sgp4_st {
  /* NORAD tle components */
  double epochTime, n0, n0dt, n0dt2, bstar, i0, raan, e0, w0, M0;

  /* derived variables */
  double n0dp, a0dp, q0ms4, theta, theta_2, zeta_4, beta0_2, eta, C1, 
    C1_2, C1_3, C1_4, C3, C4, C5, D2, D3, D4, mdot, wdot, raandot; 
  };

typedef struct sgp4_st * SGP4;

/* create a new instance of SGP4 */
extern SGP4 SGP4_create(void);

/* destroy an instance of SGP4 */
extern void SGP4_destroy(SGP4 sgp4);

/* convert NORAD time string to TJD */
extern int SGP4_getTime(const char *str, double *t);

/* do time-invariant initializations */
extern void SGP4_init(SGP4 sgp4);

/* get a state vector at a point in time */
extern void SGP4_getStateVector(SGP4 sgp4, double t, double pos[3], 
			       double vel[3]);

/* finds the pole (orbit plane normal) */
extern void SGP4_getPole(double pos[3], double vel[3], double pole[3]);

/* sets the elements of an instance of SGP4 */
extern void SGP4_set(SGP4 sgp4, double epochTime, double n0dt, double n0dt2,
		     double bstar, double i0, double raan, double e0,
		     double w0, double M0, double n0);

/* precesses a direction vector from t1 to t2 */
void SGP4_precess(double v[3], double t1, double t2);


#endif
