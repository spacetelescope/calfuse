/* H+
 *      Title     : sgp4.c
 *      Author    : Bryce A. Roberts
 *      Date      : 23 February 1999
 *      Synopsis  : implementation of SGP4 propagation routines
 *      SCCS      : @(#)sgp4.c	1.6 05/14/00
 *      Revisions :
 *      mm/dd/yy name   description
 *      01/14/00 ma     get rid of unused function get_time (May 14 2000)
 *      12/18/03 bjg    Change calfusettag.h to calfuse.h 
* H-
 */

#include <string.h>
#include <time.h>
#include "calfuse.h"
#include "sgp4.h"

/* define some absolute constants used by SGP4 */
#define MIN_PER_DAY  (1440.0)
#define KM_PER_EARTH (6378.135)
#define TWOBYTHREE   (2.0/3.0)
#define Q0           (120.0)
#define S0           (78.0)
#define THREEPIBYTWO (PIBY2*3.0)
#define KE           (0.074366916)
#define AE           (1.0)
#define S            (AE*(1.0+S0/KM_PER_EARTH))
#define Q0MS4        pow((Q0-S0)*AE/KM_PER_EARTH, 4)
#define CK2          (0.5*J2*AE*AE)
#define CK4          (-0.375*J4*AE*AE*AE*AE)
#define A30          (-J3/AE*AE*AE)

#ifndef J2
#define J2           (1.082616E-3)
#endif

#ifndef J3
#define J3           (-0.253881E-5)
#endif

#ifndef J4
#define J4           (-1.65597E-6)
#endif

#ifndef PI
#define PI (3.141592653589793)
#endif

#ifndef TWOPI
#define TWOPI (2.0*PI)
#endif

#define RAD2DEG(x) (((x)*180.0)/PI)
#define DEG2RAD(x) ((PI/180.0)*(x))

#ifndef false
#define false (0)
#endif

#ifndef true
#define true (1)
#endif


/*
 * Name:     SGP4_isLeap
 * Purpose:  determine if a given year is a leap year 
 * Input:    int year (e.g. 1999)
 * Output:   returns 1 if leap year, 0 if not
 */
static int SGP4_isLeap(int year) {
  return(((year%4==0 && year%100!=0) || year%400==0));
  }


/*
 * Name:     SGP4_doy2cal
 * Purpose:  convert year/day-of-year to Gregorian calendar month and day
 * Input:    int year, doy
 * Output:   int *month, int *day, returns 0 for fail, 1 for success
 */
static int SGP4_doy2cal(int year, int doy, int *month, int *day) {
  /* number of days in leap and non-leap year months */
  static int daytab[2][13]={{0,31,28,31,30,31,30,31,31,30,31,30,31},
                            {0,31,29,31,30,31,30,31,31,30,31,30,31}};
  int leap;

  /* do range checking */
  if (doy<1 || doy>366)
    return(false);
  
  leap=SGP4_isLeap(year) ? 1 : 0;
  *day=doy;
  for (*month=1; *day>daytab[leap][*month]; (*month)++)
    *day-=daytab[leap][*month];
  
  return(true);
  }


/*
 * Name:     SGP4_create
 * Purpose:  creates an instance of SGP4
 * Input:    none
 * Output:   SGP4 
 */
SGP4 SGP4_create(void) {
  return((SGP4)calloc(1, sizeof(struct sgp4_st)));
  }


/*
 * Name:     SGP4_destroy
 * Purpose:  destroys an instance of SGP4
 * Input:    SGP4 sgp4
 * Output:   none
 */
void SGP4_destroy(SGP4 sgp4) {
  if (sgp4)
    free(sgp4);
  }


/*
 * Name:     SGP4_init
 * Purpose:  initializes time-invariant SGP4 variables
 * Input:    SGP4 sgp4
 * Output:   none
 */
void SGP4_init(SGP4 sgp4) {
  double a1, del1, a0, del0, perigee, sStar, eta_2, eta_3, eta_4, pinvsq, 
    temp1, temp2, temp3, x3thetam1, x1m5th, xhdot1, theta_4, zeta, zeta_2, 
    zeta_3, zeta_5, beta0, beta0_3, beta0_4, C2;

  /* recover original mean motion and semimajor axis from the elements */
  a1=pow(KE/sgp4->n0, TWOBYTHREE);
  del1=1.5 * (CK2/(a1*a1))*((3*cos(sgp4->i0)*cos(sgp4->i0)-1.0)/
    pow(1-sgp4->e0*sgp4->e0, 1.5));
  a0=a1*(1-del1/3.0-del1*del1-(134.0/81.0)*del1*del1*del1);
  del0=(del1*a1*a1)/(a0*a0);

  /* find original mean motion, n0dp (n0'') */
  sgp4->n0dp=sgp4->n0/(1+del0); 

  /* find semimajor axis, a0dp (a0'') */
  sgp4->a0dp=a0/(1-del0);

  /* find the perigee (in km) */
  perigee=(sgp4->a0dp*(1.0-sgp4->e0)-AE)*KM_PER_EARTH;
  
  /* make decisions on what value of s to use based on perigee */
  if (perigee<=156.0 && perigee>=98.0) {
    sStar=sgp4->a0dp*(1.0-sgp4->e0)-S-AE;
    sgp4->q0ms4=pow(pow(Q0MS4, 0.25) + S - sStar, 4);
    }
  else if (perigee<98.0) {
    sStar=20/KM_PER_EARTH + AE;
    sgp4->q0ms4=pow(pow(Q0MS4, 0.25) + S - sStar, 4);
    }
  else {
    sStar=S;
    sgp4->q0ms4=Q0MS4;
    }

  sgp4->theta=cos(sgp4->i0); sgp4->theta_2=sgp4->theta*sgp4->theta; 
  theta_4=sgp4->theta_2*sgp4->theta_2;
  zeta=1/(sgp4->a0dp-sStar);
  
  zeta_2=zeta*zeta, zeta_3=zeta_2*zeta; 
  sgp4->zeta_4=zeta_3*zeta;
  zeta_5=sgp4->zeta_4*zeta;

  sgp4->beta0_2=1-sgp4->e0*sgp4->e0;
  beta0=sqrt(sgp4->beta0_2); beta0_3=sgp4->beta0_2*beta0; 
  beta0_4=beta0_3*beta0;

  sgp4->eta=sgp4->a0dp*sgp4->e0*zeta;
  eta_2=sgp4->eta*sgp4->eta, eta_3=eta_2*sgp4->eta, eta_4=eta_3*sgp4->eta;
  
  C2=sgp4->q0ms4*(sgp4->zeta_4)*sgp4->n0dp*pow(1-eta_2, -7.0/2.0)*
    (sgp4->a0dp*(1+ ((3.0/2.0)*eta_2) + 4*sgp4->e0*sgp4->eta + 
    sgp4->e0*eta_3) + (3.0/2.0)*(CK2*zeta/(1-eta_2))*(-0.5 + 
    1.5*sgp4->theta_2)*(8.0+24*eta_2 + 3*eta_4));
  sgp4->C1=sgp4->bstar*C2; sgp4->C1_2=sgp4->C1*sgp4->C1; 
  sgp4->C1_3=sgp4->C1_2*sgp4->C1; sgp4->C1_4=sgp4->C1_3*sgp4->C1;
  sgp4->C3=(sgp4->q0ms4*(zeta_5)*A30*sgp4->n0dp*AE*sin(sgp4->i0))/
    (CK2*sgp4->e0);
  sgp4->C4=2*sgp4->n0dp*sgp4->q0ms4*(sgp4->zeta_4)*sgp4->a0dp*(sgp4->beta0_2)*
    pow(1-eta_2, -7.0/2.0)*( (2*sgp4->eta*(1+sgp4->e0*sgp4->eta)+
    0.5*sgp4->e0+0.5*eta_3)-2*CK2*zeta/(sgp4->a0dp*(1-eta_2))*
   (3*(1-3*sgp4->theta_2)*(1 + 1.5*eta_2 -2*sgp4->e0*sgp4->eta - 
   0.5*sgp4->e0*eta_3) + 0.75*(1-sgp4->theta_2)*(2*eta_2 - 
   sgp4->e0*sgp4->eta - sgp4->e0*eta_3)*cos(2*sgp4->w0)));
  sgp4->C5=2*sgp4->q0ms4*(sgp4->zeta_4)*sgp4->a0dp*(sgp4->beta0_2)*
    pow(1-eta_2, -7.0/2.0)*(1 + (11.0/4.0)*sgp4->eta*(sgp4->eta+sgp4->e0)+
    sgp4->e0*eta_3);

  sgp4->D2=4*sgp4->a0dp*zeta*sgp4->C1_2;
  sgp4->D3=(4.0/3.0)*sgp4->a0dp*zeta_2*(17*sgp4->a0dp + 
    sStar)*sgp4->C1_3;
  sgp4->D4=TWOBYTHREE*sgp4->a0dp*zeta_3*(221*sgp4->a0dp+
    31*sStar)*sgp4->C1_4;

  /* constants for secular effects calculation */
  pinvsq=1.0/(sgp4->a0dp*sgp4->a0dp*beta0_4);
  temp1=3.0*CK2*pinvsq*sgp4->n0dp;
  temp2=temp1*CK2*pinvsq;
  temp3=1.25*CK4*pinvsq*pinvsq*sgp4->n0dp;
  x3thetam1=3.0*sgp4->theta_2 - 1.0;
  x1m5th=1.0 - 5.0*sgp4->theta_2;
  xhdot1=-temp1*sgp4->theta;

  sgp4->mdot=sgp4->n0dp+0.5*temp1*beta0*x3thetam1+
    0.0625*temp2*beta0*(13.0 - 78.0*sgp4->theta_2+137.0*theta_4);

  sgp4->wdot=-0.5*temp1*x1m5th+0.0625*temp2*(7.0-114.0*sgp4->theta_2+
    395.0*theta_4)+temp3*(3.0-36.0*sgp4->theta_2+49.0*theta_4);
  
  sgp4->raandot=xhdot1+(0.5*temp2*(4.0-19.0*sgp4->theta_2)+2.0*temp3*
    (3.0-7.0*sgp4->theta_2))*sgp4->theta;
  }


/*
 * Name:     SGP4_getStateVector
 * Purpose:  finds position, velocity at a time
 * Input:    SGP4 sgp4, double t
 * Output:   double pos[3], double vel[3]
 */
void SGP4_getStateVector(SGP4 sgp4, double t, double pos[3], double vel[3]) {
  double mdf, wdf, raandf, del_w, del_m, mp, w, RAAN, e, a, IL,
    beta_2, n, axn, ILl, aynl, ILt, ayn, U, initial, eCosE,
    eSinE, eL_2, pL, r, rDot, rfDot, cosu, sinu, u, rk, uk, RAANk, ik,
    rDotk, rfDotk, mx, my, ux, uy, uz, vx, vy, vz, tsince;
  int i;

  /* find the time since the epoch, in minutes */
  tsince=(t-sgp4->epochTime)*1440;

  /* secular effects of atmospheric drag and gravitation */
  mdf=sgp4->M0+(sgp4->mdot*tsince);
  wdf=sgp4->w0+(sgp4->wdot*tsince);
  raandf=sgp4->raan+(sgp4->raandot*tsince);
  
  del_w=sgp4->bstar*sgp4->C3*(cos(sgp4->w0))*tsince;
  del_m=-TWOBYTHREE*sgp4->q0ms4*sgp4->bstar*sgp4->zeta_4*
    (AE/(sgp4->e0*sgp4->eta))*(pow(1+sgp4->eta*cos(mdf), 3)-
    pow(1+sgp4->eta*cos(sgp4->M0), 3));
  
  mp=mdf+del_w+del_m;
  w=wdf-del_w-del_m;
  
  RAAN=raandf-(21.0/2.0)*((sgp4->n0dp*CK2*sgp4->theta)/(sgp4->a0dp*
    sgp4->a0dp*sgp4->beta0_2))* sgp4->C1*tsince*tsince;
  e=sgp4->e0-sgp4->bstar*sgp4->C4*tsince-sgp4->bstar*sgp4->C5*
    (sin(mp)-sin(sgp4->M0));
  
  a=sgp4->a0dp*pow((1+tsince*(-sgp4->C1+tsince*(-sgp4->D2+tsince*
    (-sgp4->D3-sgp4->D4*tsince)))), 2);
  
  IL=mp+w+RAAN+sgp4->n0dp*tsince*tsince*(1.5*sgp4->C1+tsince*
    ((sgp4->D2+2*sgp4->C1_2) + tsince*(0.25*(3*sgp4->D3 + 
    12*sgp4->C1*sgp4->D2 + 10*sgp4->C1_3)+tsince*(0.2*(3*sgp4->D4+12*
    sgp4->C1*sgp4->D3+6*sgp4->D2*sgp4->D2+30*sgp4->C1_2*sgp4->D2+
    15*sgp4->C1_4)))));

  beta_2=1-e*e;
  n=KE/pow(a, 3.0/2.0);

  /* find the long-period periodic terms */
  axn=e*cos(w);
  ILl=(A30*sin(sgp4->i0))/(8*CK2*a*beta_2)*(e*cos(w))*((3+5*sgp4->theta)/
    (1+sgp4->theta));
  aynl=(A30*sin(sgp4->i0))/(4*CK2*a*beta_2);
  ILt=IL+ILl;
  ayn=e*sin(w)+aynl;
    
  /* iteratively solve Kepler's equation */
  U=fmod(ILt-RAAN, TWOPI);

  initial=U;
  for (i=0; i<10; i++) {
    U=(initial-ayn*cos(U)+axn*sin(U)-U)/(1.0-ayn*sin(U)-axn*cos(U))+U;
    }

  /* preliminary quantities for short period periodics */
  eCosE=axn*cos(U)+ayn*sin(U); 
  eSinE=axn*sin(U)-ayn*cos(U);
  eL_2=axn*axn+ayn*ayn;
  pL=a*(1-eL_2);
  r=a*(1-eCosE);
  rDot=KE*(sqrt(a)/r)*eSinE;
  rfDot=KE*sqrt(pL)/r;
  cosu=(a/r)*(cos(U)-axn+(ayn*eSinE)/(1+sqrt(1-eL_2)));
  sinu=(a/r)*(sin(U)-ayn-(axn*eSinE)/(1+sqrt(1-eL_2)));
  u=fmod(atan2(sinu, cosu)+TWOPI, TWOPI);

  /* update for short-period periodics */
  rk=r*(1-1.5*CK2*sqrt(1-eL_2)/(pL*pL)*(3*sgp4->theta_2-1))+
    (CK2/(2*pL))*(1-sgp4->theta_2)*cos(2*u);
  uk=u+-(CK2/(4*pL*pL))*(7*sgp4->theta_2-1.0)*sin(2*u);;
  RAANk=RAAN+((3*CK2*sgp4->theta)/(2*pL*pL))*sin(2*u);
  ik=sgp4->i0+((3*CK2*sgp4->theta)/(2*pL*pL))*sin(sgp4->i0)*cos(2*u);
  rDotk=rDot-((CK2*n)/pL)*(1-sgp4->theta_2)*sin(2*u);
  rfDotk=rfDot+((CK2*n)/pL)*((1-sgp4->theta_2)*cos(2*u)-
  1.5*(1-3*sgp4->theta_2));
  
  /* find the position and velocity components */
  mx=-sin(RAANk)*cos(ik);
  my=cos(RAANk)*cos(ik);
  
  ux=mx*sin(uk)+cos(RAANk)*cos(uk);
  uy=my*sin(uk)+sin(RAANk)*cos(uk);
  uz=sin(ik)*sin(uk);
  
  vx=mx*cos(uk)-cos(RAANk)*sin(uk);
  vy=my*cos(uk)-sin(RAANk)*sin(uk);
  vz=sin(ik)*cos(uk);

  pos[0]=KM_PER_EARTH*rk*ux;
  pos[1]=KM_PER_EARTH*rk*uy;
  pos[2]=KM_PER_EARTH*rk*uz;

  vel[0]=MIN_PER_DAY*KM_PER_EARTH*(rDotk*ux+rfDotk*vx)/(AE*86400.0);
  vel[1]=MIN_PER_DAY*KM_PER_EARTH*(rDotk*uy+rfDotk*vy)/(AE*86400.0);
  vel[2]=MIN_PER_DAY*KM_PER_EARTH*(rDotk*uz+rfDotk*vz)/(AE*86400.0);
  }


/*
 * Name:     SGP4_set
 * Purpose:  sets the elements of an instance of SGP4
 * Input:    SGP4 sgp4, double epochTime, double n0dt, double n0dt2,
 *             double bstar, double i0, double raan, double e0,
 *             double w0, double M0, double n0
 * Output:   none
 */
void SGP4_set(SGP4 sgp4, double epochTime, double n0dt, double n0dt2,
	      double bstar, double i0, double raan, double e0,
	      double w0, double M0, double n0) {

  /* set the epoch time */
  sgp4->epochTime=epochTime;

  /* set the first time derivative of mean motion */
  sgp4->n0dt=n0dt*(TWOPI/(MIN_PER_DAY*MIN_PER_DAY));

  /* set the second time derivative of mean motion */
  sgp4->n0dt2=n0dt2*(TWOPI/(MIN_PER_DAY*MIN_PER_DAY));
      
  /* bstar drag term */
  sgp4->bstar=bstar;

  /* inclination */
  sgp4->i0=DEG2RAD(i0);
      
  /* right ascension of the ascending node */
  sgp4->raan=DEG2RAD(raan);

  /* eccentricity */
  sgp4->e0=e0;

  /* argument of perigee  */
  sgp4->w0=DEG2RAD(w0);

  /* mean anomaly */
  sgp4->M0=DEG2RAD(M0);

  /* mean motion */
  sgp4->n0=n0*(TWOPI/MIN_PER_DAY);
  }


/*
 * Name:     SGP4_getPole
 * Purpose:  finds the pole (orbit plane normal)
 * Input:    double pos[3], double vel[3]
 * Output:   double pole[3]
 */
void SGP4_getPole(double pos[3], double vel[3], double pole[3]) {
  register double pNorm;

  /* find the cross product between position and velocity */
  pole[0]=pos[1]*vel[2] - vel[1]*pos[2];
  pole[1]=pos[2]*vel[0] - vel[2]*pos[0];
  pole[2]=pos[0]*vel[1] - vel[0]*pos[1];

  /* normalize the resulting pole vector */
  pNorm=sqrt(pole[0]*pole[0]+pole[1]*pole[1]+pole[2]*pole[2]);
  pole[0]/=pNorm; pole[1]/=pNorm; pole[2]/=pNorm;
  }


/*
 * Name:     SGP4_precess
 * Purpose:  precesses a direction vector from t1 to t2
 * Input:    double v[3], double t1, double t2
 * Output:   double v[3]
 */
void SGP4_precess(double v[3], double t1, double t2) {
  double t, st, a, b, c, temp[3], sina, sinb, sinc, cosa, cosb, cosc;

  /* convert times to years */
  t=0.001*(t2-t1)/365.25;
  st=0.001*(t1-MJD2000)/365.25;

  /* find the Euler angles */
  a=DEG2RAD(t*(23062.181+st*(139.656+0.0139*st)+
	       +t*(30.188-0.344*st+17.998*t)))/3600.0;
  b=DEG2RAD(t*t*(79.280+0.410*st+0.205*t)/3600.0)+a;
  c=DEG2RAD(t*(20043.109-st*(85.33+0.217*st)+
	       +t*(-42.665-0.217*st-41.833*t)))/3600.0;
  
  /* do the precession rotation */
  sina=sin(a); sinb=sin(b); sinc=sin(c);
  cosa=cos(a); cosb=cos(b); cosc=cos(c);
  temp[0]=(cosa*cosb*cosc-sina*sinb)*v[0]+(-cosa*sinb-sina*cosb*cosc)*v[1]+
    -cosb*sinc*v[2];
  temp[1]=(sina*cosb+cosa*sinb*cosc)*v[0]+(cosa*cosb-sina*sinb*cosc)*v[1]+
    -sinb*sinc*v[2];
  temp[2]=cosa*sinc*v[0]-sina*sinc*v[1]+cosc*v[2];
  
  /* replace existing direction vector */
  v[0]=temp[0]; v[1]=temp[1]; v[2]=temp[2];
  }
