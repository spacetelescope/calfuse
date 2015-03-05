/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    int saa(saareg *saa, double lon, double lat)
 *
 * Description: Determines if FUSE is in the SAA for given models of the SAA
 *              
 *
 * Arguments:   saareg *saa           A structure containing data that roughly
 *                                    outlines the SAA region.
 *                                    (n_points: number of points
 *                                     lat:    the latitude array
 *                                     lon: the longitude array)
 *                                    (currently we have only one model).
 *              double lon,lat        Geocentric longitude and latitude of FUSE
 *
 * Returns:     int                   0 if not in SAA
 *                                    1 if in SAA
 *                    
 * History:     07/17/98  E. Murphy   Begin work.
 *              06/22/99        peb   Tidied code, added FITS_ wrappers,
 *                                    removed hardcoded variables.
 *              07/02/99        peb   Moved SAA data to a structure.
 *              08/25/99        emm   Fixed bug in setting maxlon
 *              08/28/03 v1.3   bjg   Adopt Bryce's algorithm for determining
 *					when FUSE is in the SAA.
 *              08/28/03 v1.4   bjg   v1.3 didn't compile but v1.4 does
 *              08/29/03 v1.5   bjg   slightly modified getDir
 *
 ****************************************************************************/
 
#include "calfuse.h"

int getDir(double lon0,double lat0,double lon1,double lat1,double lon2,double lat2)
{

double dx1, dx2, dy1, dy2, t1, t2;
    
  dx1=lon1-lon0; dy1=lat1-lat0;
  dx2=lon2-lon0; dy2=lat2-lat0;

  t1=dx1*dy2; t2=dy1*dx2;

  if (t1>t2) 
    return(-1);
  if (t1<t2) 
    return(1);
  

  return(1);

}


int
saa(saareg *saa, double lon, double lat)
{
int i;
int dir=0;
int oldDir=0;

for (i=0; i<saa->n_points-1; i++) {
    dir=getDir(lon,lat,(double)(saa->lon[i]),(double)(saa->lat[i]),(double)(saa->lon[i+1]),(double)(saa->lat[i+1]));
    if (dir*oldDir<0)
      return(0);
    else
      oldDir=dir;
    }
    
  return(1);
}
