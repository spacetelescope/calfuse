      DOUBLE PRECISION FUNCTION sla_GMSTA (DATE, UT)
*+
*     - - - - - -
*      G M S T A
*     - - - - - -
*
*  Conversion from Universal Time to Greenwich mean sidereal time,
*  with rounding errors minimized.
*
*  double precision
*
*  Given:
*    DATE    d      UT1 date (MJD: integer part of JD-2400000.5))
*    UT      d      UT1 time (fraction of a day)
*
*  The result is the Greenwich mean sidereal time (double precision,
*  radians, in the range 0 to 2pi).
*
*  There is no restriction on how the UT is apportioned between the
*  DATE and UT arguments.  Either of the two arguments could, for
*  example, be zero and the entire date+time supplied in the other.
*  However, the routine is designed to deliver maximum accuracy when
*  the DATE argument is a whole number and the UT lies in the range
*  0 to 1 (or vice versa).
*
*  The algorithm is based on the IAU 1982 expression (see page S15 of
*  the 1984 Astronomical Almanac).  This is always described as giving
*  the GMST at 0 hours UT1.  In fact, it gives the difference between
*  the GMST and the UT, the steady 4-minutes-per-day drawing-ahead of
*  ST with respect to UT.  When whole days are ignored, the expression
*  happens to equal the GMST at 0 hours UT1 each day.
*
*  In this routine, the entire UT1 (the sum of the two arguments DATE
*  and UT) is used directly as the argument for the standard formula.
*  The UT1 is then added, but omitting whole days to conserve accuracy.
*
*  See also the routine sla_GMST, which accepts the UT as a single
*  argument.  Compared with sla_GMST, the extra numerical precision
*  delivered by the present routine is unlikely to be important in
*  an absolute sense, but may be useful when critically comparing
*  algorithms and in applications where two sidereal times close
*  together are differenced.
*
*  Called:  sla_DRANRM
*
*  P.T.Wallace   Starlink   13 April 1998
*
*  Copyright (C) 1998 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      DOUBLE PRECISION DATE,UT

*  Seconds of time to radians
      DOUBLE PRECISION S2R
      PARAMETER (S2R=7.272205216643039903848712D-5)

      DOUBLE PRECISION D1,D2,T
      DOUBLE PRECISION sla_DRANRM


*  Julian centuries since J2000.
      IF (DATE.LT.UT) THEN
         D1=DATE
         D2=UT
      ELSE
         D1=UT
         D2=DATE
      END IF
      T=(D1+(D2-51544.5D0))/36525D0

*  GMST at this UT1.
      sla_GMSTA=sla_DRANRM(S2R*(24110.54841D0+
     :                         (8640184.812866D0+
     :                         (0.093104D0
     :                         -6.2D-6*T)*T)*T
     :                         +86400D0*(MOD(D1,1D0)+MOD(D2,1D0))))

      END
