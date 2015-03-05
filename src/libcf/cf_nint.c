/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    n = cf_nint(x)
 *
 * Description: Converts double to nearest integer.
 *
 * Variables:   double    x	A double
 *
 * Return:      int       n	Nearest integer
 *
 * History:     08/25/03    wvd  v1.1   Begin work
 *              02/09/04    wvd  v1.2   Add cf_nlong()
 *              02/17/04    wvd  v1.3   Test for overflow of int, long
 *              02/18/04    wvd  v1.4   Change format of error message.
 *					Change from values.h to limits.h
 *
 *************************************************************************/

#include <limits.h>
#include "calfuse.h"

int
cf_nint(double x)

{
	if (x > INT_MAX || x < INT_MIN)
	    cf_if_error("Cannot convert %10.4e to an integer", x);

	if (x < 0.)	return (int) (x - 0.5);
	else		return (int) (x + 0.5);
}


long
cf_nlong(double x)

{
	if (x > LONG_MAX || x < LONG_MIN)
	    cf_if_error("Cannot convert %10.4e to a long", x);

        if (x < 0.)     return (long) (x - 0.5);
        else            return (long) (x + 0.5);
}
