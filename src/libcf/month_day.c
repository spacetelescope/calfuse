/* month_day: return month,day from day of year */
/*	07/21/04    1.4   wvd	Add () to defn of leap */

#include "calfuse.h"

void month_day(int year, int yearday, int *pmonth, int *pday)
{
    static char daytab[2][13]={
	{0,31,28,31,30,31,30,31,31,30,31,30,31},
	{0,31,29,31,30,31,30,31,31,30,31,30,31}
    };
    int i, leap;

    leap = ((year%4 == 0 && year%100 != 0) || year%400 == 0);
    for (i=1; yearday>daytab[leap][i]; i++)
	yearday-=daytab[leap][i];
    *pmonth=i;
    *pday=yearday;
}
