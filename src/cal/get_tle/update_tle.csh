#!/bin/csh -f
#******************************************************************************
#*              Johns Hopkins University
#*              Center For Astrophysical Sciences
#*              FUSE
#*****************************************************************************
#*
#* Synopsis:    update_tle
#*
#* Description: Shell script for automatically downloading the latest
#*              orbital elements from the GSFC OIG, and placing these
#*              elements in the FUSE.TLE file.  It will also update 
#*              the cvz_ram_tool.html calculator.
#*		
#*		All messages are to stdout or stderr.
#*
#* Arguments:   None
#*
#* Returns:     Exit codes:
#*			0		successful execution
#* 
#* History:     07/27/99        emm     Begin work.
#*              10/23/00        mlr     changed path from v1.4 to current
#******************************************************************************/

set tlestat=0

# Step 1
#/usr/local/fusesw/calfuse/v1.4/src/cal/get_tle/get_tle.pl
#/usr/local/fusesw/calfuse/current/src/cal/get_tle/get_tle.pl
#/usr/local/j2sdk1.4.2_06/bin/java -classpath /usr/local/fusesw/calfuse/current/src/cal/get_tle GetTLE
/usr/local/java/bin/java -classpath /usr/local/fusesw/calfuse/current/src/cal/get_tle GetTLE
set cfstat=$status

# Step 2
if !({$cfstat}) then
#    /usr/local/fusesw/calfuse/v1.4/src/cal/get_tle/add_tle.pl
    /usr/local/fusesw/calfuse/current/src/cal/get_tle/add_tle.pl
    set cfstat=$status
endif

# Step 3
if !({$cfstat}) then
#    /usr/local/fusesw/calfuse/v1.4/src/cal/get_tle/make_cvzramtool.pl
    /usr/local/fusesw/calfuse/current/src/cal/get_tle/make_cvzramtool.pl
    set cfstat=$status
endif

# Step 4
if !({$cfstat}) then
#    /usr/local/fusesw/calfuse/v1.4/src/cal/get_tle/make_orbit.pl
    /usr/local/fusesw/calfuse/current/src/cal/get_tle/make_orbit.pl
    set cfstat=$status
endif

exit($cfstat)

