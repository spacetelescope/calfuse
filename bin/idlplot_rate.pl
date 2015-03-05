#!/usr/bin/env perl
use FileHandle;

# ***************************************************
# add_tle.pl
# 
# This Perl module will read in the latest five orbital elements
# from the file five.tle (which was created by get_tle.pl) and
# add any new orbital elements to the file FUSE.TLE.  The file
# FUSE.TLE is in descending order (i.e. the most recent elements
# are first).  In order to prepend the new TLE onto the old list
# I found it was easiest to store everything in a temporary file
# TEMP.TLE and rewrite FUSE.TLE.
#
# Author: Ed Murphy
#
# History:  Written  July 27, 1999
#	11/21/06   wvd	Add process ID ($$) to name of
#			batch file.
#
# ***************************************************

# Define the file names.  old_maintle_filename is not actually opened,
# but is used in a system call at the end of the program. 

if (@ARGV == 0) {
    print "You must enter the rootname of the observation.\n";
    print "Exiting.\n";

     } else {

     $batch_filename = $ARGV[0] . $$."_idl.bat";

     # Open the output batch file.
     open (BAT_OUTFILE, ">$batch_filename") || die "Cannot open $batch_filename";
     print BAT_OUTFILE "!path='$ENV{CF_IDLDIR}:'+!path\n";
     print BAT_OUTFILE "cf_plot_rate3,'" . $ARGV[0] . "'\n";
     print BAT_OUTFILE "exit\n";
 
     close (BAT_OUTFILE);

     system("idl $batch_filename > /dev/null");

     system("rm $batch_filename");

     }

### end of Perl script
