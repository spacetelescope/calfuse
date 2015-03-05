#!/usr/bin/env perl
use FileHandle;

# ******************************************************
# idlplot_corr.pl
# 
# This Perl module will run the idl script "cf_obsplot.pro"
#
# Author: Van Dixon
#
# History:  Written  July 3, 2001
#
#	10/20/05   wvd	Add call to .run cf_obsplot.pro
#	11/21/06   wvd	Add process ID ($$) to name of
#			batch file.
#	08/08/08   wvd	Add airglow argument for 900+
#			files.
#
# ******************************************************

if (@ARGV == 0) {
    print "You must enter the rootname of the observation.\n";
    print "Exiting.\n";

     } else {

     $batch_filename = $ARGV[0] . $$."_idl.bat";

     # Open the output batch file.
     open (BAT_OUTFILE, ">$batch_filename") || die "Cannot open $batch_filename";
     print BAT_OUTFILE "!path='$ENV{CF_IDLDIR}:'+!path\n";
     print BAT_OUTFILE ".run cf_obsplot.pro\n";
     if (@ARGV == 2) {
	print BAT_OUTFILE "cf_obsplot,'" . $ARGV[0] . "', airglow=1\n";
     } else {
	print BAT_OUTFILE "cf_obsplot,'" . $ARGV[0] . "'\n";
     }
     print BAT_OUTFILE "exit\n";
 
     close (BAT_OUTFILE);

     system("idl $batch_filename > /dev/null");

     system("rm $batch_filename");

     }

### end of Perl script
