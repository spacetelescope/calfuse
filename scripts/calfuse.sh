#!/bin/sh
#******************************************************************************
#*              Johns Hopkins University 
#*              Center For Astrophysical Sciences
#*              FUSE
#******************************************************************************
#*
#* Synopsis:    calfuse.sh file_name 
#*
#* Description: Shell script for processing FUSE Ver 3.2 time-tagged exposures.
#*		All messages are written to stdout or stderr.
#*
#* Arguments:   char    file_name	File name to process is 1st command-
#*                                      line argument
#*
#* Returns:     Exit codes:
#*			0		Successful execution
#*			
#*
#* Environment variables:  CF_CALDIR    Path to directory which contains the
#*                                      calibration files.
#*                         CF_IDLDIR    Path to CalFUSE IDL directory
#*                         CF_PIPELINE  Flag to determine if CALFUSE
#*                                      is running as part of the JHU 
#*                                      pipeline.
#*
#* History:     09/05/07   1.1   bot    Adapted from tcsh script
#*
#*****************************************************************************/

# Delete files after processing?  (Default is no.)
#DELETE_IDF=1	  		# Delete intermediate data file
#DELETE_BPM=1			# Delete bad-pixel map

idf=`echo $1 | sed -e "s/raw/idf/g"`
froot=`echo $1 | sed -e "s/raw.fit//g"`
logfile=${froot}.trl
ttag=`echo $froot | grep -c ttag`


# Put a timestamp in the log file (the OPUS trailer file).
if [ $ttag = 1 ]; then
    echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: Begin TTAG file $1"
    echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: Begin TTAG file $1" >> $logfile 2>&1
else
    echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: Begin HIST file $1"
    echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: Begin HIST file $1" >> $logfile 2>&1
fi

cfstat=1

# Step 1  --  Generate Intermediate Data File
if [ $ttag = 1 ]; then
    cf_ttag_init        $1    $idf    >> $logfile 2>&1
    cfstat=$?
else
    cf_hist_init        $1    $idf    >> $logfile 2>&1
    cfstat=$?
fi

# Step 2  --  Convert to FARF
if [ $cfstat = 0 ]; then
    cf_convert_to_farf        $idf    >> $logfile 2>&1
    cfstat=$?
fi

# Step 3  --  Screen photons
if [ $cfstat = 0 ]; then
    cf_screen_photons         $idf    >> $logfile 2>&1
    cfstat=$?
fi

# Step 4  --  Remove motions
if [ $cfstat = 0 ]; then
    cf_remove_motions         $idf    >> $logfile 2>&1
    cfstat=$?
fi

if [ $cfstat = 0 ]; then
	if [ ${CF_IDLDIR:-""} != "" ]; then
	    idlplot_rate.pl   $froot  >> $logfile 2>&1
            idlplot_spex.pl   $froot  >> $logfile 2>&1
	fi
fi

# Step 5  --  Assign wavelength
if [ $cfstat = 0 ]; then
    cf_assign_wavelength      $idf    >> $logfile 2>&1
    cfstat=$?
fi

# Step 6  --  Flux calibrate
if [ $cfstat = 0 ]; then
    cf_flux_calibrate         $idf    >> $logfile 2>&1
    cfstat=$?
fi

# Step 7 -- Create a bad-pixel file
if [ $cfstat = 0 ]; then
    cf_bad_pixels             $idf    >> $logfile 2>&1
    cfstat=$?
fi

# Step 8  --  Extract spectra
if [ $cfstat = 0 ]; then
    cf_extract_spectra        $idf    >> $logfile 2>&1
    cfstat=$? 
fi

# Step 8a --  Delete _bursts.dat file
if [ ${CF_PIPELINE:-""} = "" ]; then
    rm -f `echo $1 | sed -e "s/ttagfraw.fit/_bursts.dat/g"`
fi

# Step 8b --  Delete IDF file
if [ ${DELETE_IDF:-""} != "" ]; then
    echo "NOTE: Deleting intermediate data file."
    rm -f $idf
fi

# Step 8c --  Delete bad pixel map (bpm) file
if [ ${DELETE_BPM:-""} != "" ]; then
    echo "NOTE: Deleting bad pixel map (bpm) file."
    rm -f `echo $1 | sed -e "s/raw/bpm/g"`
fi

if [ $cfstat = 0 ]; then
 if [ $ttag = 1 ]; then
  echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: End   TTAG file $1"
  echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: End   TTAG file $1"  >> $logfile 2>&1
 else
  echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: End   HIST file $1"
  echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: End   HIST file $1"  >> $logfile 2>&1
 fi
else
  echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: Error processing $1"
  echo `date '+%Y %b %e %T'` "calfuse.sh-1.15: Error processing $1" >> $logfile 2>&1
fi

exit $cfstat

#******************************************************************************
