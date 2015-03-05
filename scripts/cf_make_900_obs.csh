#!/usr/bin/env tcsh -f

#******************************************************************************
#*              Johns Hopkins University 
#*              Center For Astrophysical Sciences
#*              FUSE
#******************************************************************************
#*
#* Synopsis:    cf_make_900_obs.csh association_file
#*
#* Description: This routine is derived from cf_make_all_obs.csh, but modified
#*		to consider only 900+ (airglow) files.  Output is a single
#*		quick-look image of the airglow spactrum.
#*
#*  		Extracted spectra are combined, not IDF files.
#*
#*              When an expected data set is missing, the script stops, cleans
#*              the directory, and returns 1.
#*
#* History:     08/08/08   1.00	 wvd   Create separate quick-look image
#*					for 900+ (airglow) exposures.
#*					Always combine extracted spectra.  If
#*					cf_xcorr fails, combine with no shift.
#* 		08/08/08   1.01	 wvd   Rename airglow quick-look file to
#*					M112580100000airgttagf.gif
#*
#*****************************************************************************/

# Delete files after processing?  (Default is no.)
#set DELETE_IDF                 # Delete intermediate data files
set DELETE_BPM                 # Delete bad-pixel map files

# Set program path
set rm = "rm -f"

set cf_xcorr = cf_xcorr
set cf_combine = cf_combine
set cf_pack = cf_pack
set cf_nvo = cf_nvo
set idl_obsplot = idl_obsplot.pl
set modhead = modhead

# Init var list
set detector = (1a 2b)
set channel = (lif sic)
set resolution = (2 3 4)
set obsmod = (hist ttag)

# Determine the root name and the program ID
set asnf = $1
set rn = ${asnf:s/000asnf.fit//}
set pid = `echo $rn | awk '{print substr($1, 1, 4)}'`

# Determine the object class
set tmp_file = `ls ${rn}*fcal.fit | awk '{if (NR == 1) print}'`
set tmp_buf = `$modhead $tmp_file OBJCLASS`
set objclass = $tmp_buf[2]

# Clean tmp files that the script will create (safe)
$rm tmp_xcorr.res tmp_combine.lis 
$rm tmp_good_exp.lis tmp_exp.lis tmp_seg_dn.lis 
$rm DN_${rn}*.fit

foreach om ($obsmod)
    foreach res ($resolution)    
        foreach chan ($channel)    
            foreach det ($detector)
            	
		set ignore_exp_stat = ''
		if ($om == hist) set ignore_exp_stat = -a
		if ($pid == S100) set ignore_exp_stat = -a
		if ($pid == M106) set ignore_exp_stat = -a
		
		# Find exposures that match the current segment
                set seg = $det$chan$res$om
                set readfiles = 0
                ls ${rn}9[0-9][0-9]${seg}fcal.fit |& grep -v 000${seg} > tmp_exp.lis  # Keep only "9xx" exposures
                if ($? == 0) then 
                	set readfiles = 2
			set ignore_exp_stat = -a
                endif
                
                if ($readfiles >= 1) then  # There are one or more exposures
                
                    echo " "
                    echo "*** Processing: $seg ***"
                    # [1a][lif] -> [1b][lif], [2b][lif] -> [2a][lif] etc...
                    if ($det == 1a) set det2 = 1b
                    if ($det == 2b) set det2 = 2a
                    set seg2 = $det2$chan$res$om
                
                    echo "----- cf_xcorr input   -----"
                    cat tmp_exp.lis
                    echo "----------------------------"
                    $cf_xcorr tmp_exp.lis tmp_xcorr.res  # Compute shift and sigma_shift
                    echo "----- cf_xcorr results -----"
                    cat tmp_xcorr.res
                    echo "----------------------------"

                    awk '{if ($3 >= 0) {print $6, $2} else {print $6}}' tmp_xcorr.res > tmp_good_exp.lis
                    
                        #
                        # --- Path 1: Optimize resolution ---
                        #
                    
                        echo "Optimize resolution..."
                        set n_good = `cat tmp_good_exp.lis | wc -l`
                    
                        # Extract [Day + Night] spectra
                        set s = DN_${rn}900${seg}fcal.fit
                        awk '{print "'$rn'"$1"'$seg'fcal.fit",$2}' tmp_good_exp.lis > tmp_combine.lis  # Combine $seg [dn]
			# echo "----- Combining Files ----- "
			# cat tmp_combine.lis
			# echo "Output: " $s
                        $cf_combine -k $ignore_exp_stat tmp_combine.lis $s
                        # $modhead "${s}[1]" NUM_EXP $n_good
                        $modhead "${s}[1]" COMBMETH XCORR
                        echo $s >> tmp_seg_dn.lis
                    
                        set s = DN_${rn}900${seg2}fcal.fit
                        awk '{print "'$rn'"$1"'$seg2'fcal.fit",$2}' tmp_good_exp.lis > tmp_combine.lis  # Combine $seg2 [dn]
			# echo "----- Combining Files ----- "
			# cat tmp_combine.lis
			# echo "Output: " $s
                        $cf_combine -k $ignore_exp_stat tmp_combine.lis $s
                        # $modhead "${s}[1]" NUM_EXP $n_good
                        $modhead "${s}[1]" COMBMETH XCORR
                        echo $s >> tmp_seg_dn.lis
                
                endif
            end
            
            $rm tmp_xcorr.res tmp_combine.lis
            $rm tmp_good_exp.lis tmp_exp.lis 
            
        end
        
        # Pack the 8 [detector][channel] pairs together
        if (-e tmp_seg_dn.lis) then
    
            set fcal_all = ${rn}00900all$res${om}fcal.fit  # Final output name
            $rm $fcal_all  # Clean (safe)
            set n_segs = `cat tmp_seg_dn.lis | wc -l`
    
            if (!($n_segs == 8)) then
    
                @ mseg = 8 - $n_segs
                echo "ERROR: $mseg (day + night) segments are missing"
                $rm tmp_seg_dn.lis
                goto crash

            else
    
                $cf_pack tmp_seg_dn.lis $fcal_all
                $rm tmp_seg_dn.lis
    
                # Plot figures, delete unwanted files
                $idl_obsplot {$rn} airglow
		mv ${rn}00900spec${om}f.gif ${rn}00000airg${om}f.gif
		$rm ${rn}00900lif*.gif ${rn}00900sic*.gif
		$rm $fcal_all
        
            endif
        endif
    end
end

# Clean [dn] files. 
$rm DN_${rn}*.fit

# Delete IDF files
if $?DELETE_IDF then
    echo "NOTE: Deleting intermediate data files."
    $rm ${rn}*idf.fit
endif

# Delete bad-pixel-map (bpm) files
if $?DELETE_BPM then
    echo "NOTE: Deleting bad pixel map (bpm) files."
    $rm ${rn}*bpm.fit
endif

exit(0)

crash:  # Procedure when script crashes

# Clean directory
$rm tmp_xcorr.res tmp_combine.lis 
$rm tmp_good_exp.lis tmp_exp.lis tmp_seg_dn.lis 
$rm DN_*.fit

# Return 1
exit(1)
