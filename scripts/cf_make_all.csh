#!/usr/bin/env tcsh -f

#******************************************************************************
#*              Johns Hopkins University 
#*              Center For Astrophysical Sciences
#*              FUSE
#******************************************************************************
#*
#* Synopsis:    cf_make_all.csh
#*
#* Description: Creates 3 files from a Calfuse output (properly run):
#*
#*                1) An "all" file containing 1 extension per detector (= 8).
#*                   Each extension contains a combined spectrum from the set
#*                   of exposures, using "Xcorr" or "Background" method.
#*                2) In the case of TTAG data, creates an "ano" file. Same as
#*                   the "all" file but considering "night only" exposure time.
#*                3) A National Virtual Observatory "nvo" file. One extension
#*                   containing wavelengths that span the whole FUSE range.
#*
#*              The Xcorr method consists in co-adding spectra, the latter
#*              being corrected for a possible shift. The Background method
#*              consists in combining all the IDF files.
#*
#*              The Xcorr test is performed on 4 (detector, channel) pairs, the
#*              method of other pairs are given by them:
#*                LiF 1a -> Lif 1b
#*                LiF 2b -> LiF 2a
#*                SiC 1a -> Sic 1b
#*                Sic 2b -> Sic 2a
#*              In the script, the left segments are referred as to $seg (or 
#*              $det) and the right segments are referred as to $seg2 (or
#*              $det2).
#*  
#*              When an expected data set is missing, the script stops, cleans
#*              directory, and returns 1.
#*
#* History:     04/15/05   1.0   tc    First release
#*
#*****************************************************************************/

# Set program path
set rm = "rm -f"

set cf_xcorr = cf_xcorr
set cf_combine = cf_combine
set cf_pack = cf_pack
set cf_nvo = cf_nvo
set modhead = modhead

#set cf_xcorr = /home/vela/civeit/Work/CalFuse/Xcorr/New/cf_xcorr
#set cf_combine = /home/vela/civeit/Work/CalFuse/Shiftexp/cf_combine
#set cf_pack = /home/vela/civeit/Work/CalFuse/Pack/cf_pack
#set cf_nvo = /home/vela/civeit/Work/CalFuse/Nvo/cf_nvo
#set modhead = /home/vela/civeit/local/bin/modhead

set cf_obsplot = /data1/fuse/calfuse/v3.1/idl/cf_obsplot.pro


# Init var list
set detector = (1a 2b)
set channel = (lif sic)
set resolution = (2 3 4)
set obsmod = (hist ttag)

# Determine the root name
ls [A-Z][0-9][0-9][0-9][0-9][0-9][0-9][0-9]000asnf.fit >& tmp_asnf
if ($? == 0) then
	set asnf = `cat tmp_asnf`
	set rn = ${asnf:s/000asnf.fit//}
	$rm tmp_asnf
else
	echo "ERROR: asnf file is missing"
	$rm tmp_asnf
	goto crash
endif

# Clean tmp files that the script will create (safe)
$rm tmp_xcorr.res tmp_bpm1.lis tmp_bpm2.lis tmp_combine.lis tmp_night_exp.lis
$rm tmp_all_night_exp.lis tmp_good_exp.lis tmp_exp.lis tmp_seg_dn.lis tmp_seg_no.lis
$rm DN_${rn}*.fit NO_${rn}*.fit

foreach om ($obsmod)
	foreach res ($resolution)	
		foreach chan ($channel)	
			foreach det ($detector)
			
				# Find exposures that match the current segment
				set seg = $det$chan$res$om
				ls ${rn}[0-8][0-9][0-9]${seg}fcal.fit |& grep -v 000 > tmp_exp.lis  # Reject EXP "9xx" and "000"
			 
				if ($? == 0) then  # There are one or more exposures
				
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

					awk '{if ($5 > 0) print $6, $4, $5}' tmp_xcorr.res > tmp_all_night_exp.lis
					awk '{if ($3 >= 0) print $6,$2}' tmp_xcorr.res > tmp_good_exp.lis
					awk '{if (($3 >= 0) && ($5 > 0)) print $6,$2, $4, $5}' tmp_xcorr.res > tmp_night_exp.lis
					set time_good = `awk 'BEGIN{$t = 0} {if ($3 >= 0) $t = $t + $4} END{print $t}' tmp_xcorr.res`
					set time_bad = `awk 'BEGIN{$t = 0} {if ($3 < 0) $t = $t + $4} END{print $t}' tmp_xcorr.res`
					echo "Xcorr time: $time_good  -  Background time: $time_bad"
				
					if ($time_good > 2 * $time_bad) then
					
						#
						# --- Path 1: Optimize resolution ---
						#
					
						echo "Optimize resolution..."
						set n_good = `cat tmp_good_exp.lis | wc -l`
					
						# Extract [Day + Night] spectra
						set s = DN_${rn}000${seg}fcal.fit
						awk '{print "'$rn'"$1"'$seg'fcal.fit",$2}' tmp_good_exp.lis > tmp_combine.lis  # Combine $seg [dn]
						$cf_combine tmp_combine.lis $s -k
						# $modhead "${s}[1]" NUM_EXP $n_good
						$modhead "${s}[1]" COMBMETH XCORR
						echo $s >> tmp_seg_dn.lis
					
						set s = DN_${rn}000${seg2}fcal.fit
						awk '{print "'$rn'"$1"'$seg2'fcal.fit",$2}' tmp_good_exp.lis > tmp_combine.lis  # Combine $seg2 [dn]
						$cf_combine tmp_combine.lis $s -k
						# $modhead "${s}[1]" NUM_EXP $n_good
						$modhead "${s}[1]" COMBMETH XCORR
						echo $s >> tmp_seg_dn.lis
				
						set n_night = `cat tmp_night_exp.lis | wc -l`
						if ($om == ttag && $n_night > 0) then  # Create and combine night only files
						
							echo "*** Creating night only files ***"											
							set exp_nums = `awk '{print $1}' tmp_night_exp.lis`
						
							foreach exp ($exp_nums)
							
							    # Create night-only BPM and FCAL files if they do not already exist.
								set bno1 = NO_$rn$exp$det${om}fbpm.fit		# Bad-pixel maps
								set bno2 = NO_$rn$exp$det2${om}fbpm.fit
								set cno1 = NO_$rn$exp${seg}fcal.fit		# Extracted spectra
								set cno2 = NO_$rn$exp${seg2}fcal.fit

								set etime = `egrep "^$exp" tmp_night_exp.lis | awk '{print $3}'`
								set ntime = `egrep "^$exp" tmp_night_exp.lis | awk '{print $4}'`
								set ratio = `egrep "^$exp" tmp_night_exp.lis | awk '{printf "%.0f", 0.5+$4/$3*10.}'`

								if (!(-e $cno1)) then
								    if ($etime == $ntime) then
									echo "$cno1 is a symbolic link to $rn$exp${seg}fcal.fit"
									ln -s $rn$exp${seg}fcal.fit $cno1
								    else

								set idf_file = $rn$exp$det${om}fidf.fit
								if (!(-e $bno1)) then
								    if ($ratio > 9) then 
									echo "$bno1 is a symbolic link to $rn$exp$det${om}fbpm.fit"
									ln -s $rn$exp$det${om}fbpm.fit $bno1
								    else
									echo "Creating BPM: $bno1 ..."
									cf_bad_pixels -n $bno1 $idf_file
								    endif
								endif
							
									echo "Creating (LiF + SiC) FCAL: $cno1 ..."
									cf_extract_spectra -n $bno1 -r NO_$rn$exp $idf_file  # Existence of $bno1 is not required
								    endif
								endif

								if (!(-e $cno2)) then 
								    if ($etime == $ntime) then
									echo "Using a symbolic link to $cno2"
									ln -s $rn$exp${seg2}fcal.fit NO_$rn$exp${seg2}fcal.fit
								    else

								set idf_file = $rn$exp$det2${om}fidf.fit
								if (!(-e $bno2)) then 
								    if ($ratio > 9) then 
									echo "Using a symbolic link to $bno2"
									ln -s $rn$exp$det2${om}fbpm.fit NO_$rn$exp$det2${om}fbpm.fit
								    else
									echo "Creating BPM: $bno2 ..."
									cf_bad_pixels -n $bno2 $idf_file
								    endif
								endif
							
									echo "Creating (LiF + SiC) FCAL: $cno2 ..."
									cf_extract_spectra -n $bno2 -r NO_$rn$exp $idf_file  # Existence of $bno2 is not required
								    endif
								endif
							
							end
						
							# Combine exposures into a single spectrum.
							set s = NO_${rn}000${seg}fcal.fit
							awk '{print "NO_'$rn'"$1"'$seg'fcal.fit",$2}' tmp_night_exp.lis > tmp_combine.lis  # Combine $seg [no]
							$cf_combine tmp_combine.lis $s -k
							# $modhead "${s}[1]" NUM_EXP $n_night
							$modhead "${s}[0]" DAYNIGHT NIGHT
							$modhead "${s}[1]" COMBMETH XCORR
							echo $s >> tmp_seg_no.lis
						
							set s = NO_${rn}000${seg2}fcal.fit
							awk '{print "NO_'$rn'"$1"'$seg2'fcal.fit",$2}' tmp_night_exp.lis > tmp_combine.lis  # Combine $seg2 [no]
							$cf_combine tmp_combine.lis $s -k
							# $modhead "${s}[1]" NUM_EXP $n_night
							$modhead "${s}[0]" DAYNIGHT NIGHT
							$modhead "${s}[1]" COMBMETH XCORR
							echo $s >> tmp_seg_no.lis
						
						endif
				
					else
				
						#
						# --- Path 2: Optimize background ---
						#
					
						echo "Optimize background..."
					
						# Combine IDF files
						set idf1_all = DN_${rn}000$det${om}fidf.fit  # Same for [dn] and [no]
						set idf2_all = DN_${rn}000$det2${om}fidf.fit
										
						if (!(-e $idf1_all)) then
							echo "Creating IDF: $idf1_all ..."
							ls ${rn}[0-8][0-9][0-9]$det${om}fidf.fit |& grep -v 000 > tmp_idf.lis  # Reject EXP "9xx" and "000"
							if ($? == 0) then  # IDF files exist
							
								set idf_lis = `awk '{printf "%s ",$1}' tmp_idf.lis`
								idf_combine -c $idf1_all $idf_lis  # Create combined IDF file for $seg
								$rm tmp_idf.lis
								
							else
							
								echo "ERROR: IDF files are missing"
								$rm -f tmp_idf.lis
								goto crash
							
							endif
						endif

						if (!(-e $idf2_all)) then
							echo "Creating IDF: $idf2_all ..."
							ls ${rn}[0-8][0-9][0-9]$det2${om}fidf.fit |& grep -v 000 > tmp_idf.lis  # Reject EXP "9xx" and "000"
							if ($? == 0) then  # IDF files exist
							
								set idf_lis = `awk '{printf "%s ",$1}' tmp_idf.lis`
								idf_combine $idf2_all $idf_lis  # Create combined IDF file for $seg2
								$rm tmp_idf.lis
								
							else
							
								echo "ERROR: IDF files are missing"
								$rm -f tmp_idf.lis
								goto crash
							
							endif
						endif
					
						# Get the number of (valid) combined IDF files
						set tmp_buf = `$modhead $idf1_all NSPEC`
						set n_comb1 = $tmp_buf[3]
						set tmp_buf = `$modhead $idf1_all SPEC001`
						set idf1_1  = $tmp_buf[3]
						set tmp_buf = `$modhead $idf2_all NSPEC`
						set n_comb2 = $tmp_buf[3]
						set tmp_buf = `$modhead $idf2_all SPEC001`
						set idf2_1  = $tmp_buf[3]
					
						# Combine BPM files
						set bpm1_all = DN_${rn}000$det${om}fbpm.fit
						set bpm2_all = DN_${rn}000$det2${om}fbpm.fit
					
						if (!(-e $bpm1_all)) then 
							echo "Combine all BPM (from IDF): $bpm1_all ..."
							bpm_combine $bpm1_all $idf1_all  # Create combined BPM file [dn] for $seg
						endif
					
						if (!(-e $bpm2_all)) then 
							echo "Combine all BPM (from IDF): $bpm2_all ..."
							bpm_combine $bpm2_all $idf2_all  # Idem for $seg2
						endif
					
						# Extract [Day + Night] spectra
						set s = DN_${rn}000${seg}fcal.fit
						set xxx = DN_${rn}xxx${seg}fcal.fit
						if (!(-e $xxx)) cf_extract_spectra -r DN_${rn}xxx $idf1_all  # Avoid overwrite LiF | SiC
						echo $xxx > tmp_combine.lis
						$cf_combine tmp_combine.lis $s -k  # Does not combine, just update keywords
						$modhead "${s}[0]" NSPEC $n_comb1
						$modhead "${s}[0]" SPEC001 $idf1_1
						$modhead "${s}[1]" COMBMETH BACKGRND
						echo $s >> tmp_seg_dn.lis
					
						set s = DN_${rn}000${seg2}fcal.fit
						set xxx = DN_${rn}xxx${seg2}fcal.fit
						if (!(-e $xxx)) cf_extract_spectra -r DN_${rn}xxx $idf2_all  # Avoid overwrite LiF | SiC
						echo $xxx > tmp_combine.lis
						$cf_combine tmp_combine.lis $s -k  # Does not combine, just update keywords
						$modhead "${s}[0]" NSPEC $n_comb2
						$modhead "${s}[0]" SPEC001 $idf2_1
						$modhead "${s}[1]" COMBMETH BACKGRND
						echo $s >> tmp_seg_dn.lis
					
						set n_night = `cat tmp_all_night_exp.lis | wc -l`
						if ($om == ttag && $n_night > 0) then  # Create and combine night only files
						
							echo "*** Creating night only files ***"	
							set exp_nums = `awk '{print $1}' tmp_all_night_exp.lis`
							$rm tmp_bpm1.lis tmp_bpm2.lis
							
							foreach exp ($exp_nums)
							
								# Create bpm night only files (if they do not exist yet)
								set bno1 = NO_$rn$exp$det${om}fbpm.fit
								set bno2 = NO_$rn$exp$det2${om}fbpm.fit
							
								set idf_file = $rn$exp$det${om}fidf.fit
								if (!(-e $bno1)) then
									echo "Creating BPM: $bno1 ..."
									cf_bad_pixels -n $bno1 $idf_file
								endif
								if (-e $bno1) echo $bno1 >> tmp_bpm1.lis  # If valid, add in $seg list

								set idf_file = $rn$exp$det2${om}fidf.fit
								if (!(-e $bno2)) then 
									echo "Creating BPM: $bno2 ..."
									cf_bad_pixels -n $bno2 $idf_file
								endif
								if (-e $bno2) echo $bno2 >> tmp_bpm2.lis  # If valid, add in $seg2 list
	
							end

						
							# Combine BPM for $seg and extract spectra
							
							set bpm1_all = NO_${rn}000$det${om}fbpm.fit
							set n_bpm = `cat tmp_bpm1.lis | wc -l`
							echo "Number of valid BPM files: $n_bpm"
						
							if ($n_bpm > 0) then
						
								echo $n_bpm > tmp_bpm.lis
								cat tmp_bpm1.lis >> tmp_bpm.lis
							
								if (!(-e $bpm1_all)) then
									echo "Combine all BPM (from list): $bpm1_all ..."
									bpm_combine $bpm1_all tmp_bpm.lis  # Create combined BPM file [no] for $seg
								endif
								$rm tmp_bpm.lis
							
							endif
						
							# Extract [Night only] spectra. The existence of $bpm1_all is not required
							set s = NO_${rn}000${seg}fcal.fit
							set xxx = NO_${rn}xxx${seg}fcal.fit
							if (!(-e $xxx)) cf_extract_spectra -n $bpm1_all -r NO_${rn}xxx $idf1_all  # Avoid overwrite LiF | SiC
							echo $xxx > tmp_combine.lis
							$cf_combine tmp_combine.lis $s -k  # Does not combine, just update keywords
							$modhead "${s}[0]" NSPEC $n_comb1
							$modhead "${s}[0]" SPEC001 $idf1_1
							$modhead "${s}[1]" COMBMETH BACKGRND
							echo $s >> tmp_seg_no.lis

	
							# Combine BPM for $seg2 and extract spectra
						
							set bpm2_all = NO_${rn}000$det2${om}fbpm.fit
							set n_bpm = `cat tmp_bpm2.lis | wc -l`
							echo "Number of valid BPM files: $n_bpm"
						
							if ($n_bpm > 0) then
						
								echo $n_bpm > tmp_bpm.lis
								cat tmp_bpm2.lis >> tmp_bpm.lis
							
								if (!(-e $bpm2_all)) then 
									echo "Combine all BPM (from list): $bpm2_all ..."
									bpm_combine $bpm2_all tmp_bpm.lis  # Create combined BPM file [no] for $seg2
								endif
								$rm tmp_bpm.lis
							
							endif
						
							# Extract [Night only] spectra. The existence of $bpm2_all is not required
							set s = NO_${rn}000${seg2}fcal.fit
							set xxx = NO_${rn}xxx${seg2}fcal.fit
							if (!(-e $xxx)) cf_extract_spectra -n $bpm2_all -r NO_${rn}xxx $idf2_all  # Avoid overwrite LiF | SiC
							echo $xxx > tmp_combine.lis
							$cf_combine tmp_combine.lis $s -k  # Does not combine, just update keywords
							$modhead "${s}[0]" NSPEC $n_comb2
							$modhead "${s}[0]" SPEC001 $idf2_1
							$modhead "${s}[1]" COMBMETH BACKGRND
							echo $s >> tmp_seg_no.lis

						endif
					endif
				endif
			end
			
			$rm tmp_xcorr.res tmp_bpm1.lis tmp_bpm2.lis tmp_combine.lis
			$rm tmp_all_night_exp.lis tmp_good_exp.lis tmp_exp.lis tmp_night_exp.lis
			
		end
		
		# Pack the 8 [detector][channel] pairs together ([dn] and [no] for ttag)
		if (-e tmp_seg_dn.lis) then
	
			set fcal_all = ${rn}00000all$res${om}fcal.fit  # Final output name
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
	
				# Plot figures
				echo ".compile '$cf_obsplot'" > tmp_idl_batch
				echo "cf_obsplot, '$rn'" >> tmp_idl_batch
				echo "exit" >> tmp_idl_batch
				idl tmp_idl_batch
				$rm tmp_idl_batch
		
				# Create National Virtual Observatory file
				set nvo_file = ${rn}00000nvo$res${om}fcal.fit
				$rm $nvo_file  # Clean (safe)
				$cf_nvo $fcal_all $nvo_file  # Create file
				
				if (-e tmp_seg_no.lis) then
				
					set fcal_all = ${rn}00000ano$res${om}fcal.fit  # Final output name
					$rm $fcal_all  # Clean (safe)
					set n_segs = `cat tmp_seg_no.lis | wc -l`
	
					if (!($n_segs == 8)) then
	
						@ mseg = 8 - $n_segs
						echo "ERROR: $mseg (night only) segments are missing"
						$rm tmp_seg_no.lis
						goto crash
						
					else
					
						$cf_pack tmp_seg_no.lis $fcal_all
						$rm tmp_seg_no.lis
				
					endif
					
				endif
			endif
		endif
	end
end

# Clean [dn] and [no] files. Just keep: all, (ano) and nvo
$rm DN_${rn}*.fit NO_${rn}*.fit
exit(0)


crash:  # Procedure when script crashes

# Clean directory
$rm tmp_xcorr.res tmp_bpm1.lis tmp_bpm2.lis tmp_combine.lis tmp_night_exp.lis
$rm tmp_all_night_exp.lis tmp_good_exp.lis tmp_exp.lis tmp_seg_dn.lis tmp_seg_no.lis
$rm DN_*.fit NO_*.fit

# Return 1
exit(1)
