lib_LTLIBRARIES = libcf.la
libcf_la_SOURCES = \
       calfits.c sgp4.c eclipse.c saa.c \
        state_limb.c state_geod.c space_vel.c helio_vel.c \
        geod_mag.c pole_ang.c solar_ang.c lsrd_vel.c lsrk_vel.c \
        month_day.c read_tle.c cf_velang.c \
        set_orbit_parms.c cf_error_msg.c \
        cf_cal_file.c cf_proc_check.c cf_proc_update.c \
        cf_timestamp.c cf_fuv_init.c cf_header_io.c \
        cf_check_digitizer.c cf_nint.c \
        cf_idf_io.c cf_ids_dead_time.c cf_electronics_dead_time.c \
        cf_fifo_dead_time.c cf_apply_dead_time.c \
        cf_thermal_distort.c cf_count_rate_y_distort.c cf_time_xy_distort.c \
        cf_geometric_distort.c cf_pha_x_distort.c \
        cf_active_region.c cf_find_spectra.c cf_identify_channel.c \
        cf_calculate_ycent_motion.c cf_source_aper.c\
        cf_grating_motion.c cf_fpa_position.c cf_read_fpa_pos.c \
        cf_make_mask.c cf_mirror_motion.c \
        cf_satellite_jitter.c cf_calculate_y_centroid.c \
        cf_target_count_rate.c \
        cf_screen_jitter.c cf_screen_limb_angle.c cf_screen_saa.c \
        cf_screen_high_voltage.c cf_screen_burst.c cf_screen_airglow.c \
        cf_screen_bad_pixels.c cf_set_user_gtis.c \
        cf_set_photon_flags.c cf_set_good_time_intervals.c \
        cf_modify_hist_times.c cf_screen_pulse_height.c \
        cf_convert_to_ergs.c cf_extraction_limits.c \
        cf_astigmatism.c cf_dispersion.c cf_doppler_and_heliocentric.c \
        cf_apply_filters.c cf_scale_bkgd.c \
        cf_make_wave_array.c cf_rebin_background.c \
        cf_rebin_probability_array.c cf_optimal_extraction.c \
        cf_write_extracted_spectrum.c cf_standard_or_optimal_extraction.c \
        cf_init_support.c cf_modify_hist_pha.c \
        cf_fes_proc_check.c cf_fes_proc_update.c
libcf_la_LDFLAGS = -lm -lgfortran
libcf_la_CFLAGS = $(cfitsio_CFLAGS) -I$(srcdir)/include -DCFORTRAN -Df2cFortran
libcf_la_LIBADD = $(cfitsio_LIBS) ../slalib/libslarandom.la ../slalib/libsla.la
#AM_CFLAGS = -I$(srcdir)/include
#AM_CFLAGS = -I$(pkgincludedir)/include

