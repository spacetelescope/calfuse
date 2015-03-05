/*****************************************************************************
 *              Johns Hopkins University
 *              Center For Astrophysical Sciences
 *              FUSE
 *****************************************************************************
 *
 * Synopsis:    #include "calfuse.h"
 *
 * Description: Master include file for FUSE calibration pipeline processing
 *		system. Global shared variables and structures are defined
 *              here, as well as symbols and default values of parameters.
 *
 * Arguments:   none
 *
 * Returns:     none
 *
 * History:     11/01/02   1.1   peb    Begin work
 *		12/10/02   1.3   wvd    Continue work
 *		12/20/02   1.4   wvd    Change flags to unsigned char in 
 *					screening routines.
 *		01/14/03   1.5   wvd    Change cf_fpa_pos to cf_read_fpa_pos
 *		01/17/03   1.6   wvd    Added new calibration files
 *					 DIGI_CAL and SPEC_CAL
 *					Added cf_check_digitizer and DIGI_COR
 *		02/05/03   1.7   wvd    Added spectral extraction routines
 *					Interpolate between FLUX_CAL files.
 *					Define HC
 *		02/12/03   1.8   wvd    Added convert_to_ergs
 *					Replace FLUX_CAL with AEFF_CAL
 *		02/14/03   1.9   wvd    Added FLAT, WGTS and WORM calibration
 *					file keywords.
 *		02/24/03   1.10  wvd    Added subroutines of cf_extract
 *					and char *cf_hist_file(char *);
 *              02/28/03   1.12  peb    Added function definitions for 
 *                                      cf_rebin_and_flux_calibrate_backround,
 *                                      cf_standard_or_optimal_extraction,
 *                                      cf_optimal_extraction, and
 *                                      cf_write_extracted_spectrum.
 *              03/04/03   1.13  peb    Added astig_read_file,
 *                                      astig_check_input_image,
 *                                      astig_target_aperture
 *              03/05/03   1.14  wvd    Added cf_target_count_rate
 *              03/07/03   1.15  peb    Changed cf_error_init, added
 *                                      cf_verbose, and made pha, timeflgs, and
 *                                      loc_flgs function type consistent.
 *              03/11/03   1.16  wvd    Changed channel from unsigned char to
 *					char in a few subroutines.
 *              03/18/03   1.17  wvd    variable argument list for cf_verbose 
 *              03/25/03   1.18  peb    variable argument list for cf_if_error
 *                                      and cf_if_warning
 *		03/21/03   1.19  wvd	Add flag for photon in pot hole
 *		04/01/03   1.21  wvd	Delete cf_errmsg; obsolete
 *		04/08/03   1.22  wvd	Change definition of cf_apply_filters
 *		04/09/03   1.23  wvd	Add bad-pixel maps to calibration files.
 *		04/17/03   1.25  wvd	Add final_call to cf_identify_channel
 *					Add final_call and weight to
 *					cf_calculate_y_centroid
 *		04/18/03   1.26  wvd	Add cf_find_spectra
 *		04/21/03   1.27  wvd	Modify cf_grating_motion, 
 *					cf_mirror_motion, and
 *					cf_satellite_jitter.
 *					Define FRAME_TOLERANCE.
 *              04/28/03  v1.28  wvd    Modified cf_fuv_init not to extrapolate
 *                                      the last two calibration files forward.
 *              05/10/03  v1.29  wvd    Pass locflag to cf_set_photon_flags
 *              05/16/03  v1.30  wvd    Add cf_make_mask.
 *					Set PERFORM keywords for HIST data.
 *              05/28/03  v1.31  rdr    Modified def of cf_optimal_extraction
 *                                      and cf_write_extracted_spectrum
 *              05/30/03  v1.34  wvd    Pass WEIGHTS to cf_set_photon_flags
 *              06/09/03  v1.35  rdr    Change definition of cf_apply_filters
 *              06/11/03  v1.36  wvd	Change HV array to type short.
 *              06/11/03  v1.37  wvd	Pass datatype to cf_read_col and
 *					cf_write_col.
 *              07/16/03  v1.38  wvd	Move initialization routines to libcf
 *              07/23/03  v1.39  wvd	Add HSKP_CAL to list of cal files
 *					Increment OPUS_VERSION to 2.7
 *              08/01/03  v1.40  wvd	Add cf_apply_dead_time to pipeline,
 *					modify arguments of dead-time routines.
 *              08/04/03  v1.41  wvd	Convert count-rate arrays to shorts.
 *              08/06/03  v1.42  wvd	Delete GTI's from cf_satellite_jitter
 *              08/21/03  v1.43  wvd	Change channel array to unsigned char
 *					in subroutines of cf_remove_motion
 *              08/22/03  v1.44  wvd	Add cf_get_extraction_limits.
 *              08/22/03  v1.45  wvd	Change limits of extraction window
 *					from type int to type short.
 *              08/25/03  v1.46  wvd	Add cf_nint, change coltype in 
 *					cf_idf_io routines from string to int
 *              08/28/03  v1.47  wvd	Modify structure saareg
 *              09/10/03  v1.48  wvd	Define cf_set_user_gtis
 *					Modify args of cf_screen_burst
 *              09/15/03  v1.49  wvd	Add structure top_level_routines
 *					Delete astig_read_file and
 *					astig_target_aperture.
 *					Modify args of cf_identify_channel
 *              10/02/03  v1.50  wvd	Exchange PHA and GTI flags.
 *					Change version number to 3.0.2.
 *              10/08/03  v1.51  wvd	Change counts_out to type long
 *              10/26/03  v1.52  wvd	Change arguments of cf_find_spectra
 *					and cf_calculate_y_centroid.
 *              10/31/03  v1.53  wvd	Change channel to unsigned char
 *					throughout.
 *              11/26/03  v1.54  wvd	Change aic_rate and fec_rate to float
 *					throughout.
 *              12/21/03  v1.55  wvd	Remove underscore from idf and bpm
 *					filenames.
 *                                      Change version number to 3.0.3.
 *              02/09/04  v1.56  wvd	Employ new scheme for dead-time 
 *					correction.  Add cf_nlong() and
 *					cf_screen_fifo_overflow.
 *                                      Change version number to 3.0.4.
 *              02/27/04  v1.57  rdr    Change def of cf_thermal_distortion
 *              03/02/04  v1.58  wvd	Implement WPC array in extraction
 *					routine.
 *					Add cf_x2lambda.
 *					Change version number to 3.0.5.
 *              03/16/04  v1.59  wvd	Delete WPC array.
 *					Smooth HIST data in X.
 *					Comment out cf_astig_farf.
 *					Change version number to 3.0.6.
 *              04/05/04  v1.60  wvd	Modify cf_geometric_distort to 
 *					rescale SPECBINY only for HIST data.
 *              04/09/04  v1.61  wvd	Fix bugs in cf_optimal_extraction.
 *					Change version number to 3.0.7.
 *              04/09/04  v1.62  bjg    Define FILL_DATA and LOCATION_FILL 
 *              04/26/04  v1.63  wvd    Replace cf_rebin_and_flux_
 *					calibrate_background with
 *					cf_rebin_background.
 *					Modify args to cf_optimal_extraction
 *					and cf_find_spectra.
 *              06/02/04  v1.64  wvd    Add cf_modify_hist_times.
 *					Populate the header keywords
 *					TIME_COR, COMB_COR, and QUIK_COR
 *					Modify args to cf_calculate_y_centroid,
 *					cf_find_spectra, cf_satellite_jitter,
 *					cf_apply_filters, and
 *					cf_write_extracted_spectrum.
 *					Modify order of CALIBRATION_STEP_KEYS.
 *              08/19/04  v1.65  wvd    Add FES definitions and subroutines.
 *              10/12/04  v1.66  wvd    Change version number to 3.0.8
 *              02/01/05  v1.67  wvd    Change version number to 3.1.0
 *					Modify args to cf_screen_burst
 *              03/02/05  v1.68  wvd    Add cf_modify_hist_pha and PHAH_COR.
 *					Walk correct HIST data.
 *					Change cf_ttag_bkgd to cf_scale_bkgd
 *					and pass weights array to it.
 *					Change cf_get_extraction_limits to
 *					cf_extraction_limits; it now returns
 *					X limits of extraction window.
 *					Add cf_screen_airglow.
 *              03/22/05  v1.69  wvd    Change TIME_SUNRISE and TIME_SUNSET
 *					from floats to shorts.
 *              04/19/05  v1.70  wvd    Change version number to 3.1.1
 *              06/15/05  v1.71  wvd    BUG FIX: cf_extract_spectra always
 *					read the point-source probability
 *					array from WGTS_CAL file.  Now uses
 *					variable "extended" to determine
 *					which HDU to read.  Modify args to 
 *					cf_rebin_probability_array
 *					Change version number to 3.1.2
 *					Delete QUIK_COR from structure
 *					CALIBRATION_STEP_KEYS.
 *              08/30/05  v1.72  wvd    Define MAX_EXPTIME = 55000
 *					Delete cf_read_fpa_pos
 *					Change version number to 3.1.3
 *              09/19/05  v1.73  wvd    Reinstall cf_read_fpa_pos, as it is
 *					needed by ttag_combine.
 *              09/30/05  v1.74  wvd    Change version number to 3.1.4
 *					Pass photon array to
 *					cf_screen_fifo_overflow.
 *              11/22/05  v1.75  wvd    Add cf_screen_bad_pixels, 
 *					cf_screen_jitter, and
 *					cf_get_potholes
 *					Change version number to 3.1.5
 *              01/24/06  v1.76  wvd    Change version number to 3.1.6
 *              02/03/06  v1.77  wvd    Change version number to 3.1.7
 *              05/15/06  v1.78  wvd    Divide cf_astigmatism_and_dispersion
 *					into two separate routines.  Incorporate
 *					cf_x2lambda into cf_dispersion.
 *					Change version number to 3.1.8
 *					Delete cf_astig_farf.
 *              06/12/06  v1.79  wvd    Add pole_ang.c
 *					Change version number to 3.1.9
 *					Add -a to cf_remove_motions
 *              11/02/06  v1.80  wvd    Add cf_time_xy_distort.c
 *					Change version number to 3.2.0
 *					Add APER_COR to list of cal steps.
 *					Change cf_screen_fifo_overflow to
 *					cf_fifo_deadtime.  Modify args to
 *					it, cf_apply_dead_time, and
 *					cf_target_count_rate.  Run
 *					cf_target_count_rate on HIST data.
 *              03/07/07  v1.81  wvd    Modify arguments to space_vel.
 *              05/18/07  v1.82  wvd    Change version number to 3.2.1.
 *              09/15/07  v1.83  bot    Change version number to 3.2.2.
 *              10/16/07  v1.84  bot    Added brackets in 
 *					FES_CALIBRATION_STEP_KEYS
 *					and in CALIBRATION_FILE_KEYS.
 *              08/22/08  v1.85  wvd    Change version number to 3.2.3.
 *					Many changes to better handle
 *					bright-earth observations & 
 *					900-level airglow exposures.
 *
 ****************************************************************************/

#include "calfitsio.h"

#define CALFUSE_VERSION "3.2.3"

#define LARGEMJD        9999999999.0
#define OPUS_VERSION    2.7      /* Oldest compatible version of OPUS */
#define PI 3.1415926535897932384626433832795028841971693993751
#define RADIAN 0.017453292519943295769236907684886127134428718885417
#define C 299792.458
#define HC 1.98644746104e-8	/* erg A */
#define MU 3.986005E5		/* km^3 s^-2 */
#define RE 6371.00		/* km */
#define RS 6.960E5
#define AU 1.495978707E8
#define FRAME_TOLERANCE 0.004

#define FESPIX          266256          /* This is 516*516 */
#define FILL_DATA       21865
#define NXMAX           16384
#define NYMAX           1024

#define TEMPORAL_DAY      (0x01<<0)
#define TEMPORAL_LIMB     (0x01<<1)
#define TEMPORAL_SAA      (0x01<<2)
#define TEMPORAL_HV       (0x01<<3)
#define TEMPORAL_BRST     (0x01<<4)
#define TEMPORAL_OPUS     (0x01<<5)
#define TEMPORAL_JITR     (0x01<<6)
#define TEMPORAL_USER     (0x01<<7)

#define LOCATION_SHLD     (0x01<<0)
#define LOCATION_AIR      (0x01<<1)
#define LOCATION_STIML    (0x01<<2)
#define LOCATION_STIMR    (0x01<<3)
#define LOCATION_PHA      (0x01<<4)
#define LOCATION_BADPX    (0x01<<5)
#define LOCATION_FILL     (0x01<<6)

#define MAX_EXPTIME	55000

struct fes_keyword_tab
{
       char  name[9];
       char  value[8];
       char  proc[18];
};

#define NUM_FES_PROC_STEPS 6

#define FES_CALIBRATION_STEP_KEYS { \
  {"INIT_FES\0","PERFORM\0","cf_fes_init\0"},\
  {"MASK_FES\0","PERFORM\0","cf_fes_mask\0"},\
  {"BIAS_FES\0","PERFORM\0","cf_fes_bias\0"},\
  {"FLAT_FES\0","PERFORM\0","cf_fes_flat\0"},\
  {"UNDS_FES\0","PERFORM\0","cf_fes_undistort\0"},\
  {"FLUX_FES\0","PERFORM\0","cf_fes_flux\0"},\
}

#define NUM_FES_CAL_KEYS 5

#define FES_CALIBRATION_FILE_KEYS { \
      "MASK",1,"FCL","\0","\0","\0",0.0,0.0,LARGEMJD,0,0,0,\
      "BIAS",2,"FCL","\0","\0","\0",0.0,0.0,LARGEMJD,0,0,0,\
      "FFLT",2,"FCL","\0","\0","\0",0.0,0.0,LARGEMJD,0,0,0,\
      "FGEO",1,"FCL","\0","\0","\0",0.0,0.0,LARGEMJD,0,0,0,\
      "FFLX",2,"FCL","\0","\0","\0",0.0,0.0,LARGEMJD,0,0,0\
    }

typedef struct {
    double ra_ap;
    double dec_ap;
    double limb;
} orbital;

typedef struct {
    long   n_points;
    float *lat;
    float *lon;
} saareg;

typedef struct {
    long ntimes;        /*  The number of intervals  */
    double *start;      /*  An array of starting times (in seconds)  */
    double *stop;       /*  An array of stoping times (in seconds)  */
} GTI;                  /*  Good Time Intervals  */

struct keyword_tab
{
    char  name[9];
    char  hist_value[8];
    char  hist_proc[32];
    char  ttag_value[8];
    char  ttag_proc[32];
};

#define NTOP_LEVEL_ROUTINES 11

#define TOP_LEVEL_ROUTINES { \
    "cf_ttag_init", \
    "cf_hist_init", \
    "cf_convert_to_farf", \
    "cf_ttag_countmap", \
    "cf_ttag_gainmap", \
    "cf_remove_motions", \
    "cf_assign_wavelength", \
    "cf_screen_photons", \
    "cf_flux_calibrate", \
    "cf_bad_pixels", \
    "cf_extract_spectra" \
}

#define NUM_PROC_STEPS 40

#define CALIBRATION_STEP_KEYS { \
    {"INIT_COR", "PERFORM", "cf_hist_init", "PERFORM", "cf_ttag_init"}, \
    {"DIGI_COR", "PERFORM", "cf_check_digitizer", "PERFORM", "cf_check_digitizer"}, \
    {"IDS__COR", "PERFORM", "cf_ids_dead_time", "PERFORM", "cf_ids_dead_time"}, \
    {"ELEC_COR", "PERFORM", "cf_electronics_dead_time", "PERFORM", "cf_electronics_dead_time"}, \
    {"FIFO_COR", "OMIT", "cf_fifo_dead_time", "PERFORM", "cf_fifo_dead_time"}, \
    {"DEAD_COR", "PERFORM", "cf_apply_dead_time", "PERFORM", "cf_apply_dead_time"}, \
    {"THRM_COR", "PERFORM", "cf_thermal_distort", "PERFORM", "cf_thermal_distort"}, \
    {"RATE_COR", "PERFORM", "cf_count_rate_y_distort", "PERFORM", "cf_count_rate_y_distort"}, \
    {"TMXY_COR", "PERFORM", "cf_time_xy_distort", "PERFORM", "cf_time_xy_distort"}, \
    {"GEOM_COR", "PERFORM", "cf_geometric_distort", "PERFORM", "cf_geometric_distort"}, \
    {"PHAH_COR", "PERFORM", "cf_modify_hist_pha", "OMIT", "cf_modify_hist_pha"}, \
    {"PHAX_COR", "PERFORM", "cf_pha_x_distort", "PERFORM", "cf_pha_x_distort"}, \
    {"ACTV_COR", "PERFORM", "cf_active_region", "PERFORM", "cf_active_region"}, \
    {"LIMB_COR", "PERFORM", "cf_screen_limb_angle", "PERFORM", "cf_screen_limb_angle"}, \
    {"SAA__COR", "PERFORM", "cf_screen_saa", "PERFORM", "cf_screen_saa"}, \
    {"VOLT_COR", "PERFORM", "cf_screen_high_voltage", "PERFORM", "cf_screen_high_voltage"}, \
    {"BRST_COR", "OMIT", "cf_screen_burst", "PERFORM", "cf_screen_burst"}, \
    {"APER_COR", "PERFORM", "cf_screen_jitter", "PERFORM", "cf_screen_jitter"}, \
    {"UGTI_COR", "OMIT", "cf_set_user_gtis", "PERFORM", "cf_set_user_gtis"}, \
    {"FLAG_COR", "PERFORM", "cf_set_photon_flags", "PERFORM", "cf_set_photon_flags"}, \
    {"GTI__COR", "PERFORM", "cf_set_good_time_intervals", "PERFORM", "cf_set_good_time_intervals"}, \
    {"TIME_COR", "PERFORM", "cf_modify_hist_times", "OMIT", "cf_modify_hist_times"}, \
    {"AIRG_COR", "PERFORM", "cf_screen_airglow", "PERFORM", "cf_screen_airglow"}, \
    {"BPIX_COR", "PERFORM", "cf_screen_bad_pixels", "PERFORM", "cf_screen_bad_pixels"}, \
    {"PHA__COR", "OMIT", "cf_screen_pulse_height", "PERFORM", "cf_screen_pulse_height"}, \
    {"FIND_COR", "PERFORM", "cf_find_spectra", "PERFORM", "cf_find_spectra"}, \
    {"YMOT_COR", "OMIT", "cf_calculate_ycent_motion", "PERFORM", "cf_calculate_ycent_motion"}, \
    {"GRAT_COR", "PERFORM", "cf_grating_motion", "PERFORM", "cf_grating_motion"}, \
    {"FPA__COR", "PERFORM", "cf_fpa_position", "PERFORM", "cf_fpa_position"}, \
    {"MIRR_COR", "PERFORM", "cf_mirror_motion", "PERFORM", "cf_mirror_motion"}, \
    {"JITR_COR", "OMIT", "cf_satellite_jitter", "PERFORM", "cf_satellite_jitter"}, \
    {"YCNT_COR", "PERFORM", "cf_calculate_y_centroid", "PERFORM", "cf_calculate_y_centroid"}, \
    {"CHID_COR", "PERFORM", "cf_identify_channel", "PERFORM", "cf_identify_channel"}, \
    {"TCRT_COR", "PERFORM", "cf_target_count_rate", "PERFORM", "cf_target_count_rate"}, \
    {"ASTG_COR", "PERFORM", "cf_astigmatism", "PERFORM", "cf_astigmatism"}, \
    {"WAVE_COR", "PERFORM", "cf_dispersion", "PERFORM", "cf_dispersion"}, \
    {"DOPP_COR", "PERFORM", "cf_doppler_and_heliocentric", "PERFORM", "cf_doppler_and_heliocentric"}, \
    {"FLAT_COR", "OMIT", "cf_flat_field", "OMIT", "cf_flat_field"}, \
    {"FLUX_COR", "PERFORM", "cf_convert_to_ergs", "PERFORM", "cf_convert_to_ergs"}, \
    {"COMB_COR", "PERFORM", "cf_coadd", "PERFORM", "cf_coadd"} \
}


struct cal_file_tab
{
    char  name[5];
    int   numfiles;
    char  extension[4];
    char  filenames[3][19];
    float aftermjd[3];
    int   interp[3];
};

#define NUMCALKEYS 28

#define CALIBRATION_FILE_KEYS { \
    {"AEFF", 2, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"AIRG", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"ASTG", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"BCHR", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"BKGD", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"CHID", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"DIGI", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"ELEC", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"FLAT", 2, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"GEOM", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"GRAT", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"HSKP", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"JITR", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"MIRR", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"PARM", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"PHAH", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"PHAX", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"QUAL", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"RATE", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"TMXY", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"SAAC", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"SCRN", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"SPEC", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"STIM", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"VOLT", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"WAVE", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"WGTS", 1, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}},\
    {"WORM", 2, "CAL", {"", "", ""}, {0.0, 0.0, LARGEMJD}, {0, 0, 0}} \
}

extern int verbose_level;

void   cf_error_init(const char *, const char *, FILE *);
void   cf_verbose(int, const char *, ...);
void   cf_if_fits_error(int);
void   cf_if_warning(char *, ...);
void   cf_if_error(char *, ...);
void  *cf_malloc(size_t);
void  *cf_calloc(size_t, size_t);

void   cf_timestamp(const char *, const char *, char *);
int    cf_proc_check(fitsfile *, char *);
int    cf_proc_update(fitsfile *, char *, char *);
int    cf_fuv_init(fitsfile *);

char  *cf_cal_file(char *);
char  *cf_parm_file(char *);
int    eclipse(double *, double, double *);
double geod_mag(double, double);
double helio_vel(double, double, double);
double lsrd_vel(double, double);
double lsrk_vel(double, double);
void   month_day(int, int, int*,  int*);
void   read_tle(fitsfile *);
double pole_ang(double *, double *, double, double);
double solar_ang(double, double, double);
int    saa(saareg *, double, double);
double space_vel(double *, double, double);
void   state_geod(double *, double, double *, double *, double *);
double state_limb(double *, double, double, double, double *,int *);
void   cf_velang(fitsfile *, double);
int    astig_check_input_image(fitsfile *);

long   cf_read_col(fitsfile *, int, char *, void **);
int    cf_write_col(fitsfile *, int, char *, void *, long);
int    cf_nint (double);
long   cf_nlong (double);

int    cf_fes_proc_check(fitsfile *, char *);
int    cf_fes_proc_update(fitsfile *, char *, char *);

int    cf_add_header_keywords(fitsfile *);
long   cf_get_times(fitsfile *, double **);
int    cf_get_gti(fitsfile *, double **, double **);
int    cf_get_geocorona(fitsfile *, short **, short **, short **, short **);
int    cf_get_potholes(fitsfile *, float **, float **, float **, float **);
int    cf_timeline(fitsfile *);
int    cf_set_background_limits(fitsfile *);
long   cf_extraction_limits(fitsfile *, int, int, short **, short **,
				short*, short*);

int    cf_check_digitizer(fitsfile *);
int    cf_ids_dead_time(fitsfile *, long, float *, float *);
int    cf_electronics_dead_time(fitsfile *, long, float *, float *);
int    cf_fifo_dead_time(fitsfile *, long, float *, long, float *, float *,
			float *);
int    cf_apply_dead_time(fitsfile *, long, float *, float *,
				long, float *, float *, float *);
int    cf_thermal_distort(fitsfile *, long, float *, float *, float *,
                                unsigned char *);
int    cf_count_rate_y_distort(fitsfile *, long, float *, float *,
			       unsigned char *, long, float *, float *);
int    cf_time_xy_distort(fitsfile *, long, float *, float *, unsigned char *);
int    cf_geometric_distort(fitsfile *, long, float *, float *,
			    unsigned char *);
int    cf_modify_hist_pha(fitsfile *, long, unsigned char *, unsigned char *);
int    cf_pha_x_distort(fitsfile *, long, unsigned char *, float *,
			unsigned char *);
int    cf_active_region(fitsfile *, long, float *, float *, unsigned char *);

int    cf_find_spectra(fitsfile *, long, float *, float *, float *,
		unsigned char *, unsigned char *, unsigned char *, int);
int    cf_identify_channel(fitsfile *, long, float *, float *, unsigned char *,
			unsigned char *, int, int);
int    cf_calculate_ycent_motion(fitsfile*, long, float*, float*,
		unsigned char*, unsigned char*, long, float*, float*, float*);
int    cf_source_aper(fitsfile*, int*);
int    cf_grating_motion(fitsfile *, long, float *, float *, float *, 
			unsigned char *, long, float *, short *);
int    cf_fpa_position(fitsfile *, long, float *, unsigned char *);
int    cf_read_fpa_pos (fitsfile *, float *, float *);
int    cf_mirror_motion(fitsfile *, long, float *, float *, float *, 
			unsigned char *, long, float *, short *);
int    cf_satellite_jitter(fitsfile *, long, float *, float *, float *,
		unsigned char *, long, float *, unsigned char *);
int    cf_calculate_y_centroid(fitsfile*, long, float*, float*, float*,
			unsigned char*, unsigned char*, unsigned char*);
int    cf_target_count_rate(fitsfile *, long, float *, float *, unsigned char *,
			unsigned char *, long, float *, float *, float *);
int    cf_screen_limb_angle(fitsfile *, long, unsigned char *, float *);
int    cf_screen_saa(fitsfile *, long, unsigned char *, float *, float *);
int    cf_screen_airglow(fitsfile *, long, float *, float *, unsigned char *);
int    cf_screen_bad_pixels(fitsfile *, long, float *, float *, unsigned char *);
int    cf_screen_burst(fitsfile *, long, float *, float *, float *,
		unsigned char *, GTI *, long, float *, unsigned char *, float *,
		short *);
int    cf_screen_jitter(fitsfile *, long, float *, unsigned char *);
int    cf_set_user_gtis(fitsfile *, long, float *, unsigned char *);
int    cf_screen_high_voltage(fitsfile *, long, unsigned char *, short *);
int    cf_screen_pulse_height(fitsfile *, long, unsigned char *, unsigned char *);
int    cf_set_photon_flags(fitsfile *, long, float *, float *, unsigned char *,
			unsigned char *, long, float *, unsigned char *);
int    cf_set_good_time_intervals(fitsfile *, long, float *, unsigned char *,
			GTI *);
int    cf_modify_hist_times(fitsfile *, long, float *, GTI *);
int    cf_astigmatism(fitsfile *, long, float *, float *, unsigned char *);
int    cf_dispersion(fitsfile *, long , float *, unsigned char *, float *);
int    cf_doppler_and_heliocentric(fitsfile *, long, float *, unsigned char *,
		float *, long, float *, float *);
int    cf_convert_to_ergs(fitsfile *, long , float *, float *,
			unsigned char *, float *);
int    cf_apply_filters(fitsfile *, int, long, unsigned char *, unsigned char *,
			long, unsigned char *, long *, long *,	long *, long **);
int    cf_scale_bkgd(fitsfile *, long, float *, float *, float *, unsigned char *, 
        unsigned char *, unsigned char *, long,
	long *, long, long, int *, int *, int *, int *, int *, float **,
	int *, int *, int *, float **);
int    cf_make_mask(fitsfile *, int, long, float *, int, int, float **);
int    cf_make_wave_array(fitsfile *, int, long *, float **);
int    cf_rebin_probability_array(fitsfile *, int, int, long, float *, int *,
			float *, float **);
int    cf_rebin_background(fitsfile *, int, long, float *,
		int, int, int, int, float *, float **);
int    cf_standard_or_optimal_extraction(fitsfile *, int *);
int    cf_optimal_extraction(fitsfile *, int, int, float *, float *,
		unsigned char *, float *, long, long *, float *,
		float *, int, float, float *, long, float *, float **,
		float **, long **, float **, float **, short **);
int    cf_write_extracted_spectrum(fitsfile *, int, int, long, float *,
		float *, float *, long *, float *, float *, short *, char *);
char   *cf_hist_file(char *);
