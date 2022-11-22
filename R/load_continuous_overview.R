source(file.path(IE_root, 'R', 'continuous_IE_detection_init.R'))

reg_method = 'glm'
reg_method = 'lmp'
reg_method = 'glm_log'
ncores = 40
ncores = 36
ncores = 20
ncores = 1
ncores = floor(12 * 12/30 * 12/30)
ncores = floor(12 * 12/30)
ncores = 24

## We filter based on the p-value of beta_y, which can fairly be
## assumed to be significant.
setting_dtf <-
  compile_all_coef_overview(
    redo = F,
    redo_subanalyses = F,
    check_data_availability = F,
    skip_missing = F,
    # N_downsample_settings = 50,
    ncores = ncores,
    reg_method = reg_method,
    include_non_ds_tallies = T,
    z_normalize = F
  )

# setting_dtf$patient_inclusion_crit
# duplicated(colnames(setting_dtf))
# levels(setting_dtf$LOH_HLA)

f_setting_dtf <-
  setting_dtf %>%
  filtering(
    args_list = list(
      min_nz_patients = 100,
      min_project_size = 250,
      intercept_filter_magnitude = 1e-3,
      # max_delta_n_cov = 1,
      # force_positive_intercept = F,
      force_intercept_ge_rc = T,
      # intercept_filter_magnitude = 1e-2,
      intercept_filter_p_val = 1e-3
    )
  )

p_setting_dtf <- pretty_overview_dtf(f_setting_dtf)

test_opt_string_consistency(setting_dtf)
test_opt_string_consistency(f_setting_dtf)
test_opt_string_consistency(p_setting_dtf)
