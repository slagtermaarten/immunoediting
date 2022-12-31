# getDTthreads()
setDTthreads(1)

IE_requires_computation <- mod_time_comparator(
  # minimum_mod_time = '2021-4-21 22:36', verbose = TRUE)
  # minimum_mod_time = '2021-4-23 13:51', verbose = TRUE)
  # minimum_mod_time = '2021-4-24 11:51', verbose = TRUE)
  # minimum_mod_time = '2022-8-25 11:51', verbose = TRUE)
  # minimum_mod_time = '2022-9-13 11:51', verbose = TRUE)
  # minimum_mod_time = '2022-9-14 11:51', verbose = TRUE)
  # minimum_mod_time = '2022-9-21 14:31', verbose = TRUE)
  # minimum_mod_time = '2022-10-03 14:31', verbose = TRUE)
  # minimum_mod_time = '2022-10-23 14:31', verbose = TRUE)
  # minimum_mod_time = '2022-11-03 11:00', verbose = TRUE)
  minimum_mod_time = '2022-11-19 11:00', verbose = TRUE)
# 2022-10-04 16:00 Implemented ACF stuff, but didn't bother to reset
# timer while only a subset of analyses has run last night as I
# probably won't need the ACF stuff and will focus on permutatation
# p-values instead, so about 25% of subanalyses will have ACF now.

N_perms_global = 250

## {{{ Constants
scalar_analysis_components <- c('test_yr', 'p_val', 'p_val_no_reg',
                                'logFC', 'n_ref', 'n_test')
other_analysis_components <- c('reg_plot', 'lm', 'analysis_name',
                               'reg_plot_by_project')
analysis_components <- c(scalar_analysis_components,
  other_analysis_components)
analysis_grp_vars <- c('overlap_var', 'patient_inclusion_crit',
  'LOH_HLA', 'analysis_name', 'analysis_idx', 'project_extended')

tumor_types <- c(
  `Pan*'-'*cancer` = "pan_cancer",
  `Adrenal~gland` = "adrenal_gland",
  `Bile~duct` = "bile_duct",
  Bladder = "bladder",
  `Breast~Basal` = "breast_basal",
  `Breast~Her2` = "breast_her2",
  `Breast~LumA` = "breast_luma",
  `Breast~LumB` = "breast_lumb",
  `Breast~Normal^{'like'}` = "breast_normallike",
  Cervix = "cervix",
  `Cervix~HPV^{'-'}` = "cervix_hpv_min",
  `Cervix~HPV^{'+'}` = "cervix_hpv_plus",
  `Colon~MSI^{'H'}` = "colon_msih",
  `Colon~non*'-'*MSI^{'H'}` = "colon_non_msih",
  DLBCL = "dlbcl",
  `DLBCL~EBV^{'-'}` = "dlbcl_ebv_min",
  `DLBCL~EBV^{'+'}` = "dlbcl_ebv",
  Esophagus = "esophagus",
  `Follicular~lymphoma~EBV^{'-'}` = "follicular_lymphoma_ebv_min",
  Glioblastoma = "glioblastoma",
  `Head~and~neck` = "head_and_neck",
  `Head~and~neck~HPV^{'-'}` = "head_and_neck_hpv_min",
  `Head~and~neck~HPV^{'+'}` = "head_and_neck_hpv_plus",
  `Kidney~chromophobe` = "kidney_chromophobe",
  `Kidney~clear~cell` = "kidney_clear_cell",
  `Kidney~papillary` = "kidney_papillary",
  Liver = "liver",
  `Lower~grade~glioma` = "lower_grade_glioma",
  `Lung~adeno` = "lung_adeno",
  `Lung~squamous` = "lung_squamous",
  `Lymphoid~leukemia` = "lymphoid_leukemia",
  Melanoma = "melanoma",
  `Cutaneous melanoma` = "melanoma",
  `Multiple~myeloma` = "multiple_myeloma",
  `Myeloid~leukemia` = "myeloid_leukemia",
  `NA` = "na",
  Ovary = "ovary",
  Pancreas = "pancreas",
  Paraganglion = "paraganglion",
  Prostate = "prostate",
  Rectum = "rectum",
  Sarcoma = "sarcoma",
  `Stomach~EBV^{'-'}` = "stomach_ebv_min",
  `Stomach~EBV^{'+'}` = "stomach_ebv_plus",
  `Stomach~EBV^{'-'}~MSI^{'H'}` = "stomach_ebv_min_msih",
  Testes = "testes",
  Thyroid = "thyroid",
  `Uterine~endometrial` = "uterine_endometrial",
  Uterus = "uterus",
  `Uveal~melanoma` = "uveal_melanoma")

tumor_types_inv <- setNames(names(tumor_types), tumor_types)
### }}}


determine_p_var <- function(fill_var) {
  if (fill_var == 'rc_CYT') {
    return('p_val_CYT')
  } else if (grepl('bayesian|estimate', fill_var)) {
    return('post_prob')
  } else {
    return('p_val')
  }
}


#' Determine threshold between estimates indicating neo-antigen
#' depletion/enrichment with increased PS
#'
#'
determine_threshold <- function(fill_var) {
  p_var <- switch(fill_var,
    'rc' = 0,
    'rc_CYT' = 0,
    'yr_fractional_change' = 0,
    'effective_editing_max' = 1,
    'depletion_max' = 1,
    'depletion_full' = 1,
    'depletion_mean' = 1)
}


#' Generate file name(s) to be generated in \code{test_continuous_IE}
#' and other functions
#'
#' @param analysis_name One of: twoD_sens_analysis,
#' marty_param_titration, driv_ess_param_titration,
#' clon_param_titration, rooney_param_titration
gen_cont_IE_fn <- function(
  base_name = 'cont-IE',
  focus_allele = 'A0201',
  root_dir = rds_dir,
  sub_dir = '',
  analysis_name = 'twoD_sens_analysis',
  project_extended = NULL,
  reg_method = 'rlm',
  LOH_HLA = 'no_LOHHLA',
  hla_sim_range = NULL,
  analysis_idx = '',
  replace_idx = F,
  fill_var = 'rc',
  p_val_bound = NULL,
  overlap_var = 'mean_score',
  permute_input = as.character(c()),
  extension = 'rds',
  z_normalize = F,
  patient_inclusion_crit = '') {

  if (!exists('base_name') || missing(base_name) ||
      is.null(base_name) || is.na(base_name)) {
    base_name <- ''
  }

  if (replace_idx) {
    neo_settings <- as.list(ds_param_grid[as.integer(analysis_idx), ])
    analysis_idx <-
      complement_options(neo_settings) %>%
      gen_opts_string(no_spaces = T, no_hla = T)
  }

  analysis_idx <- prepend_hyphen(analysis_idx)
  analysis_name <- prepend_hyphen(analysis_name)
  focus_allele <- prepend_hyphen(focus_allele)
  hla_sim_range <- ifelse(
    is.null(hla_sim_range) || hla_sim_range == '',
    '', paste(hla_sim_range, collapse = '_')) %>%
    prepend_hyphen
  reg_method <- prepend_hyphen(reg_method)
  LOH_HLA <- prepend_hyphen(LOH_HLA)
  patient_inclusion_crit <- prepend_hyphen(patient_inclusion_crit)
  overlap_var <- prepend_hyphen(overlap_var)
  p_val_bound <- ifelse(eps(p_val_bound, 0.05), '',
    sprintf('-pvb-%s', format(p_val_bound,
        scientific = T, digits = 2)))
  fill_var <- ifelse(fill_var == 'rc', '', prepend_hyphen(fill_var))
  tumor_type <- ifelse(is.null(project_extended), '',
    tumor_types[project_extended])
  z_normalize <- ifelse(z_normalize == T, '', '-not_z_normalized')
  permute_input <- ifelse(length(permute_input) == 0, '',
    paste0('-permute_input=', paste0(permute_input, collapse = '-')))
  extension <- prepend_string(extension, prepend_string = '.')

  file_head <-
    tryCatch(paste0(
        base_name, tumor_type, analysis_name, focus_allele,
        overlap_var, LOH_HLA, patient_inclusion_crit, reg_method,
        p_val_bound, fill_var, hla_sim_range, z_normalize,
        analysis_idx, permute_input, extension),
      error = function(e) { print(e); browser() })
  ## Strip preceeding hyphen if it's there
  file_head <- gsub('^-', '', file_head)
  if (!is.null(root_dir) && !is.na(root_dir) && root_dir != '') {
    if (!is.null(sub_dir) && !is.na(sub_dir) && root_dir != '') {
      root_dir <- file.path(root_dir, sub_dir)
    }
    fn <- file.path(root_dir, file_head)
    dir.create(root_dir, showWarnings = F, recursive = T)
  } else {
    fn <- file_head
  }

  return(fn)
}


#' Select rows corresponding to a particular tumor type or return the
#' entire table if queried with 'Pan-cancer'
#'
#'
subset_project <- function(dtf, project_extended = NULL,
  pan_can_includes_all = TRUE) {

  if (is.null(project_extended) || is.na(project_extended)) {
    return(dtf)
  }

  pancan_mode <- project_extended %in%
    c('Pan*\'-\'*cancer', 'Combined', 'Pan-cancer', 'pan_cancer',
      'pancancer')

  if (!pancan_mode || (pancan_mode && !pan_can_includes_all)) {
    l_project_extended <- as.character(project_extended)
    setDT(dtf)
    data_subs <- dtf[project_extended == l_project_extended]
    browser(expr = length(unique(data_subs$project_extended)) > 1)
  } else {
    data_subs <- dtf
  }
  return(data_subs)
}


#' Apply Z-scaling to the y_var column, storing pre-normalization
#' stats
#'
#'
normalize_columns <- function(dtf) {
  if (null_dat(dtf)) return(NULL)
  setDT(dtf)
  varnames <- c('ol', 'CYT', 'y_var')
  varnames <- c('y_var')
  for (varn in varnames) {
    res <- scale(dplyr::pull(dtf, varn))
    dtf[, (varn) := as.vector(res)]
    dtf[, (sprintf('%s_sd', varn)) := attr(res, 'scaled:scale')]
    dtf[, (sprintf('%s_mean', varn)) := attr(res, 'scaled:center')]
  }
  return(dtf)
}


#' Reduce donor_summary object to columns that are essential to the
#' continuous IE analyses in order to save RAM while doing the
#' analyses
#'
#'
reduce_to_essential_IE_columns <- function(
  dtf, extra_columns = NULL) {
  if (null_dat(dtf)) return(NULL)
  setDT(dtf)
  candidate_cols <- c('project_extended', 'donor_id', 'ol',
    'CYT', 'c', 'i', 'mean_score', 'mean_score_AB', 'mean_score_C',
    'weight', 'IE_essentiality_impaired', 'CYT', 'hla_allele_status',
    'hla_allele_status_b', 'y_var', 'y_var_sd') %>%
    c(extra_columns)
  dtf <- tryCatch(
    dtf[, intersect(colnames(dtf), candidate_cols), with = F],
    error = function(e) { print(dtf); browser() })
  return(dtf)
}


#' Determine the dependent variable in any \code{donor_summary} like
#' object
#'
#' The first variable starting with a 'c' is generally taken as the
#' variable defining the amount of neo-antigens created
#'
determine_continuous_IE_yval <- function(dtf, replace_zeroes = T) {
  ## c_name will 'c_muts.missense_mutation' for regular
  ## donor_summaries
  c_name <- grep('^c.*', names(dtf), value = T) %>%
    { grep('silent_mutation', ., invert = T, value = T) } %>%
    .[1]
  i_name <- gsub('c_(.*)', 'i_\\1', c_name)
  if (!is.na(c_name) && c_name != 'NA') {
    num <- dplyr::pull(dtf, c_name)
    denum <- dplyr::pull(dtf, i_name)
    dtf$y_var <- num / denum
    if (replace_zeroes) {
      dtf$y_var[which(is.na(num) & !is.na(denum))] <- 0
    } else {
      all(dtf$y_var > 0, na.rm = T)
    }
    dtf$c <- ifelse(is.na(num), 0, num)
    dtf$i <- denum
    dtf$weight <- dplyr::pull(dtf, i_name)
    attr(dtf, 'y_label') <- 'Neo-antigen yield rate'
  } else {
    if ('B' %in% colnames(dtf)) {
      dtf$y_var <- dtf$B
      dtf$weight <- 1
      attr(dtf, 'y_label') <- 'Predicted / expected neo-antigens'
    } else {
      perm_browser()
    }
  }
  return(dtf)
}
if (F) {
  with(new.env(), {
    ds <- get_donor_summary()
    ds_f <- determine_continuous_IE_yval(ds)
    all(ds$c_muts.missense_mutation /
      ds$i_muts.missense_mutation == ds_f$y_var,
        na.rm = T
    ) %>% stopifnot()
    # all(1 == ds_f$y_var, na.rm = T) %>% stopifnot()
  })
}


make_names_df_friendly <- function(v) {
  tolower(gsub('\\.|\\s', '_', v))
}


ds_ol_stats <- function(dtf) {
  if (null_dat(dtf)) return(NULL)
  dtf$ol_b <- cut(dtf$ol, seq(0, 1, by = .1), labels = F)
  # h_dist <- table(dtf$ol_b)
  # h_dist <- set_names(h_dist, paste('d_', names(h_dist), sep = ''))
  h_stats <- dtf %>%
    group_by(ol_b) %>%
    dplyr::mutate(
      w_i = sum(i, na.rm = T),
      w_c = sum(c, na.rm = T),
      w_y_var = w_c / w_i
    ) %>%
    summarize(n = n(), across(c(w_y_var, w_i, w_c, y_var, i, c),
      c('median' = median, 'mad' = mad))) %>%
    dplyr::select(!matches('w_.*mad'))
  return(h_stats)
}


test_data_sufficiency <- function(dtf,
  min_pts = 20,
  min_points_lq = 5,
  min_points_uq = 5) {

  if (null_dat(dtf)) return(F)
  setDT(dtf)
  discretized_ps <- dtf[, .N, cut(ol, breaks = seq(0, 1, by = .25))]
  discretized_ps <- discretized_ps[!is.na(cut)]

  lq <- '(0,0.25]'
  uq <- '(0.75,1]'

  ## st := simple test
  st <- function(v) !is.null(v) && length(v) > 0 && v
  out <-
    !is.null(dtf) &&
    nrow(dtf) >= min_pts &&
    st(with(discretized_ps, N[cut == lq] > min_points_lq)) &&
    st(with(discretized_ps, N[cut == uq] > min_points_uq)) &&
    # any(dtf$y_var > 0) &&
    T

  if (is.na(out)) browser()
  return(out)
}


compute_dtf_stats <- function(dtf) {
  out <- list(
    'n_patients' = dtf[, .N],
    'n_nz_patients' = dtf[c > 0, .N],
    'mean_yr' = mean(dtf[, c/i])
  )
  return(out)
}


filter_low_TMB_patients <- function(dtf) {
  ## Subselect patients that have a fighting chance of getting at
  ## least one neo-antigen with the current filtering settings. Use
  ## a simple model to estimate that TMB threshold
  if (null_dat(dtf)) return(NULL)
  dtf$i <- tryCatch(round(dtf$i), error = function(e) { browser() })
  dtf$c <- round(dtf$c)
  # dtf$c_q <- dtf$c / max(dtf$c)
  # dtf$i_q <- dtf$i / max(dtf$i)
  nz_mod <- fit_nb_mod(dtf, c ~ 1 + i)
  # test_plot({ par(mfrow = c(2, 2)); plot(nz_mod) }, w = 17.4, h = 20)

  ## Test how many variants are needed for 1 (0 on the log scale of
  ## the coefficients) expected mutation
  min_req_muts <- coef(nz_mod) %>%
    { -1 * .['(Intercept)'] / .['i'] } %>%
    unname
  l_dtf <- dtf[i > min_req_muts, ]

  nz_report <- list(
    nz_converged = nz_mod$converged,
    nz_intercept = coef(nz_mod)[1],
    nz_rc = coef(nz_mod)[2],
    nz_theta = nz_mod$theta
  )
  attr(dtf, 'nz_report') <- nz_report
  return(dtf)
}


replace_values <- function(vec, reps) {
  out <- vec
  for (i in seq_along(reps)) {
    idxs <- setdiff(which(vec == names(reps)[i]), NA)
    out <- replace(out, idxs, unname(reps[i]))
  }
  return(out)
}
# replace_values(c('a', 'b', 'c', 'c'), c('c' = 'd', 'e' = 'f', 'g' = 'a'))
# replace_values(c('f', 'b', 'c', 'c'), c('c' = 'd', 'e' = 'f', 'g' = 'a'))


broom2row <- function(broom_out,
  term_rename = c('(Intercept)' = 'offset', 'i' = 'intercept',
    'i:ol' = 'rc', 'i:CYT' = 'CYT_rc', 'i:CYT:ol' = 'CYT_ps_rc')) {
  broom_out$term <- replace_values(broom_out$term, term_rename)

  out <- tidyr::pivot_wider(broom_out,
    values_from = !term, names_from = term)
  ## Make naming consistent with earlier work
  colnames(out) <-
    stringr::str_replace_all(colnames(out), '\\.', '_')
  colnames(out) <-
    stringr::str_replace_all(colnames(out), 'estimate_', '')
  colnames(out) <-
    stringr::str_replace_all(colnames(out), 'p_value', 'p_val')
  colnames(out) <- replace_values(colnames(out),
    c('p_val_rc' = 'p_val', 'df' = ''))

  return(out)
}


apply_intercept_scale <- function(
  dtf,
  intercept_scale = NULL,
  find_intercept_method = 'rlm') {

  if (!is.null(intercept_scale)) {
    ## Adjust the data to be more in the expected range
    if (find_intercept_method == 'lm') {
      orig_b0 <- tryCatch(coef(lm(y_var ~ ol, dtf))[1],
        error = function(e) { NA })
    } else if (find_intercept_method == 'rlm') {
      orig_b0 <- tryCatch(fit_rlm_model(dtf)$intercept,
        error = function(e) { NA })
    }
    if (is.na(orig_b0)) return(NULL)
    dtf$y_var <- dtf$y_var * (intercept_scale / orig_b0)
    attr(dtf, 'orig_b0') <- orig_b0
  }
  return(dtf)
}


get_repertoire_overlap <- function(ro = NULL, hla_alleles = 'A0201',
                                   ncores = 36) {
  if (is.null(ro)) {
    ro <- compute_PS_overview(hla_alleles = hla_alleles,
                              mail_notify = T,
                              verbose = F,
                              ncores = ncores)
  } else {
    ## ro <- ro
  }
  return(ro)
}


recover_tumor_type <- function(
  tumor_type = NULL,
  project_extended = NULL) {
  if (!is.null(tumor_type) && tumor_type %in% tumor_types) {
    return(tumor_type)
  } else if (is.null(tumor_type) && !is.null(project_extended)) {
    return(tumor_types[as.character(project_extended)])
  }
}
# recover_tumor_type(project_extended = 'Rectum')


test_opt_string_consistency <- function(dtf) {
  setDT(dtf)
  if (!'opt_string' %in% colnames(dtf)) {
    # message('Skipping STS-opt string consistency check')
    return(NULL)
  }
  if (any(c('TRUE', 'FALSE') %in% dtf$sts_filtering)) {
    stopifnot(dtf[sts_filtering == TRUE, 
      all(grepl('STS', opt_string))]) 
    stopifnot(dtf[sts_filtering == FALSE, 
      !all(grepl('STS', opt_string))])
  } else {
    stopifnot(dtf[sts_filtering == 'STS',
      all(grepl('STS', opt_string))]) 
    stopifnot(dtf[sts_filtering == 'no STS', 
      !all(grepl('STS', opt_string))])
  }
  return(TRUE)
}


prettify_focus_allele <- function(dtf) {
  setDT(dtf)
  if ('focus_allele' %in% colnames(dtf)) {
    if (!is.factor(dtf$focus_allele)) {
      dtf$focus_allele <- factor(dtf$focus_allele)
    }
    if (!any(grepl('\\*', levels(dtf$focus_allele)))) {
      new_levs <- levels(dtf$focus_allele) %>%
          { setNames(., ppHLA(.)) }
      dtf[, focus_allele :=
        forcats::fct_recode(focus_allele, !!!new_levs)]
      test_opt_string_consistency(dtf)
      # levels(dtf$focus_allele)
    }
  }
  return(dtf)
}


prettify_APR <- function(dtf) {
  setDT(dtf)
  if ('percentile_rank' %in% colnames(dtf)) {
    if (!is.factor(dtf$percentile_rank)) {
      dtf$percentile_rank <- factor(dtf$percentile_rank)
    }
    if (!any(grepl('APR', levels(dtf$percentile_rank)))) {
      perc_rank_print <- purrr::map_chr(
        auto_name(levels(dtf$percentile_rank)),
        ~glue::glue('APR={.x}')) %>%
        { setNames(names(.), .) }
      dtf[, percentile_rank :=
        forcats::fct_recode(percentile_rank, !!!perc_rank_print)]
    }
  }
  test_opt_string_consistency(dtf)
  return(dtf)
}


prettify_STS <- function(dtf) {
  setDT(dtf)
  if ('sts_filtering' %in% colnames(dtf)) {
    if (!is.factor(dtf$sts_filtering)) {
      dtf$sts_filtering <- factor(dtf$sts_filtering)
    }
    if (!any(grepl('STS', levels(dtf$sts_filtering)))) {
      new_levs <- c('STS' = 'TRUE', 'no STS' = 'FALSE')
      dtf[, sts_filtering_n :=
        forcats::fct_recode(sts_filtering, !!!new_levs)]
      dtf[, .N, keyby = .(sts_filtering, sts_filtering_n)]
      dtf[, sts_filtering := sts_filtering_n]
      dtf[, sts_filtering_n := NULL]
    }
  }
  test_opt_string_consistency(dtf)
  return(dtf)
}


prettify_LOHHLA <- function(dtf) {
  setDT(dtf)
  if ('LOH_HLA' %in% colnames(dtf)) {
    if (!is.factor(dtf$LOH_HLA)) {
      dtf$LOH_HLA <- factor(dtf$LOH_HLA)
    }
    if (any(grepl('no_', levels(dtf$LOH_HLA)))) {
      new_levs <- c(
        'no LOH HLA' = 'no_LOHHLA', 
        'LOHHLA' = 'LOH HLA',
        'strict LOH HLA' = 'strict_LOHHLA'
      )
      dtf[, LOH_HLA :=
        forcats::fct_recode(LOH_HLA, !!!new_levs)]
    }
  }
  test_opt_string_consistency(dtf)
  return(dtf)
}


prettify_overlap_var <- function(dtf) {
  setDT(dtf)
  if ('overlap_var' %in% colnames(dtf)) {
    if (!is.factor(dtf$overlap_var)) {
      dtf$overlap_var <- factor(dtf$overlap_var)
    }
    if (any(grepl('_', levels(dtf$overlap_var)))) {
      new_levs <- c(
        'HLA A, B, C' = 'mean_score', 
        'HLA A, B' = 'mean_score_AB'
      )
      dtf[, overlap_var :=
        forcats::fct_recode(overlap_var, !!!new_levs)]
    }
  }
  test_opt_string_consistency(dtf)
  return(dtf)
}


prettify_variant_selection <- function(dtf) {
  setDT(dtf)

  # levels(f_setting_dtf$analysis_name)

  if ('analysis_name' %in% colnames(dtf)) {
    if (!is.factor(dtf$analysis_name)) {
      dtf$analysis_name <- factor(dtf$analysis_name)
    }
    if (any(grepl('_', levels(dtf$analysis_name)))) {
      new_levs <- 
        with(display_settings$analysis_name, setNames(breaks, labels))
      dtf[, analysis_name :=
        forcats::fct_recode(analysis_name, !!!new_levs)]
      # levels(dtf$analysis_name)
    }
  }
  test_opt_string_consistency(dtf)

  return(dtf)
}


format_overview_res <- function(
  dtf,
  reg_method = NULL,
  z_normalize = NULL,
  ds_frac = NULL,
  iter = NULL) {

  setDT(dtf)

  if (all(c('intercept', 'rc') %in% colnames(dtf))) {
    dtf[, yr_fractional_change := rc / intercept]
  } else if ('estimate' %in% colnames(dtf)) {
    dtf[, yr_fractional_change := estimate]
  }

  if ('hla_sim_range' %in% colnames(dtf) &&
      is.factor(dtf$hla_sim_range)) {
    dtf[, hla_sim_range :=
    c('all', '<= 1')[as.integer(sapply(hla_sim_range, is.null)) + 1]]
    dtf[, hla_sim_range := as.factor(hla_sim_range)]
  }

  if ('analysis_name' %in% colnames(dtf)) {
    dtf[, analysis_name :=
      factor(analysis_name, levels = analysis_names)]
  }

  if (!is.null(reg_method)) {
    dtf[, 'reg_method' := reg_method]
  }

  if (!is.null(z_normalize)) {
    dtf[, 'z_normalize' := z_normalize]
  }

  if (!is.null(ds_frac)) {
    dtf[, 'ds_frac' := ds_frac]
  }

  if (!is.null(iter)) {
    dtf[, 'iter' := iter]
  }

  if (all(c('scale', 'intercept') %in% colnames(dtf))) {
    if (all(dtf$intercept == 1)) {
      dtf[, 'norm_scale' := scale]
    } else {
      dtf[, 'norm_scale' := scale / intercept]
    }
  }

  if (all(c('perm_delta_mean_c_mean',
        'perm_delta_mean_c_ci_l') %in%
      colnames(dtf)) &&
      !'perm_delta_mean_c_SE' %in% colnames(dtf)) {
    dtf[, 'perm_delta_mean_c_SE' :=
      (perm_delta_mean_c_mean - perm_delta_mean_c_ci_l) / 1.96]
  }

  if ('baseline_yr' %in% colnames(dtf)) {
    dtf[, baseline_yr := NULL]
  }

  if ('offset' %in% colnames(dtf)) {
    ## offset is \beta_0
    ## intercept is \beta_r
    ## rc is \beta_h
    dtf[, delta := rc / (offset / median_TMB + intercept)]
  }

  if ('perm_delta_pq' %in% colnames(dtf)) {
    dtf <- dtf[!is.na(perm_delta_pq)]
  }

  if (T) {
    if (all(c('delta', 'perm_delta_median') %in% colnames(dtf))) {
      dtf[, 'delta_n' := delta - perm_delta_median]
    }
  } else {
    if (all(c('delta', 'perm_delta_mean') %in% colnames(dtf))) {
      dtf[, 'delta_n' := delta - perm_delta_mean]
    }
  }

  stopifnot('analysis_idx' %in% colnames(dtf))
  dtf <- cbind(dtf, ds_param_grid[dtf[, analysis_idx], ])
  dtf$VE_threshold %<>% friendly_factor
  dtf$expression_threshold %<>% friendly_factor
  dtf$sts_filtering %<>% friendly_factor
  dtf$percentile_rank %<>% friendly_factor
  dtf$processing_threshold %<>% friendly_factor

  if ('delta_n' %in% colnames(dtf)) {
    l_id_vars <- setdiff(id_vars, c('tumor_type', 'focus_allele'))
    dtf[, 'delta_n_cov' := sd(delta_n) / abs(mean(delta_n)),
      by = l_id_vars]
  }

  dtf <- format_coef_overview_(dtf)
  dtf <- recode_RNA_expression(dtf)

  return(dtf)
}


pretty_overview_dtf <- function(dtf) {
  dtf <- prettify_focus_allele(dtf)
  dtf <- prettify_APR(dtf)
  dtf <- prettify_STS(dtf)
  dtf <- prettify_variant_selection(dtf)
  dtf <- prettify_LOHHLA(dtf)
  dtf <- prettify_overlap_var(dtf)
  return(dtf)
}


variabilize_project_extended <- function(tumor_type) {
  # tumor_type %>%
  #   tolower %>%
  #   { gsub('\'|\\*|\\+|\\{|\\}|\\^', '', .) } %>%
  #   { gsub('-', '_', .) } %>%
  #   { gsub('\\~', '_', .) } %>%
  #   { gsub('_{2,}', '_', .) } %>%
  #   { gsub('_$', '', .) }
  tumor_types[tumor_type]
}


recover_project_extended <- function(pe) {
    c("Pan*'-'*cancer", donor_summary[, levels(project_extended)]) %>%
    auto_name %>%
    sapply(variabilize_project_extended) %>%
    { setNames(names(.), .) } %>%
    { .[pe] }
}
# recover_project_extended('stomach_ebv')
# c("Pan*'-'*cancer", donor_summary[, levels(project_extended)]) %>%
#   sapply(variabilize_project_extended)


cat_padded_msg <- function(msg, width = 60) {
  cat(paste0(
    '\n', paste(rep('=', width), collapse = ''), '\n',
    stringr::str_pad(paste('', msg, ''),
      width = width, side = 'both', pad = '='), '\n',
    paste(rep('=', width), collapse = ''), '\n',
    collapse = '')
  )
}


print_overview_stats <- function(dtf, 
  plot_fishtails = NULL, stage_id = '') {
  bayesian_analysis <-
    all(c('evid_ratio', 'intercept_estimate') %in% colnames(dtf))

  latest_stats <- attr(dtf, 'latest_stats') %||% NULL

  if (!is.null(stage_id) && stage_id != '') {
    cat_padded_msg(msg = stage_id)
  }

  if ('focus_allele' %in% colnames(dtf)) {
    cat('\nAllele frequency stats\n')
    print(dtf[, table(focus_allele)] %>% { . / sum(.) })
  }

  cut_breaks <- c(-Inf, -1e-16, 1e-16, Inf)
  if ('intercept' %in% colnames(dtf)) {
    cat('\nlm frequency stats\n')
    setDT(dtf)
    # dtf[, 'exp_result' := intercept > 0 & rc < 0]
    dtf[, .N, keyby = .(
      'intercept' = cut(intercept, breaks = cut_breaks),
      'rc' = cut(rc, breaks = cut_breaks)
    )] %>%
    .[, 'frac' := N / sum(N)] %>% print
  }

  if ('scale' %in% colnames(dtf)) {
    cat('\nSummary of scale\n')
    print(dtf[, summary(scale)])
  }

  if ('n_patients' %in% colnames(dtf)) {
    cat('\nSummary of N patients\n')
    print(dtf[, summary(n_patients)])
  }

  ## This assumes that all levels are observed
  N_possible_analyses <- dtf %>%
    dplyr::select(focus_allele, overlap_var, patient_inclusion_crit,
                  project_extended, LOH_HLA, analysis_name,
                  analysis_idx) %>%
    as.list %>%
    map(~unique(.x)) %>%
    { do.call(expand.grid, .) } %>%
    nrow

  N_unique_observed_wo_allele <- dtf %>%
    dplyr::distinct(overlap_var, patient_inclusion_crit,
                    project_extended,
                    LOH_HLA, analysis_name, analysis_idx) %>%
    nrow

  cat(sprintf('%d/%d (%s) analyses identified', nrow(dtf),
      N_possible_analyses,
      scales::percent(nrow(dtf) / N_possible_analyses)), '\n')

  cat(sprintf('%d/%d (%s) analyses identified (grouping HLA-alleles)',
      N_unique_observed_wo_allele,
      N_possible_analyses / length(focus_hlas),
      scales::percent(N_unique_observed_wo_allele /
        (N_possible_analyses / length(focus_hlas)) )
      ), '\n')

  if (F && 'yr_fractional_change' %in% colnames(dtf)) {
    cat('\nSummary of yr_fractional_change\n')
    print(dtf[, summary(yr_fractional_change)])
  }

  if ('perm_delta_mean_pc' %in% colnames(dtf)) {
    cat('\nSummary of perm_delta_mean_pc\n')
    print(dtf[, summary(perm_delta_mean_pc)])
  }

  if (F && 'rc' %in% colnames(dtf)) {
    cat('\nSummary of rc\n')
    print(dtf[, summary(rc)])
  }

  if ('log2_intercept_cov' %in% colnames(dtf)) {
    cat('\nSummary of -log2_intercept_cov\n')
    print(dtf[, summary(-log2_intercept_cov)])
  }

  if ('post_prob' %in% colnames(dtf)) {
    cat('\nSummary of post_prob\n')
    dtf[, .('N'=.N, 'expected'=nrow(dtf)/20),
        keyby=cut(post_prob, breaks = 20)] %>%
      dplyr::mutate(fd = log2(N) - log2(expected)) %>%
      print
  }

  if ('delta_CI_h' %in% colnames(dtf)) {
    cat('\nSummary of delta_CI_l\n')
    print(dtf[, summary(delta_CI_l)])
    cat('\nSummary of delta_CI_h\n')
    print(dtf[, summary(delta_CI_h)])
    cat('\nSummary of delta_SE\n')
    print(dtf[, summary(delta_SE)])
    cat('\n', '# CI_h < 0: ', dtf[delta_CI_h < 0, .N])
    cat('\n', '# CI_l > 0: ', dtf[delta_CI_l > 0, .N])
  }

  if (F && !is.null(latest_stats)) {
    print(latest_stats, n = 1e9)
  }

  # dtf[, mean(exp_result), list(analysis_idx, focus_allele)] %>%
  #   dplyr::pull(V1) %>%
  #   hist(breaks = 100)
  # library(fgsea)

  # # ordering <- dtf[, .SD, by = analysis_grp_vars]
  # ordering <- dtf[order(, .SD, by = analysis_grp_vars]
  # ordering <- dtf[, .('fraction_exp' = mean(exp_result)),
  #         list(analysis_idx)] %>%
  #   .[order(fraction_exp)]
  # # dtf[, mean(exp_result), by = list(get(analysis_grp_vars[6]))]
  # mod <- glm(exp_result ~ n_patients, data = dtf, family = binomial())
  # predictions <- predict(mod, dtf, type = 'response')
  # library(ROCR)
  # pred <- prediction(predictions, dtf$exp_result)
  # perf <- performance(pred, 'acc')
  # which.max(perf@x.values[[1]])
  # perf@x.values[[1]][which.max(perf@y.values[[1]])]
  # plot(perf)

  # pacman::p_load('MLmetrics')
  # boxplot(as.integer(dtf$exp_result), predictions)
  # predict(mod, dtf, type = 'response')
  # hist(predictions, breaks = 100)
  # threshold = 0.05
  # F1_Score(
  #   y_true = dtf$exp_result,
  #   y_pred = as.numeric(predictions > threshold))
  # plot(glm(exp_result ~ n_patients, data = dtf, family = binomial()))

  if (bayesian_analysis && !is.null(plot_fishtails)) {
    plot_id <- paste0(
      # tumor_types[unique(dtf$project_extended)], '-',
      gsub(' ', '_', tolower(stage_id)),
      paste0(sort(unique(as.character(dtf$focus_allele))), collapse = '_'),
      collapse = '_')
    p <- plot_bayesian_error_vs_effect(dtf, grouping = plot_fishtails)
    o_fn <- file.path(img_loc,
      glue::glue('bayesian_cont_IE-error_vs_effect-{plot_id}.png'))
    if (plot_fishtails != 'project') {
      h <- 7; w = 17.4;
    } else {
      h <- 25; w = 2 * 17.4;
    }
    ggsave(o_fn, p, height = h, width = w, unit = 'cm')
  }
  return(invisible())
}


assess_coef_composition <- function(dtf) {
  grp_vars <- c('focus_allele', 'VE_threshold',
    'patient_inclusion_crit', 'expression_threshold', 'LOH_HLA',
    'analysis_name', 'sts_filtering', 'percentile_rank',
    'overlap_var', 'project_extended') %>%
    intersect(colnames(dtf))

  out <- purrr::map(auto_name(grp_vars), ~dtf[, .N, by = .x])
  return(out)
}


test_dplyr <- function() {
  res <- tryCatch({
    by_cyl <- mtcars %>% dplyr::group_by(cyl)
    by_cyl %>% dplyr::filter(disp == max(disp))
  }, error = function(e) NULL)
  # res <- tryCatch(stop(simpleError('a')), error = function(e) NULL)
  !is.null(res)
}
stopifnot(test_dplyr())


hmp_wrapper <- function(v) {
  pacman::p_load('harmonicmeanp')
  p.hmp(v, L = length(v))
}


#' Summarize analyses of multiple alleles to one representative allele
#'
#' Sort the alleles with respect to the \code{sort_var} and slice out
#' the middle one. The different methods pertain to different ways of
#' dealing with cases where there are an even number of alleles to
#' pick a representative allele from.
#'
#' aggressive: pick the lowest (most consistent with
#' immunoediting/neo-antigen depletion) careful: select middle two
#' analyses and average their results conservative: pick the highest
#' (most consistent with neo-antigen enrichment)
#'
pick_representative_allele <- function(dtf,
  sort_var = 'depletion_full', method = 'DT_careful') {

  dtf <- format_overview_res(dtf)

  if (!sort_var %in% colnames(dtf)) {
    mystop(msg = glue('{sort_var} not in input'))
  }

  analysis_id_cols <- colnames(all_cont_IE_settings) %>%
    setdiff('focus_allele') %>%
    c('analysis_idx', 'project_extended') %>%
    unique()

  if (grepl('dplyr', method)) {
    sub_method <- gsub('dplyr_', '', method)

    if (sub_method == 'conservative') {
      get_proper_idx <- function(N) {
        switch(N, '1' = 1, '2' = 2, '3' = 2, '4' = 3, '5' = 3)
      }
    } else if (sub_method == 'aggressive') {
      get_proper_idx <- function(N) {
        switch(N, '1' = 1, '2' = 1, '3' = 2, '4' = 2, '5' = 3)
      }
    }

    ## arrange sorts in ascending fashion, so this returns the 3rd item for a
    ## group of 5 but the 2nd item for a group of 4
    rc_sel <- dtf %>%
      group_and_count_alleles(sort_var = sort_var) %>%
      dplyr::arrange(.data[[sort_var]], .by_group = T) %>%
      dplyr::filter(!is.na(.data[[sort_var]])) %>%
      dplyr::filter(row_number() == get_proper_idx(n()))
  } else if (method == 'DT_careful') {
    setDT(dtf)
    rc_sel <-
      dtf %>%
      .[!is.na(get(sort_var))] %>%
      .[order(get(sort_var)),
        cbind(.SD, N_alleles = .N),
        by = analysis_id_cols] %>%
      .[, i := 1:.N, by = analysis_id_cols]

    if ('opt_string' %in% colnames(rc_sel)) {
      ## Strip allele away from opt string, since we're integrating it
      ## out
      rc_sel[, opt_string :=
        gsub('[^ ]*\\s{1}(.*)', '\\1', opt_string)]
    }
    rc_sel_even <- rc_sel[N_alleles %% 2 == 0]
    rc_sel_uneven <- rc_sel[N_alleles %% 2 == 1]

    ## Seperate processing for analyses with an even and uneven number of
    ## alleles
    if (nrow(rc_sel_even) > 0) {
      numeric_columns <- purrr::map_chr(rc_sel, class) %>%
        purrr::keep(~.x %in% c('numeric', 'integer')) %>%
        names
      p_val_columns <- grep('p_val', numeric_columns, value = T)
      numeric_columns <- setdiff(numeric_columns,
        c(p_val_columns, 'i', 'analysis_idx'))

      ## Select the middle two rows for each group
      rc_sel_even <- rc_sel_even[i %in% c(N_alleles / 2, N_alleles / 2 + 1)]

      ## Integrate the middle two analyses
      ## Simply average numeric, non p-value columns;
      ## HMP average p-values
      df_a <- rc_sel_even[, lapply(.SD, mean),
        by = analysis_id_cols, .SDcols = c(numeric_columns)]
      df_b <- rc_sel_even[, lapply(.SD, hmp_wrapper),
        by = analysis_id_cols, .SDcols = p_val_columns]
      rc_sel_even <- merge(df_a, df_b, by = analysis_id_cols)
    }

    if (nrow(rc_sel_uneven) > 0) {
      rc_sel_uneven <-
        rc_sel_uneven[i == ceiling(N_alleles / 2), .SD,
          by = analysis_id_cols]
    }

    rc_sel <- rbind(rc_sel_even, rc_sel_uneven, fill = T)
    rc_sel[, i := NULL]
    rc_sel[, focus_allele := NULL]
    na_cols <- map_lgl(rc_sel, ~all(is.na(.x))) %>%
      which %>% names
    if (length(na_cols) > 0) {
      rc_sel[, !na_cols]
    }
  } else if (method == 'DT_aggressive') {
    ## This method will pick the 2nd with n == 4 and 3rd with n == 5
    setDT(dtf)
    rc_sel <-
      dtf[, cbind(.SD[order(get(sort_var))][max(ceiling(.N / 2), 1)],
        'N_alleles' = .N),
          keyby = analysis_id_cols]
  }
  N_unique_alleles <- length(unique(dtf$focus_allele))

  if (any((rc_sel %>% pull(N_alleles) %>% max()) > N_unique_alleles)) {
    stop('Some IE analyses seem to have been summarized incorrectly')
  }

  setDT(rc_sel)
  return(rc_sel)
}


group_and_count_alleles <- function(dtf, sort_var = 'depletion_max') {
  res <- dplyr::group_by(dtf,
    overlap_var, patient_inclusion_crit, LOH_HLA,
    analysis_name, analysis_idx, project_extended) %>%
    dplyr::filter(!is.na(.data[[sort_var]])) %>%
    dplyr::mutate(N_alleles = uniqueN(focus_allele))
  return(res)
}


if (F) {
  test_dat <-
  structure(list(overlap_var = c("mean_score", "mean_score", "mean_score",
  "mean_score", "mean_score"), patient_inclusion_crit = c("", "",
  "", "", ""), LOH_HLA = c("LOHHLA", "LOHHLA", "LOHHLA", "LOHHLA",
  "LOHHLA"), analysis_name = c("clon_param_titration", "clon_param_titration",
  "clon_param_titration", "clon_param_titration", "clon_param_titration"
  ), analysis_idx = c(1L, 1L, 1L, 1L, 1L), project_extended = c("Bladder",
  "Bladder", "Bladder", "Bladder", "Bladder"), focus_allele = c("A0201",
  "A1101", "B0702", "B2705", "B4001"), rc = c(0.0037288312085828,
  -0.0022698896321099, 0.00259134228709325, 0.000380860442528955,
  0.00113542117110707), p_val = runif(5)), row.names = c(NA, -5L), groups = structure(list(
      overlap_var = "mean_score", patient_inclusion_crit = "",
      LOH_HLA = "LOHHLA", analysis_name = "clon_param_titration",
      analysis_idx = 1L, project_extended = "Bladder", .rows = structure(list(
          1:5), ptype = integer(0), class = c("vctrs_list_of",
      "vctrs_vctr", "list"))), row.names = 1L, class = c("tbl_df",
  "tbl", "data.frame"), .drop = TRUE), class = c("grouped_df",
  "tbl_df", "tbl", "data.frame"))
  altered_median <- pick_median(test_dat$rc)
  pick_median <- function(x) sort(x) %>% { .[ceiling(length(.) / 2)] }
  median(test_dat$rc, na.rm = T) == mean(sort(test_dat$rc)[c(2,3)])
  stopifnot(pick_representative_allele(test_dat, sort_var = 'rc',
      method = 'dplyr_aggressive')$rc == altered_median)
  stopifnot(pick_representative_allele(test_dat, sort_var = 'rc',
      method = 'dplyr_aggressive')$N_alleles == 4)


  test_dat <-
  structure(list(overlap_var = c("mean_score", "mean_score", "mean_score",
  "mean_score", "mean_score"), patient_inclusion_crit = c("", "",
  "", "", ""), LOH_HLA = c("LOHHLA", "LOHHLA", "LOHHLA", "LOHHLA",
  "LOHHLA"), analysis_name = c("clon_param_titration", "clon_param_titration",
  "clon_param_titration", "clon_param_titration", "clon_param_titration"
  ), analysis_idx = c(1L, 1L, 1L, 1L, 1L), project_extended = c("Bladder",
  "Bladder", "Bladder", "Bladder", "Bladder"), focus_allele = c("A0201",
  "A1101", "B0702", "B2705", "B4001"), rc = c(0.0037288312085828,
  -0.0022698896321099, 0.00259134228709325, 0.000380860442528955,
  0.00113542117110707), p_val = runif(5)), row.names = c(NA, -5L), groups = structure(list(
      overlap_var = "mean_score", patient_inclusion_crit = "",
      LOH_HLA = "LOHHLA", analysis_name = "clon_param_titration",
      analysis_idx = 1L, project_extended = "Bladder", .rows = structure(list(
          1:5), ptype = integer(0), class = c("vctrs_list_of",
      "vctrs_vctr", "list"))), row.names = 1L, class = c("tbl_df",
  "tbl", "data.frame"), .drop = TRUE), class = c("grouped_df",
  "tbl_df", "tbl", "data.frame"))
  vanilla_median <- pick_median(test_dat$rc)
  stopifnot(pick_representative_allele(test_dat, sort_var = 'rc',
      method = 'dplyr_aggressive')$rc == vanilla_median)

  # pick_representative_allele(test_dat, sort_var = 'rc', method = 'DT_careful')
}


cont_IE_fn_to_args <- function(
  fn = paste0('cont-IE-twoD_sens_analysis-A0201-mean_score',
    '-no_LOHHLA-FDR1-rlm-not_z_normalized-38.rds')) {
  fn <- gsub('cont-IE-|.rds', '', fn)
  args <- setNames(as.list(strsplit(fn, '-')[[1]]),
    c('analysis_name', 'focus_allele', 'overlap_var', 'LOH_HLA',
      'patient_inclusion_crit', 'reg_method',
      'z_normalize', 'analysis_idx'))
  args$z_normalize <- args$z_normalize == 'not_z_normalized'
  args$analysis_idx <- as.integer(args$analysis_idx)
  return(args)
}


drop_high_VE <- function(dtf) {
  setDT(dtf)
  dtf <- dtf[VE_threshold != 'VE=5']
  dtf <- dtf[VE_threshold != '5']
  dtf[, VE_threshold := droplevels(VE_threshold)]
  return(dtf)
}
