IE_requires_computation <- mod_time_comparator(
  minimum_mod_time = '2021-4-02 06:58', verbose = TRUE)

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


#' Generate file name(s) to be generated in \code{test_continuous_IE} and
#' other functions
#'
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
  extension = 'rds',
  z_normalize = F,
  patient_inclusion_crit = '') {

  if (!exists('base_name') || missing(base_name) ||
      is.null(base_name) || is.na(base_name)) {
    base_name <- ''
  }

  if (replace_idx) {
    analysis_idx <-
      complement_options(c(as.list(ds_param_grid[as.integer(analysis_idx), ]))) %>%
      gen_opts_string(no_spaces = T, no_hla = T)
  }
  analysis_idx <- prepend_hyphen(analysis_idx)
  analysis_name <- prepend_hyphen(analysis_name)
  focus_allele <- prepend_hyphen(focus_allele)
  hla_sim_range <- ifelse(is.null(hla_sim_range) || hla_sim_range == '',
    '', paste(hla_sim_range, collapse = '_')) %>%
    prepend_hyphen
  reg_method <- prepend_hyphen(reg_method)
  LOH_HLA <- prepend_hyphen(LOH_HLA)
  patient_inclusion_crit <- prepend_hyphen(patient_inclusion_crit)
  overlap_var <- prepend_hyphen(overlap_var)
  extension <- prepend_string(extension, prepend_string = '.')
  p_val_bound <- ifelse(eps(p_val_bound, 0.05), '',
    sprintf('-pvb-%s', format(p_val_bound, scientific = T, digits = 2)))
  fill_var <- ifelse(fill_var == 'rc', '', prepend_hyphen(fill_var))
  tumor_type <- ifelse(is.null(project_extended), '',
    tumor_types[project_extended])
  z_normalize <- ifelse(z_normalize == T, '', '-not_z_normalized')

  file_head <-
    tryCatch(paste0(base_name, tumor_type,
                    analysis_name, focus_allele, overlap_var, LOH_HLA,
                    patient_inclusion_crit, reg_method, p_val_bound,
                    fill_var, hla_sim_range, z_normalize, analysis_idx,
                    extension),
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
    data_subs <- dtf %>%
      dplyr::filter(project_extended == {{project_extended}})
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
    'CYT',
    'mean_score', 'mean_score_AB', 'mean_score_C', 'weight',
    'IE_essentiality_impaired', 'CYT', 'hla_allele_status',
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
determine_continuous_IE_yval <- function(dtf) {
  ## c_name will 'c_muts.missense_mutation' for regular
  ## donor_summaries
  c_name <- grep('^c.*', names(dtf), value = T) %>%
    { grep('silent_mutation', ., invert = T, value = T) } %>%
    .[1]
  i_name <- gsub('c_(.*)', 'i_\\1', c_name)
  if (!is.na(c_name) && c_name != 'NA') {
    dtf$y_var <- dplyr::pull(dtf, c_name) / dplyr::pull(dtf, i_name)
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


test_data_sufficiency <- function(dtf,
  min_pts = 20,
  min_points_lq = 3,
  min_points_uq = 3) {

  if (is.null(dtf)) return(F)
  setDT(dtf)
  discretized_ps <- dtf[, .N, cut(ol, breaks = seq(0, 1, by = .25))]

  lq <- '(0,0.25]'
  uq <- '(0.75,1]'

  ## st := simple test
  st <- function(v) !is.null(v) && length(v) > 0 && v
  out <-
    !is.null(dtf) &&
    nrow(dtf) >= min_pts &&
    st(with(discretized_ps, N[cut == lq] > min_points_lq)) &&
    st(with(discretized_ps, N[cut == uq] > min_points_uq)) &&
    any(dtf$y_var > 0)
  return(out)
}


#' Take in a \code{donor_summary}(-like) object and compute the relationship
#' between yield rate and HLA repertoire overlap
#'
#'
compute_continuous_IE_statistics <- function(
  dtf,
  reg_method = 'rlm',
  focus_allele = 'A0201',
  partition_vars = c(),
  bayesian_args = list(N_iter = 4000),
  debug = F) {

  if (null_dat(dtf)) return(NULL)

  ## Partition dtf in all relevant subgroups
  p_dtf <- partition_dtf(
    dtf = dtf,
    partition_vars = partition_vars
  )
  p_dtf <- p_dtf[purrr::map_lgl(p_dtf, test_data_sufficiency)]
  if (length(p_dtf) == 0) return(NULL)

  if (reg_method %in% c('rlm_one_group', 'rlm')) {
    f <- fit_rlm_model
  } else if (grepl('bayesian', reg_method)) {
    # ls(brm_scaffold)
    ## Example input: reg_method = 'bayesian_unbiased', strip off
    ## 'bayesian_'
    model_name <- gsub('[^_]*_(.*)', '\\1', reg_method)
    if (model_name == 'bayesian') {
      model_name <- 'biased'
    }
    f <- suppressMessages(pryr::partial(fit_bayesian_regression,
        model_name = model_name))
  } else if (reg_method == 'lm') {
    f <- lm_fit_model
  } else if (reg_method == 'gam') {
    f <- gam_fit_model
  }
  lm_results <- purrr::map(p_dtf, f)

  ## Append Wilcoxon-test results
  wilcox_results <- purrr::map(p_dtf, fit_wilcox_model)
  out <- purrr::map(seq(p_dtf),
    ~c(lm_results[[.x]], wilcox_results[[.x]])
  )
  return(out)
}


fit_wilcox_model <- function(dtf) {
  if (is.null(dtf)) return(NULL)
  dtf[, 'ol_quartile' :=  cut(ol, breaks = seq(0, 1, by = .25))]
  dtf[as.integer(ol_quartile) %in% c(1, 4)]
  wc <- suppressWarnings(wilcox.test(y_var ~ ol_quartile,
    data = dtf[as.integer(ol_quartile) %in% c(1, 4)]))
  return(list(
    'wilcox_stat' = unname(wc$statistic),
    'wilcox_p' = wc$p.value
  ))
}


#' Fit robust regression model with Huber error model
#'
#'
fit_rlm_model <- function(dtf) {
  if (!test_data_sufficiency(dtf)) return(NULL)

  unnormalized_mod <- fit_rlm_(dtf, simple_output = T)
  if (is.null(unnormalized_mod)) return(NULL)
  dtf$y_var <- dtf$y_var / coef(unnormalized_mod)[1]
  normalized_mod <- fit_rlm_(dtf, simple_output = F)
  # orig_pars <- unnormalized_mod[c('intercept', 'rc')] %>%
  #   { setNames(., paste('orig_', names(.), sep = '')) }
  # stopifnot(maartenutils::eps(
  #   1/unnormalized_mod$intercept * unnormalized_mod$rc,
  #   normalized_mod$rc, 1e-3))
  normalized_mod %>%
    modifyList(unnormalized_mod[c('intercept', 'rc')])
}


fit_rlm_ <- function(dtf, simple_output = F) {
  dtf <- setDT(dtf)[is.finite(y_var) & is.finite(ol)]

  lc <- tryCatch(suppressWarnings(
      MASS::rlm(y_var ~ ol, data = dtf,
                x.ret = T, model = T, maxit = 1000)),
    error = function(e) { NULL })
  if (simple_output) {
    return(lc)
  }

  quants <- c(.1, .25, .5, .75, .9)
  delta_names <- 
    paste0('delta_', c('CI_l', 'mean', 'CI_h'), sep = '') %>%
    c('delta_SE')
  NMAD_eval_locs <- seq(0.1, .5, by = .1)
  NMAD_names <- paste('NMAD_prob_', 10 * NMAD_eval_locs, sep = '')

  if (is.null(lc)) {
    AFDP_def <- map(quants, ~NA_real_) %>%
      setNames(paste0('AFDP_', 100*quants, sep = ''))
    delta_def <- map(auto_name(delta_names), ~NA_real_)
    NMAD_def <- map(auto_name(NMAD_names), ~NA_real_)
    return(list(
      'intercept' = NA_real_,
      'p_val_intercept' = NA_real_,
      'rc' = NA_real_,
      'p_val' = NA_real_,
      'df' = NA_real_,
      't_val' = NA_real_,
      'n_patients' = dtf[, .N],
      'converged' = F,
      'yr_fractional_change' = NA_real_,
      'scale' = NA_real_,
      'norm_scale' = NA_real_
    ) %>% c(AFDP_def) %>% c(delta_def) %>% c(NMAD_def)) 
  }

  ## Compute absolute fractional difference between predictions (AFDP)
  ## and observations
  non_zero_idx <- dtf[, y_var > 0 | y_var < 0]
  obs <- dtf[, y_var][non_zero_idx]
  pred <- predict(lc)[non_zero_idx]
  dev <- abs((pred - obs) / obs)
  AFDP_stats <- as.list(quantile(dev, quants)) %>%
    setNames(paste0('AFDP_', 100*quants, sep = ''))

  coef_mat <- summary(lc)$coefficients
  coef_mat <- cbind(coef_mat,
    'p_value' = sapply(seq_along(coef_mat[, 3]), function(i) {
      tryCatch(sfsmisc::f.robftest(lc, var = i)$p.value,
        error = function(e) { NA })
    }),
    'df' = summary(lc)$df[1:nrow(coef_mat)]
  )

  delta <- predict(lc, 
    newdata = data.frame(ol = 1), se.fit = T)
  delta <- as.list(delta$fit + c(-1.96, 0, 1.96) * delta$se.fit - 1) %>%
    append(list(delta$se.fit)) %>%
    setNames(delta_names)

  intercept <- predict(lc, 
    newdata = data.frame(ol = 0), se.fit = T)
  ## What is the probability that NMAD is smaller than or equal to the
  ## the upper boundaries defined in NMAD_eval_locs?
  NMAD_probs <- pnorm(NMAD_eval_locs, mean = lc$s / intercept$fit, 
    sd = intercept$se.fit) %>%
    as.list %>%
    setNames(NMAD_names)

  return(list(
    'lm' = lc,
    'intercept' = coef_mat[1, 1],
    'p_val_intercept' = coef_mat[1, 4],
    'rc' = coef_mat[2, 1],
    'p_val' = coef_mat[2, 4],
    'df' = coef_mat[2, 5],
    't_val' = coef_mat[2, 3],
    'n_patients' = dtf[, .N],
    'converged' = lc$converged,
    'yr_fractional_change' = coef_mat[2, 1] / coef_mat[1, 1],
    'scale' = lc$s,
    'norm_scale' = lc$s / coef_mat[1, 1]
  ) %>% c(AFDP_stats) %>% c(delta) %>% c(NMAD_probs))
}


fit_lm_model <- function(dtf) {
  dtf <- dtf[is.finite(y_var) & is.finite(ol)]

  lc <- tryCatch(lm(y_var ~ ol, data = dtf),
    error = function(e) { print(e); NULL })
  coef_mat <- summary(lc)$coefficients
  if (is.null(lc) || dim(coef_mat)[1] == 1) return(NULL)

  res <- list(
    'lm' = lc,
    'intercept' = coef_mat[1, 1],
    'p_val_intercept' = coef_mat[1, 4],
    'rc' = tryCatch(coef_mat[2, 1], error = function(e) { browser() }),
    'p_val' = coef_mat[2, 4],
    't_val' = coef_mat[2, 3],
    'n_patients' = dtf[, .N])

  return(res)
}


gam_fit_model <- function(dtf) {
  dtf <- dtf[is.finite(y_var) & is.finite(ol) &
    is.finite(CYT) & is.finite(weight)]

  gam_fit <- tryCatch(mgcv::gam(y_var ~ s(ol), bs = 'cs',
      data = dtf), error = function(e) { NULL })

  if (is.null(gam_fit)) return(NULL)

  preds <- predict(gam_fit, data.frame(ol = c(0, 1)), se.fit = T) %>%
    as.data.frame

  intercept <- tryCatch(preds[1, 1],
    error = function(e) { print(e); browser() })

  rc <- tryCatch(preds[2, 1] - preds[1, 1],
    error = function(e) { browser() })

  ## Compute probability that rc is of the other sign than its estimate.
  ## I.e. P(e > 0) if e < 0
  ## TODO: find better/more principled way of combining error measurements.
  ## Summing SEs is only correct for summation and subtraction of random
  ## variables
  p_val <- pnorm(0, mean = rc, sd = sum(preds[, 2]))
  # ## Compute probability that rc is of the other sign than its estimate.
  # ## I.e. P(e > 0) if e < 0
  # ## From Wikipedia: https://en.wikipedia.org/wiki/Ratio_distribution
  # ## Díaz-Francés, Eloísa; Rubio, Francisco J. (2012-01-24). "On the existence of a normal approximation to the distribution of the ratio of two independent normal random variables". Statistical Papers. Springer Science and Business Media LLC. 54 (2): 309–323. doi:10.1007/s00362-012-0429-2. ISSN 0932-5026.
  if (rc < 0) {
    p_val <- 1 - p_val
  }

  res <- list(
    'lm' = gam_fit,
    'intercept' = intercept,
    'rc' = rc,
    'p_val' = p_val,
    'n_patients' = dtf[, .N])

  return(res)
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


recover_tumor_type <- function(tumor_type = NULL, project_extended = NULL) {
  if (!is.null(tumor_type) && tumor_type %in% tumor_types) {
    return(tumor_type)
  } else if (is.null(tumor_type) && !is.null(project_extended)) {
    return(tumor_types[as.character(project_extended)])
  }
}
# recover_tumor_type(project_extended = 'Rectum')


#' Master function for continuous IE detection. Perform and plot IE analyses for
#' all samples combined and for each level of column \code{project_extended}
#' separately
#'
#' @param reg_plot Include regression plots of the combined cohort and all
#' cohorts (levels of \code{project_extended}) separately
#'
test_continuous_IE <- function(
  analysis_name = 'twoD_sens_analysis',
  project_extended = NULL,
  analysis_idx = 1,
  focus_allele = 'A0201',
  redo = F,
  return_grob = F,
  reg_method = 'rlm',
  patient_inclusion_crit = 'none',
  ncores = 1,
  check_res = F,
  partition_vars = c(),
  stats_idx = 1,
  hla_sim_range = NULL,
  z_normalize = F,
  LOH_HLA = 'no_LOHHLA',
  verbose = F,
  overlap_var = 'mean_score',
  include_call = F,
  ds_frac = NULL,
  return_res = T) {

  ## No point in doing an intercept vs. slope comparison as is being
  ## done in the current Bayesian analysis when z_normalize == TRUE
  if (z_normalize && reg_method == 'bayesian') {
    stop('Z normalization and the Bayesian reg_method should',
      ' not be used in together')
  }

  if (include_call) 
    f_args <- as.list(environment())

  o_fn <- gen_cont_IE_fn(
    base_name = 'cont-IE',
    focus_allele = focus_allele,
    analysis_name = analysis_name,
    reg_method = reg_method,
    hla_sim_range = hla_sim_range,
    project_extended = project_extended,
    sub_dir = 'cont-IE-test',
    p_val_bound = NULL,
    analysis_idx = analysis_idx,
    z_normalize = z_normalize,
    LOH_HLA = LOH_HLA,
    patient_inclusion_crit = patient_inclusion_crit,
    overlap_var = overlap_var
  )

  ## Don't cache downsampled sub-analyses
  if (is.null(ds_frac) && !IE_requires_computation(o_fn) && 
      !redo && !check_res) {
    if (return_res) {
      ret_val <- tryCatch(
        readRDS(o_fn),
        error = function(e) { file.remove(o_fn); NULL })
      if (!is.null(ret_val)) return(ret_val)
    } else {
      return(NULL)
    }
  }

  prep <- prep_cont_IE_analyses(
    focus_allele = focus_allele,
    redo = redo,
    analysis_name = analysis_name,
    project_extended = project_extended,
    hla_sim_range = hla_sim_range,
    z_normalize = z_normalize,
    LOH_HLA = LOH_HLA,
    analysis_idx = analysis_idx,
    overlap_var = overlap_var,
    patient_inclusion_crit = patient_inclusion_crit,
    ds_frac = ds_frac
  )

  if (is.null(prep)) return(NULL)

  if (null_dat(prep$dtf)) {
    ret <- list(
      'analysis_name' = glue::glue('{focus_allele}_{analysis_name}'),
      'stats' = NA
    )
  } else {
    stats_by_project <- plyr::llply(auto_name(prep$projects), 
      function(pe) {
      data_subs <- subset_project(prep$dtf, pe)
      if (null_dat(data_subs)) return(NULL)
      tryCatch({
        compute_continuous_IE_statistics(
          dtf = data_subs,
          reg_method = reg_method,
          focus_allele = focus_allele,
          partition_vars = partition_vars
        )
      }, error = function(e) { print(e); NULL })
    }, .parallel = (ncores > 1))

    ret <- list(
      'analysis_name' = 
        (attr(prep$dtf, 'analysis_name') %||%
          glue::glue('{focus_allele}_{analysis_name}')),
      'stats' = stats_by_project
    )
  }

  if (is.null(ds_frac) && !check_res) {
    saveRDS(ret, o_fn)
  } else if (check_res) {
    stored_res <- readRDS(o_fn)
    file.mtime(o_fn)
    new_res <- coef(ret[['stats']]$`Pan*'-'*cancer`[[1]]$lm)
    old_res <- coef(stored_res[['stats']]$`Pan*'-'*cancer`[[1]]$lm)
    return(all(eps(new_res, old_res, 1e-16)))
  }

  if (include_call) {
    attr(ret, 'call') <- f_args
    attr(ret, 'prep') <- prep
  }
  return(ret)
}


#'
#'
#' @param test_object Output from test_continuous_IE()
wrapper_plot_repertoire_overlap_vs_yield_rate <- function(
  test_object,
  colour_vars = c(),
  projects = attr(test_object, 'prep')$projects,
  shape_var = NULL,
  point_alpha = .8,
  stats_idx = 1,
  ncores = 1,
  return_grob = F, ...) {

  if (is.null(attr(test_object, 'call')) ||
      is.null(attr(test_object, 'prep'))) {
    stop('Compute test object with include_call = TRUE')
  }

  reg_plots <- plyr::llply(maartenutils::auto_name(projects), function(pe) {
    data_subs <- subset_project(attr(test_object, 'prep')$dtf, pe,
      pan_can_includes_all = TRUE)
    if (null_dat(data_subs)) browser()
    tryCatch(plot_PS_vs_yr(
        dtf = data_subs,
        reg_method = attr(test_object, 'call')$reg_method,
        return_grob = return_grob,
        shape_var = shape_var,
        hla_allele = attr(test_object, 'call')$focus_allele,
        point_alpha = point_alpha,
        colour_vars = colour_vars,
        partition_vars = attr(test_object, 'call')$partition_vars,
        lm = test_object$stats[[pe]][[stats_idx]][['lm']], ...),
      error = function(e) { print(e) })
  }, .parallel = (ncores > 1))
}



#'  Prepare data for multiple types of downstream analyses. Filter patients and
#'  prepare the required regression variables.
#'
#'
prep_cont_IE_analyses <- function() {
  ## Determine donor_summary file to load in
  obj_fn <- gen_cont_IE_ds_fn(
    hla_allele = focus_allele,
    analysis_name = analysis_name,
    analysis_idx = analysis_idx,
    overlap_var = overlap_var
  )
  if (is.null(obj_fn) || is.na(obj_fn) || length(obj_fn) == 0)
    return(NULL)

  ## Load in the data and compute the ordinal variable, save attributes
  dtf <- readRDS(obj_fn)

  # analysis_name = 'rooney_param_titration'
  if (analysis_name == 'rooney_param_titration') {
    cond_setnames(dtf, 'analyis_name', 'analysis_name')
    dtf <- dtf[grepl('Rooney UMF', analysis_name)]
  }

  attr(dtf, 'y_type') <- if (!grepl('rooney', analysis_name))
    'neo-antigen yield rate' else 'observed/expected neo-antigens'
  attr(dtf, 'x_label') <- sprintf('%s presentation score',
    quickMHC::ppHLA(focus_allele))
  attr(dtf, 'y_label') <- sprintf('%s %s',
    quickMHC::ppHLA(focus_allele), attr(dtf, 'y_type'))

  atts <- attributes(dtf)
  dtf <- determine_continuous_IE_yval(dtf)

  ## Extract the  tumor projects present in this object
  dtf <- subset_project(dtf, project_extended)

  ## Downsample patients
  if (!is.null(ds_frac)) {
    dtf <- dtf[sample(1:nrow(dtf), ceiling(nrow(dtf)*ds_frac), 
      replace = T), ]
  }

  ## Annotate with required additional variables
  if ('IE_essentiality_impaired' %in% partition_vars ||
      any(patient_inclusion_crit %in% c('strict_TR', 'TR'))) {
    dtf <- setDT(dtf)
    FDR_thresh <- switch(patient_inclusion_crit, 
      'strict_TR' = .01, 'TR' = .1
    )
    dtf <- 
      annotate_donor_by_IE_essentiality(dtf, 
        FDR_thresh = FDR_thresh, verbose = verbose)
    dtf[, mean(IE_essentiality_impaired)]
  }
  if (any(patient_inclusion_crit %in% c('strict_TR', 'TR'))) {
    dtf <- dtf[IE_essentiality_impaired == F]
  }

  ## Create column name "ol" for presentation score
  repertoire_overlap_dat <-
    get_repertoire_overlap(hla_alleles = focus_allele)
  HLA_LOH_l <- ifelse(grepl('no', LOH_HLA), 'no_LOHHLA', 'LOHHLA')
  ro_s <- repertoire_overlap_dat[LOH_HLA == HLA_LOH_l]
  dtf <- tryCatch(
    merge_hla(dtf, 
      focus_allele = as.character(focus_allele),
      repertoire_overlap_dat = ro_s,
      overlap_var = as.character(overlap_var)
    ), error = function(e) { print(e); browser() })
  dtf <- dtf[!is.na(ol)]

  ## Annotate CYT, currently not used in model
  dtf <- merge_CYT(dtf)

  ## Subset away non-completely HLA-loss annotated pts if required
  if ('LOHHLA_complete' %in% colnames(repertoire_overlap_dat) &&
      LOH_HLA == 'strict_LOHHLA') {
    dtf <- maartenutils::controlled_merge(dtf,
      repertoire_overlap_dat[hla == focus_allele,
        .(donor_id, LOH_HLA, LOHHLA_complete, LOHHLA_AB_complete)],
      by_cols = 'donor_id', maintain_attr = 'ds_opts')
    if (overlap_var == 'mean_score') {
      dtf <- dtf[LOHHLA_complete == T]
    } else if (overlap_var == 'mean_score_AB') {
      dtf <- dtf[LOHHLA_AB_complete == T]
    }
  }

  ## Trim to pts with valid dependent variable (y_var)
  dtf <- dtf[!is.na(y_var) & y_var != Inf]
  if (!is.null(hla_sim_range) && length(hla_sim_range) == 2) {
    dtf <- dtf[ol >= hla_sim_range[1] & ol <= hla_sim_range[2]]
  }
  if (null_dat(dtf) || all(is.na(dtf$y_var)) || all(is.na(dtf$ol))) {
    return(NULL)
  }

  ## Reduce object to crucial columns only in order to make it fit in
  ## memory more easily
  dtf <- reduce_to_essential_IE_columns(dtf)

  ## Z-normalize if so desired (no!)
  if (z_normalize) dtf <- normalize_columns(dtf)

  ## Restore all attributes the data.table operations have removed
  ## (sigh)
  for (att in setdiff(names(atts),
       c('class', 'row.names', '.internal.selfref', 'sorted',
         'data.table', 'names'))) {
    attr(dtf, att) <- atts[[att]]
  }

  observed_projects <-
    sort(as.character(unique(dtf$project_extended)))
  observed_projects <-
    c(sort(as.character(unique(dtf$project_extended))),
      'Pan*\'-\'*cancer')
  out <- list(
    'projects' = observed_projects,
    'dtf' = dtf,
    'obj_fn' = obj_fn
  )
  return(out)
}
formals(prep_cont_IE_analyses) <- formals(test_continuous_IE)
formals(prep_cont_IE_analyses)$min_points_lq <-
  formals(prep_cont_IE_analyses)$min_points_uq <- 3


#' Compile overview for the combination of a set of IE analysis
#' settings and a single focus allele. Rows correspond to different
#' neo-antigen pipeline settings
#'
#' @value \code{data.table} Containing one row per combination of all
#' immunoediting settings (minus \code{focus_allele}) and a tumor type
#' for the focus allele specified with parameter \code{focus_allele}.
#'
prep_continuous_param_grid <- function(
  focus_allele = 'A0201',
  prep = NULL,
  reg_method = 'rlm',
  LOH_HLA = F,
  patient_inclusion_crit = '',
  analysis_name = 'twoD_sens_analysis',
  proj_order = NULL,
  analysis_idxs = main_analysis_idxs,
  redo = T,
  # hla_sim_range = NULL,
  legend_pos = 'bottom',
  ncores = 16,
	z_normalize = F,
  ds_frac = NULL,
  stats_idx = 1,
  return_res = F,
  overlap_var = 'mean_score',
  plot_ranges = NULL,
  plot_sig_var = 'signif_effect.adj',
  ...) {

  dtf <- plyr::llply(analysis_idxs, function(analysis_idx) {
    dtf <- test_continuous_IE(
      focus_allele = focus_allele,
      reg_method = reg_method,
      return_res = T,
      redo = redo,
      LOH_HLA = LOH_HLA,
      ds_frac = ds_frac,
      z_normalize = z_normalize,
      analysis_name = analysis_name,
      analysis_idx = analysis_idx,
      patient_inclusion_crit = patient_inclusion_crit,
      overlap_var = overlap_var
    )
    if (null_dat(dtf)) return(NULL)
    projects <- names(dtf$stats)
    if (length(projects) == 0) { return(NULL) }
    if (is.null(dtf) || is.null(dtf$stats) || !is.list(dtf$stats)) 
      return(NULL)

    ## Extract stats
    ## lol := list of lists
    lol <- purrr::map(dtf$stats, `[[`, stats_idx)
    good_idx <- {
      !sapply(lol, is.null) & 
      sapply(lol, function(x) !is.null(names(x)))
    } %>% which()
    if (is.null(good_idx)) return(NULL)
    lol <- lol[good_idx]
    if (is.null(lol)) {
      return(NULL)
    }

    ## `test_continuous_IE` can return bare Bayesian model fits or
    ## lists of summary stats already extracted from these fits. Both
    ## scenarios are accounted for here.
    if (all(purrr::map(lol, class) == 'brmsfit')) {
      if (T) {
        o_fn <- gen_cont_IE_fn(
          base_name = 'cont-IE',
          focus_allele = focus_allele,
          analysis_name = analysis_name,
          reg_method = reg_method,
          sub_dir = 'cont-IE-test',
          p_val_bound = NULL,
          analysis_idx = analysis_idx,
          z_normalize = z_normalize,
          LOH_HLA = LOH_HLA,
          patient_inclusion_crit = patient_inclusion_crit,
          overlap_var = overlap_var
        )
        warning(o_fn)
        file.remove(o_fn)
        return(NULL)
      }
      IE_tests <- purrr::map(lol, perform_bayesian_IE_test)
    } else {
      IE_tests <- purrr::map(lol, ~ .x[names(.x) != 'lm'])
    }

    ret_val <- tryCatch(rbindlist(IE_tests, fill = T) %>%
      .[, 'project_extended' := projects[good_idx]] %>%
      .[, 'analysis_idx' := analysis_idx],
      error = function(e) { print(e); NULL })
    return(ret_val)
  }, .parallel = (ncores > 1)) %>% rbindlist(fill = T)

  if (null_dat(dtf) || ncol(dtf) == 2) return(NULL)
  dtf[project_extended == 'Combined',
    project_extended := 'Pan*\'-\'*cancer']

  pipeline_param_titration_grid_dat <- dtf %>%
    .[, analysis_name := sapply(analysis_idx, function(l_idx)
        gen_opts_string(complement_options(
            c(as.list(ds_param_grid[l_idx, ]),
              'hla_allele' = focus_allele))))]

  f_stats_idx <- ifelse(stats_idx == 1, '', 
    glue::glue('-SI_{stats_idx}'))

  return(list(
    'pipeline_param_titration_grid_dat' =
      pipeline_param_titration_grid_dat,
    'f_stats_idx' = f_stats_idx))
}


format_overview_res <- function(
  dtf, 
  reg_method = NULL, 
  z_normalize = NULL,
  ds_frac = NULL,
  iter = NULL) {
  if ('patient_inclusion_crit' %in% colnames(dtf) &&
      is.list(dtf$patient_inclusion_crit)) {
    dtf$patient_inclusion_crit <- unlist(dtf$patient_inclusion_crit)
  }

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
    dtf[, 'z_normalize' := F]
  }

  if (!is.null(ds_frac)) {
    dtf[, 'ds_frac' := ds_frac]
  }

  if (!is.null(iter)) {
    dtf[, 'iter' := iter]
  }

  if (all(c('scale', 'intercept') %in% colnames(dtf))) {
    dtf[, 'norm_scale' := scale / intercept]
  }

  return(dtf)
}


#' Compile all continuous IE summary statistics into one object
#'
#'
compile_all_coef_overview <- function(
  redo = F,
  redo_subanalyses = F,
  ncores = 32,
  debug_func = F,
  reg_method = 'rlm',
  ds_frac = NULL,
  iter = NULL,
  analysis_idxs = main_analysis_idxs,
  z_normalize = F,
  hla_alleles = focus_hlas,
  include_non_ds_tallies = (idxs_name != 'main_analysis')) {

  fa_flag <- paste0('-focus_hlas_', paste(hla_alleles, 
      collapse = '_'))
  if (fa_flag == '-focus_hlas_A0201_A1101_B0702_B2705_B4001') {
    fa_flag <- ''
  }
  idxs_name <- attr(analysis_idxs, 'name') %||%
    gsub('_idxs$', '', deparse(substitute(analysis_idxs)))
  analysis_idxs_flag <- paste0('-analysis_idxs=', idxs_name)
  if (!is.null(ds_frac))
    stopifnot(!is.null(iter))
  all_coef_fn <- file.path(rds_dir,
    glue::glue('all_coef_overviews\\
      {fa_flag}\\
      {make_flag(z_normalize)}\\
      {analysis_idxs_flag}\\
      {make_flag(reg_method)}\\
      {make_flag(ds_frac)}\\
      {make_flag(iter)}\\
      {make_flag(include_non_ds_tallies)}\\
      .rds'))

  if (!IE_requires_computation(all_coef_fn) && !redo) {
    return(
      format_overview_res(readRDS(all_coef_fn),
        reg_method = reg_method,
        ds_frac = ds_frac,
        iter = iter,
        z_normalize = z_normalize
      )
    )
  }

  message('Compiling all_coef_overview')
  all_settings <-
    as.list(all_cont_IE_settings) %>%
    purrr::map(unique) %>%
    { .[names(.) != 'focus_allele'] } %>%
    { purrr::exec(tidyr::expand_grid, !!!., 
        focus_allele = hla_alleles) }

  if (!include_non_ds_tallies) {
    all_settings %<>% filter(analysis_name == 'twoD_sens_analysis')
  }

  ## Make sure donor summaries are available for all hla_alleles and
  ## analysis_idxs
  cat('Preparing all tallies:\n')
  tidyr::expand_grid(
    hla_allele = hla_alleles,
    analysis_idx = analysis_idxs) %>%
    plyr::a_ply(1, function(r) {
      l_opts <- default_opts
      l_opts[['verbose']] <- F
      l_opts[['parallel']] <- F
      l_opts[['hla_allele']] <- r[['hla_allele']]
      create_donor_summary(
        analysis_idx = r[['analysis_idx']],
        base_opts = l_opts,
        requires_computation = ds_requires_computation,
        ncores = ncores,
        precompute = F,
        hla_allele = r[['hla_allele']],
        redo = F,
        donor_ids = av_donors
      )

      if (include_non_ds_tallies) {
        gen_non_ds_tallies(
          hla_alleles = r[['hla_allele']],
          analysis_idxs = r[['analysis_idx']],
          ncores = ncores,
          requires_computation = ds_requires_computation,
          verbose = T
        )
      }
    }, .parallel = F, .progress = 'text')

  cat('Performing all regressions:\n')
  i <- 0
  ## Compile overview of all analyses for these patient selection
  ## criteria
  all_coef_overviews <- plyr::adply(all_settings, 1, function(r) {
    if (debug_func) {
      i <<- i + 1
      print(i)
      if (i == 1) browser() else return(NULL)
    }
    patient_inclusion_crit <- r[['patient_inclusion_crit']] %>%
      { as.character(unlist(.)) }

    ## Create single grid for all tumor types and this particular set
    ## of IE settings
    dtf <- prep_continuous_param_grid(
      focus_allele = r[['focus_allele']],
      LOH_HLA = r[['LOH_HLA']],
      analysis_idxs = analysis_idxs,
      # hla_sim_range = hla_sim_range,
      redo = redo_subanalyses,
      patient_inclusion_crit = patient_inclusion_crit,
      ds_frac = ds_frac,
      ncores = ncores,
      reg_method = reg_method,
      simplify_pvals = F,
      z_normalize = z_normalize,
      analysis_name = r[['analysis_name']],
      overlap_var = r[['overlap_var']],
      plot_sig_var = plot_sig_var,
      project_ranking = 'sum_fill_var'
    )

    if (is.null(dtf)) {
      return(NULL)
    }

    ret_val <- dtf$pipeline_param_titration_grid_dat
    setnames(ret_val, 'CYT_rc', 'rc_CYT', skip_absent = T)
    setnames(ret_val, 'analysis_name', 'opt_string',
      skip_absent = T)
    setnames(ret_val, 'l_analysis_name', 'analysis_name',
      skip_absent = T)
    return(ret_val)
  }, .parallel = F, .inform = F, .progress = 'text')

  saveRDS(all_coef_overviews, all_coef_fn)
  mymessage(
    msg = sprintf('Wrote all_coef_overviews to %s', all_coef_fn),
    instance = 'compile_all_coef_overview'
  )
  return(
    format_overview_res(
      all_coef_overviews,
      iter = iter,
      ds_frac = ds_frac,
      reg_method = reg_method,
      z_normalize = z_normalize
    )
  )
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


print_overview_stats <- function(dtf, stage_id = '',
  plot_fishtails = NULL) {
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
    dtf[, 'exp_result' := intercept > 0 & rc < 0]
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

  ## This assumes that all levels are observed
  N_possible_analyses <- dtf %>%
    dplyr::select(focus_allele, overlap_var, patient_inclusion_crit,
                  project_extended, LOH_HLA, analysis_name, analysis_idx) %>%
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

  if ('yr_fractional_change' %in% colnames(dtf)) {
    cat('\nSummary of yr_fractional_change\n')
    print(dtf[, summary(yr_fractional_change)])
  }

  if ('rc' %in% colnames(dtf)) {
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

  if (!is.null(latest_stats)) {
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


#' V2 of regression coefficient heatmap, which is the average of
#' results for multiple focus alleles and is restricted to one tumor
#' type. After summarization, statistical power is optionally computed
#' for each individual analyis.
#'
#'
prep_pan_IE_heatmap <- function(
  tumor_type = 'Melanoma',
  redo = F,
  hla_allele = 'integrated',
  fill_var = 'depletion_full',
  ncores = 1,
  return_res = F,
  include_rooney = T,
  sort_var = 'depletion_full',
  p_val_bound = .25,
  include_non_ds_tallies = T,
  PPV = .45,
  z_normalize = F,
  analysis_idxs = main_analysis_idxs,
  es = '0.8',
  force_old = F,
  include_power = F,
  reg_method = 'rlm',
  use_PS_accuracy = F,
  integration_method = 'union',
  pick_allele_method = 'DT_careful',
  filter_intercept = 'positive',
  pp_hl_threshold = .01,
  filter_error = 0,
  plot_fishtails = NULL,
  adaptive_p_val_bound = T, ...) {

  mymessage(msg = sprintf('Preparing %s', tumor_type))
  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  }

  dtf <- compile_all_coef_overview(
    redo = F,
    ncores = ncores,
    analysis_idxs = analysis_idxs,
    reg_method = reg_method,
    include_non_ds_tallies = include_non_ds_tallies,
    z_normalize = z_normalize,
  ) %>% setDT()
  dtf <- dtf[analysis_idx %in% analysis_idxs]
  print_overview_stats(dtf, 'before any filtering')

  ## As this is an 'analysis' object, a subset of rows will be
  ## specifically for the pan-can setting. For 'data' objects, all
  ## rows are required when the tumor type is 'pan-cancer'
  dtf <- subset_project(dtf = dtf, project_extended = tumor_type,
    pan_can_includes_all = FALSE)
  if (null_dat(dtf)) { return(NULL) }

  if (integration_method == 'union' &&
      'hla_sim_range' %in% colnames(dtf)) {
    dtf <- dtf[hla_sim_range != '<= 1']
  }
  ## Filter out currently undesired analyses
  if (!include_rooney) { dtf <- dtf[!grepl('rooney', analysis_name)] }

  single_hla_mode <- !is.null(hla_allele) && hla_allele != 'integrated'
  if (single_hla_mode) {
    dtf <- dtf[focus_allele == hla_allele]
  }

  print_overview_stats(dtf, 'after analysis type-filtering')

  if (!grepl('bayesian', reg_method)) {
    if (filter_intercept == 'positive') {
      ## An intercept of 0 indicates an underlying anomaly if
      ## Z-normalization == TRUE; we cannot have less than 0
      ## neo-antigens per mutation in the real
      ## world
      dtf <- dtf[intercept > 0]
    } else if (filter_intercept == 'negative') {
      dtf <- dtf[intercept < 0]
    }
    if (null_dat(dtf)) { return(NULL) }
    dtf[, 'yr_fractional_change' := rc / intercept]

    if (!single_hla_mode) {
      message('Picking representative (central) alleles')
      ## Summarize the five focus alleles down to one representative
      ## analysis
      dtf <- pick_representative_allele(dtf, sort_var = sort_var,
        method = pick_allele_method)
      print_overview_stats(dtf, 'after picking representative allele')
    }
    p_var_name <- determine_p_var(fill_var)
    adj_p_var_name <- sprintf('%s.adj', p_var_name)
    signif_effect_name <- sprintf('%s_signif_effect', p_var_name)
    adj_signif_effect_name <-
      sprintf('%s_signif_effect.adj', p_var_name)
    fv_threshold <- determine_threshold(fill_var)

    ## Compute FDR-corrected p-values
    setDT(dtf)
    dtf[, (adj_p_var_name) := p.adjust(p_val, method = 'BY')]
    stopifnot(c(p_var_name, adj_p_var_name, fill_var) %in% colnames(dtf))
    stopifnot(all(dtf[, get(adj_p_var_name) >= get(p_var_name)], na.rm = T))
    dtf[, (signif_effect_name) :=
      ifelse(get(p_var_name) <= p_val_bound,
        ifelse(get(fill_var) <= fv_threshold,
          'L', 'G'), '')]
    dtf[, (adj_signif_effect_name) :=
      ifelse(get(adj_p_var_name) <= p_val_bound,
        ifelse(get(fill_var) <= fv_threshold,
          'L', 'G'), '')]
    ## Report on the number of 'significant' findings
    dplyr::count(dtf, .data[[signif_effect_name]]) %>%
      print
    dplyr::count(dtf, .data[[adj_signif_effect_name]]) %>%
      print
  } else {
    dtf <- filter_bayesian_coef_overview(
      dtf = dtf,
      pick_allele_method = pick_allele_method,
      plot_fishtails = plot_fishtails,
      filter_error = filter_error
    )
  }
  if (null_dat(dtf)) { return(NULL) }

  return(dtf)
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


perform_filter_coef_overview <- function(dtf, code, 
  print_messages = T) {

  code <- rlang::enquo(code)
  comp_before <- assess_coef_composition(dtf)
  bools <- rlang::eval_tidy(code, data = dtf)
  bools[is.na(bools)] <- F
  ## Overall fold decrease
  N_before <- nrow(dtf)
  overall <- mean(bools) %>% log2
  dtf <- dtf[bools, ]
  N_after <- nrow(dtf)
  comp_after <- assess_coef_composition(dtf)

  if (!print_messages) return(dtf)

  stats <- map_dfr(names(comp_before), function(n) {
      out <- merge(comp_before[[n]], comp_after[[n]], by = n)
      out$log2_fc <- log2(out[, 3]) - log2(out[, 2])
      out$diff_evenness <- compute_evenness(-out$log2_fc)
      out$level <- out[[n]]
      out[[n]] <- NULL
      out <- tibble('variable' = n, out)
      return(out)
    }) %>%
    dplyr::mutate(variable = as.character(variable)) %>%
    rbind(
      tibble(
        variable = 'overall',
        level = NA,
        'diff_evenness' = NA,
        'N.x' = N_before, 'N.y' = N_after,
        log2_fc = overall
      )) %>%
    dplyr::mutate(fc = 2^log2_fc)

  # return(list('dtf' = dtf, 'stats' = stats))
  attr(dtf, 'latest_stats') <- stats
  return(dtf)
}


filter_coef_overview <- function(
  dtf,
  print_messages = T,
  force_positive_intercept = F,
  intercept_filter_p_val = NULL,
  intercept_filter_magnitude = 1e-16,
  rc_filter_magnitude = NULL,
  scale_filter = NULL,
  AFDP_50_filter = NULL,
  AFDP_75_filter = NULL,
  AFDP_90_filter = NULL,
  adaptive_scale_filter = FALSE,
  min_patients = NULL,
  min_project_size = NULL) {

  tmp <- dtf[
    maartenutils::eps(scale, 0, 1e-16) & 
    intercept > 1e-3]
  if (nrow(tmp) > 0) stop('Implement zero error filtering!')

  dtf <- perform_filter_coef_overview(
    dtf = dtf,
    print_messages = print_messages,
    code = !is.na(rc) & !is.na(intercept) & converged
  )
  if (print_messages)
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'Filtering non-converged sub-analyses'
    )

  if (!is.null(intercept_filter_p_val)) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = p_val_intercept <= intercept_filter_p_val
    )
    if (print_messages)
      print_overview_stats(
        dtf = dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After intercept p-val filtering'
      )
  }

  if (!is.null(intercept_filter_magnitude)) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(intercept) >= intercept_filter_magnitude
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After intercept magnitude filtering'
      )
  }

  if (force_positive_intercept) {
    dtf <- perform_filter_coef_overview(
      print_messages = print_messages,
      dtf = dtf, code = intercept > 0
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Force positive intercept'
      )
  }

  if (!is.null(scale_filter)) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(norm_scale) <= scale_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After scale filtering'
      )
  }

  if (!is.null(AFDP_50_filter)) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = AFDP_50 <= AFDP_50_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After AFDP_50 filtering'
      )
  }

  if (!is.null(AFDP_75_filter)) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = AFDP_75 <= AFDP_75_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After AFDP_75 filtering'
      )
  }

  if (!is.null(AFDP_90_filter)) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = AFDP_90 <= AFDP_90_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After AFDP_90 filtering'
      )
  }

  if (!is.null(adaptive_scale_filter) && 
      dtf[, any(intercept < 0)]) {
    ## Determine what scale is not tolerable
    thresh <- dtf[intercept < 0, min(abs(norm_scale))]
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(norm_scale) < thresh
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After adaptive scale filtering'
      )
  }

  if (!is.null(rc_filter_magnitude)) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(rc / intercept) <= rc_filter_magnitude
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After relative rc magnitude filtering'
      )
  }

  if (!is.null(min_patients) && 
      dtf[, any(n_patients < min_patients)]) {
    dtf <- perform_filter_coef_overview(
      dtf = dtf,
      print_messages = print_messages,
      code = n_patients >= min_patients
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Filtering minimum number of pts per analysis'
      )
  }

  if (!is.null(min_project_size)) {
    project_counts <- dtf[, .N, by = project_extended]
    allowed_projects <-
      project_counts[N >= min_project_size, project_extended]
    dtf <- dtf[project_extended %in% allowed_projects]
  }

  dtf[, tumor_type := tumor_types[as.character(project_extended)]]
  ## Order projects by mean yr_fractional_change
  if (F) {
    pe <- dtf[,
      wa_var(.SD, varname = 'yr_fractional_change',
        error_var = 'norm_scale'),
      by = .(tumor_type, project_extended)] %>%
      .[order(V1), .(project_extended, tumor_type, V1)]
  } else {
    pe <- dtf[, mean(yr_fractional_change),
      by = .(tumor_type, project_extended)] %>%
      .[order(V1), .(project_extended, tumor_type, V1)]
  }
  if (print_messages)
    print(pe)
  dtf[, project_extended := factor(project_extended,
    levels = pe$project_extended)]
  dtf[, tumor_type := factor(tumor_type, levels = pe$tumor_type)]

  return(dtf)
}


filter_bayesian_coef_overview <- function(
  dtf,
  pick_allele_method = 'est_error',
  intercept_filter_error = NULL,
  rc_filter_error = NULL,
  delta_filter_error = NULL,
  single_hla_mode = data.table::uniqueN(dtf$project_extended) == 1,
  filter_intercept = T,
  plot_fishtails = NULL) {

  print_overview_stats(dtf,
    plot_fishtails = plot_fishtails,
    stage_id = 'Before any filtering'
  )

  ## Filter away low-confidence analyses
  if (dtf[, any(!sampling_convergence)]) {
    dtf <- dtf[sampling_convergence == T]
    print_overview_stats(
      dtf = dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After filtering for sampling convergence'
    )
  }

  if (!is.null(intercept_filter_error)) {
    dtf <- dtf[-log2_intercept_cov >= intercept_filter_error]
    print_overview_stats(
      dtf = dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After filtering for intercept error'
    )
  }

  if (!is.null(rc_filter_error)) {
    dtf <- dtf[-log2_rc_error >= rc_filter_error]
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After filtering for rc error')
  }

  if (!is.null(delta_filter_error)) {
    dtf <- dtf[-log2_est_error >= delta_filter_error]
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After filtering for delta error')
  }

  if (!single_hla_mode && pick_allele_method == 'est_error') {
    s_dtf <- copy(dtf[order(log2_est_error), .SD[1],
      by = analysis_grp_vars])
    dtf <- s_dtf
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After picking representative alleles')
  }

  if ('intercept' %in% colnames(dtf) &&
      filter_intercept && dtf[intercept < 0, .N] > 0) {
    dtf <- dtf[intercept >= 0]
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After removing negative intercept sub-analyses')
  }
  return(dtf)
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
#' Sort the alleles with respect to the \code{sort_var} and slice out the middle
#' one. The different methods pertain to different ways of dealing with cases
#' where there are an even number of alleles to pick a representative allele
#' from.
#'
#' aggressive: pick the lowest (most consistent with immunoediting/neo-antigen
#' depletion)
#' careful: select middle two analyses and average their results
#' conservative: pick the highest (most consistent with neo-antigen enrichment)
#'
pick_representative_allele <- function(dtf, sort_var = 'depletion_full',
  method = 'DT_careful') {

  dtf <- format_overview_res(dtf)

  if (!sort_var %in% colnames(dtf)) {
    mystop(msg = glue('{sort_var} not in input'))
  }

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
    group_vars <- c('overlap_var', 'patient_inclusion_crit', 'LOH_HLA',
      'analysis_name', 'analysis_idx', 'project_extended')
    rc_sel <- dtf %>%
      .[!is.na(get(sort_var))] %>%
      .[order(get(sort_var)), cbind(.SD, N_alleles = .N), by = group_vars] %>%
      .[, i := 1:.N, by = group_vars]

    if ('opt_string' %in% colnames(rc_sel)) {
      ## Strip allele away from opt string, since we're integrating it out
      rc_sel[, opt_string := gsub('[^ ]*\\s{1}(.*)', '\\1', opt_string)]
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
        by = group_vars, .SDcols = c(numeric_columns)]
      df_b <- rc_sel_even[, lapply(.SD, hmp_wrapper),
        by = group_vars, .SDcols = p_val_columns]
      rc_sel_even <- merge(df_a, df_b, by = group_vars)
    }

    if (nrow(rc_sel_uneven) > 0) {
      rc_sel_uneven <-
        rc_sel_uneven[i == ceiling(N_alleles / 2), .SD, by = group_vars]
    }

    rc_sel <- rbind(rc_sel_even, rc_sel_uneven, fill = T)
    rc_sel[, i := NULL]
    rc_sel[, focus_allele := NULL]
    na_cols <- map_lgl(rc_sel, ~all(is.na(.x))) %>% which %>% names
    if (length(na_cols) > 0) {
      rc_sel[, !na_cols]
    }
  } else if (method == 'DT_aggressive') {
    setDT(dtf)
    analysis_id_cols <- colnames(all_cont_IE_settings) %>%
      setdiff('focus_allele') %>%
      c('analysis_idx', 'project_extended')
    rc_sel <-
      dtf[, cbind(.SD[order(rc)][max(ceiling(.N / 2), 1)], 'N_alleles' = .N),
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
  res <- dplyr::group_by(dtf, overlap_var, patient_inclusion_crit, LOH_HLA,
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
