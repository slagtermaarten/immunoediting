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
  minimum_mod_time = '2022-11-03 11:00', verbose = TRUE)
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
    'n_nz_patients' = dtf[c > 0, .N]
  )
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
    f <- suppressMessages(pryr::partial(fit_rlm_model,
        remove_low_TMB_patients = F))
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
  } else if (reg_method == 'glm') {
    f <- glm_fit_model
  } else if (reg_method == 'glm_log') {
    f <- pryr::partial(glm_fit_model, log_transform = T)
  } else if (reg_method == 'glm_log_CYT') {
    f <- pryr::partial(glm_fit_model, log_transform = TRUE,
      fit_CYT = TRUE)
  } else if (reg_method == 'lmp') {
    f <- lmp_fit_model
  } else {
    f <- get(reg_method)
  }
  lm_results <- purrr::map(p_dtf, f)

  if (TRUE && !grepl('wilcox', reg_method)) {
    ## Append Wilcoxon-test results if that's not all that was asked
    ## for
    wilcox_results <- purrr::map(p_dtf, fit_wilcox_model)
    out <- purrr::map(seq(p_dtf),
      ~c(lm_results[[.x]], wilcox_results[[.x]])
    )
  }

  if (TRUE) {
    ## Append general stats
    results <- purrr::map(p_dtf, compute_dtf_stats)
    out <- purrr::map(seq(p_dtf),
      ~c(lm_results[[.x]], results[[.x]])
    )
  }

  return(out)
}


fit_wilcox_model <- function(dtf) {
  if (is.null(dtf)) return(NULL)
  setDT(dtf)[, 'ol_quartile' := cut(ol, breaks = seq(0, 1, by = .25))]
  dtf <- dtf[as.integer(ol_quartile) %in% c(1, 4)]

  wc <- tryCatch(
    suppressWarnings(wilcox.test(y_var ~ ol_quartile,
        data = dtf, conf.int = T)),
    error = function(e) { NULL })
  if (is.null(wc)) return(NULL)

  return(list(
    'wilcox_estimate' = unname(wc$estimate),
    'wilcox_ci_l' = wc$conf.int[1],
    'wilcox_ci_h' = wc$conf.int[2],
    'wilcox_p' = wc$p.value,
    'wilcox_n_patients' = dtf[, .N],
    'wilcox_n_patients_l' = dtf[as.integer(ol_quartile) == 1, .N],
    'wilcox_n_patients_h' = dtf[as.integer(ol_quartile) == 4, .N],
    'wilcox_stat' = unname(wc$statistic),
    'wilcox_method' = wc$method
  ))
}


fit_weighted_t_test  <- function(dtf) {
  if (is.null(dtf)) return(NULL)
  setDT(dtf)[, 'ol_quartile' :=  cut(ol, breaks = seq(0, 1, by = .25))]
  dtf <- dtf[as.integer(ol_quartile) %in% c(1, 4)]

  test <- tryCatch(weights::wtd.t.test(
    x = dtf[as.integer(ol_quartile) == 1, y_var],
    y = dtf[as.integer(ol_quartile) == 4, y_var],
    weight = dtf[as.integer(ol_quartile) == 1, c],
    weighty = dtf[as.integer(ol_quartile) == 4, c],
    samedata = FALSE
  ), error = function(e) { NULL })

  if (is.null(test)) return(NULL)

  out <- unlist(test[c('coefficients', 'additional')])
  names(out) <- names(out) %>%
    str_replace_all('\\.', '_') %>%
    str_replace_all(' ', '') %>%
    tolower %>%
    { paste('wt_', ., sep = '') }

  return(out)
}


rlm_high_TMB_patients <- function(...) {
  # dots <- as.list(...)
  # fit_rlm_model(..., remove_low_TMB_patients = F)
  fit_rlm_model(..., debug = T, remove_low_TMB_patients = T)
}


#' Fit robust regression model with Huber error model
#'
#'
fit_rlm_model <- function(
  dtf, perform_sanity_checks = F,
  weighted = T, include_lm = F, remove_low_TMB_patients = T,
  debug = F, N_perms = 100) {

  if (null_dat(dtf)) return(NULL)
  start_time <- Sys.time()
  if (!test_data_sufficiency(dtf))
    return(tibble_row(message = 'data_insufficient'))

  if (remove_low_TMB_patients) {
    l_dtf <- filter_low_TMB_patients(dtf)
    nz_report <- attr(l_dtf, 'nz_report')
    if (!test_data_sufficiency(l_dtf)) {
      return(
        tibble_row(message = 'data_insufficient_after_NZ_filtering'))
    }
  } else {
    l_dtf <- dtf
    nz_report <- list()
  }

  unnormalized_mod <- fit_rlm_(
    dtf = l_dtf,
    weighted = weighted,
    simple_output = T
  )
  if (is.null(unnormalized_mod))
    return(tibble_row(message = 'rlm_model_fit_failed'))

  orig_pars <- coef(unnormalized_mod) %>%
    as.list %>% setNames(c('intercept', 'rc'))

  pre_trans_scale_avg_norm <-
    unnormalized_mod$s / median(dtf$y_var, na.rm = T)

  l_dtf$y_var <- l_dtf$y_var / orig_pars[['intercept']]
  normalized_mod <- fit_rlm_(l_dtf, weighted = weighted,
    simple_output = F, include_lm = F)

  if (is.null(normalized_mod))
    return(tibble_row(message = 'normalized_rlm_model_fit_failed'))

  if (!is.null(N_perms) && N_perms > 0) {
    perm_res <- map_dfr(1:N_perms, function(n) {
      ord <- sample(1:nrow(l_dtf))
      p_dtf <- l_dtf[, .(ol = ol[ord], c, i, y_var)]
      p_unnormalized_mod <- fit_rlm_(p_dtf, weighted = weighted,
        simple_output = T)
      if (is.null(p_unnormalized_mod)) return(NULL)
      p_dtf$y_var <- p_dtf$y_var / coef(p_unnormalized_mod)[1]
      fit_rlm_(p_dtf, include_lm = F, weighted = weighted)
    })

    num_cols <- map_lgl(perm_res, is.numeric) %>%
      which %>% names %>% { . }

    ## Compute stats over all numerical values of the permutation
    ## distribution
    perm_stats <- perm_res %>%
      dplyr::select(any_of(num_cols)) %>%
      summarise(across(everything(),
          c('median' = median, 'mad' = mad, 'mean' = mean),
          na.rm = T)) %>%
      { . }

    # perm_stats[['delta_mean_median']]
    delta_mean_c_mean <-
      normalized_mod$delta_mean - perm_stats[['delta_mean_mean']]
    perm_var <- 1/N_perms * sum(perm_res$delta_SE^2)
    ## SE of 'corrected' distribution is sum of original and
    ## permutation dist
    delta_mean_c_SE <- sqrt(normalized_mod$delta_SE^2 + perm_var)

    ## _pq is observation quantile: fraction of permutation
    ## observations below the observed value for given statistic
    perm_stats %<>%
      append(
        qnorm(p = c(.025, .5, .975),
          mean = delta_mean_c_mean, sd = delta_mean_c_SE) %>%
        set_names(c('delta_mean_c_ci_l', 'delta_mean_c_mean',
            'delta_mean_c_ci_h')) %>%
        as.list()
      )

    for (cn in num_cols) {
      acn <- glue('{cn}_pq')
      perm_stats[[acn]] <-
        mean(normalized_mod[[cn]] >= perm_res[[cn]])
    }

    if (F) {
      for (cn in c('delta_mean', 'delta_CI_l', 'delta_CI_h')) {
        if (F) {
          ccn <- glue('{cn}_pc')
          perm_stats[[ccn]] <-
            normalized_mod[[cn]] - median(perm_res[[cn]])
        } else {
          ccn <- glue('{cn}_pc') %>% paste0(c('', '_ci_l', '_ci_h'))
          norm_dist <- normalized_mod[[cn]] - perm_res[[cn]]
          # hist(norm_dist, breaks = 100)
          perm_stats %<>% append(
            DescTools::MeanCI(norm_dist) %>%
              setNames(ccn) %>%
              as.list
          )
        }
      }
    }

    if (F) {
      # hist(perm_res$delta_mean, breaks = 100)
      delta_mean_sum <- summary(perm_res$delta_mean) %>%
        set_names(gsub('\\.| ', '', tolower(names(.))))
      d1 <- abs(delta_mean_sum['median'] - delta_mean_sum['1stqu'])
      d2 <- abs(delta_mean_sum['median'] - delta_mean_sum['3rdqu'])
      ## The scale in which the median operates
      med_scale <- floor(log10(abs(delta_mean_sum['median'])))
      ## Test whether
      perm_stats[['symmetry_test']] <-
        maartenutils::eps(d1, d2, 10^med_scale)
    }

    perm_stats[c('delta_mean_median', 'delta_mean_c_mean',
      'delta_mean_pc')]

    perm_stats <- perm_stats %>%
      set_names(paste('perm', names(.), sep = '_')) %>%
      { . }
    # stopifnot(!all(duplicated(names(perm_stats))))

    if (debug) print(perm_stats, width = 3000)
  } else {
    perm_stats <- list()
  }

  if (perform_sanity_checks) {
    trans_test <-
      (abs(pre_trans_scale_avg_norm) == Inf &&
      abs(normalized_mod$scale_avg_norm) == Inf) ||
      maartenutils::eps(
        pre_trans_scale_avg_norm,
        normalized_mod$scale_avg_norm, 1e-3)
    browser(expr = !trans_test)
    stopifnot(trans_test)

    trans_test <- maartenutils::eps(
      orig_pars$rc / orig_pars$intercept,
      normalized_mod$rc, 1e-3)
    browser(expr = !trans_test)
    stopifnot(trans_test)
  }

  out <- normalized_mod %>%
    modifyList(list(message = 'OK')) %>%
    modifyList(orig_pars) %>%
    ## Ensure scale_avg_norm is not affected by 1/intercept scaling
    ## in any way
    modifyList(list('scale_avg_norm' = pre_trans_scale_avg_norm)) %>%
    modifyList(perm_stats) %>%
    map(unname)

  # if (abs(out$delta_mean) > .1 && debug) {
  #   browser()
  # }
  # h_dist <- ds_ol_stats(l_dtf)
  # print_plot(with(h_dist, plot(ol_b, w_y_var_median)))

  end_time <- Sys.time()
  out$elapsed_time <- format(end_time - start_time, format = '%M')

  return(out)
}


fit_rlm_ <- function(dtf, simple_output = F, include_AFDP = F,
  weighted = T, include_lm = !simple_output) {
  dtf <- setDT(dtf)[is.finite(y_var) & is.finite(ol)]

  if (!test_data_sufficiency(dtf)) return(NULL)

  lc <- tryCatch(suppressWarnings(
      MASS::rlm(
        formula = y_var ~ ol,
        weights = if (weighted) dtf$weight else NULL,
        data = dtf,
        wt.method = 'case',
        x.ret = T, model = T, maxit = 1000)),
    error = function(e) { NULL })

  if (simple_output) {
    return(lc)
  }

  quants <- c(.1, .25, .5, .75, .9)
  delta_names <-
    paste0('delta_', c('CI_l', 'mean', 'CI_h'), sep = '') %>%
    c('delta_SE')
  NMADR_eval_locs <- quants %>%
    { setNames(., paste('NMADR_q', 100 * ., sep = '')) }

  if (is.null(lc)) {
    # delta_def <- map(auto_name(delta_names), ~NA_real_)
    # NMADR_def <- map(NMADR_eval_locs, ~NA_real_)
    # out <- list(
    #   'message' = 'model_fit_failed',
    #   'intercept' = NA_real_,
    #   'p_val_intercept' = NA_real_,
    #   'rc' = NA_real_,
    #   'p_val' = NA_real_,
    #   'df' = NA_real_,
    #   't_val' = NA_real_,
    #   'n_patients' = dtf[, .N],
    #   'converged' = F,
    #   'yr_fractional_change' = NA_real_,
    #   'scale' = NA_real_,
    #   'scale_avg_norm' = NA_real_,
    #   'norm_scale' = NA_real_
    # ) %>% c(delta_def) %>% c(NMADR_def)

    # if (include_AFDP) {
    #   AFDP_def <- map(quants, ~NA_real_) %>%
    #     setNames(paste0('AFDP_', 100*quants, sep = ''))
    #   out <- c(out, AFDP_def)
    # }
    # return(out)
    return(list('message' = 'model_fit_failed'))
  }

  coef_mat <- summary(lc)$coefficients
  coef_mat <- cbind(coef_mat,
    'p_value' = sapply(seq_along(coef_mat[, 3]), function(i) {
      tryCatch(sfsmisc::f.robftest(lc, var = i)$p.value,
        error = function(e) { NA })
    }),
    'df' = summary(lc)$df[1:nrow(coef_mat)]
  )

  if (is.null(coef_mat[2,1]) || is.na(coef_mat[2,1])) {
    browser()
  }
  tryCatch(1 / coef_mat[2, 1], error = function(e) { browser() })

  delta <- tryCatch(
    predict(lc, newdata = data.frame(ol = 1), se.fit = T),
    error = function(e) { NULL })
  if (is.null(delta)) {
    return(list(message = 'predict_step_failed'))
  }

  delta <-
    { delta$fit + c(-1.96, 0, 1.96) * delta$se.fit - 1 } %>%
    as.list() %>%
    append(list(delta$se.fit)) %>%
    setNames(delta_names)

  ## Median of absolute residuals after RLM convergence, identical to
  ## lc$s
  # MADR <- median(abs(residuals(lc))) / 0.6745
  intercept <- predict(lc,
    newdata = data.frame(ol = 0), se.fit = T)
  ## Evaluate F_{NMADR}^{-1}(q)
  NMADR_probs <- qnorm(NMADR_eval_locs,
    mean = lc$s / intercept$fit,
    sd = intercept$se.fit) %>%
    as.list
  norm_scale <- lc$s / coef_mat[1, 1]
  stopifnot(maartenutils::eps(
      norm_scale, NMADR_probs[['NMADR_q50']], 1e-3))

  out <- list(
    'message' = 'OK',
    'intercept' = coef_mat[1, 1],
    'p_val_intercept' = coef_mat[1, 4],
    'rc' = coef_mat[2, 1],
    'p_val' = coef_mat[2, 4],
    'df' = coef_mat[2, 5],
    't_val' = coef_mat[2, 3],
    'n_patients' = dtf[, .N],
    'converged' = lc$converged,
    'scale' = lc$s,
    'scale_avg_norm' = lc$s / median(dtf$y_var),
    'norm_scale' = norm_scale
  ) %>% c(delta) %>% c(NMADR_probs)

  if (include_lm)
    out <- append(list('lm' = lc))

  if (include_AFDP) {
    ## Compute absolute fractional difference between predictions
    ## (AFDP)
    ## and observations
    non_zero_idx <- dtf[, y_var > 0 | y_var < 0]
    obs <- dtf[, y_var][non_zero_idx]
    pred <- predict(lc)[non_zero_idx]
    dev <- abs((pred - obs) / obs)
    AFDP_stats <- as.list(quantile(dev, quants)) %>%
      setNames(paste0('AFDP_', 100*quants, sep = ''))
    out <- c(out, AFDP_stats)
  }

  return(out)
}


fit_glm_ <- function(dtf, simple_output = FALSE,
  fit_PS = TRUE, 
  fit_offset = TRUE,
  fit_CYT = F, 
  add_resid_ac = FALSE, 
  diagnose = FALSE,
  do_test_data_sufficiency = TRUE) {

  dtf <- setDT(dtf)[is.finite(i) & is.finite(c) & is.finite(ol)]
  if (do_test_data_sufficiency && !test_data_sufficiency(dtf)) 
    return(NULL)

  # if (method == 'glm') {
  # } else {
  #   fit_func <- pryr::partial(lmPerm::lmp, scale = F)
  # }
  fit_func <- stats::glm

  if (fit_PS) {
    if (fit_CYT) {
      if (fit_offset) {
        fit <- tryCatch(fit_func(c~1+i+i:ol+CYT, data = dtf),
          error = function(e) { NULL })
      } else {
        fit <- tryCatch(fit_func(c~0+i+i:ol+CYT, data = dtf),
          error = function(e) { NULL })
      }
    } else {
      if (fit_offset) {
        fit <- tryCatch(fit_func(c~1+i+i:ol, data = dtf),
          error = function(e) { NULL })
      } else {
        fit <- tryCatch(fit_func(c~0+i+i:ol, data = dtf),
          error = function(e) { NULL })
      }
    }
  } else {
    if (fit_CYT) {
      if (fit_offset) {
        fit <- tryCatch(fit_func(c~1+i+CYT, data = dtf),
          error = function(e) { NULL })
      } else {
        fit <- tryCatch(fit_func(c~0+i+CYT, data = dtf),
          error = function(e) { NULL })
      }
    } else {
      if (fit_offset) {
        fit <- tryCatch(fit_func(c~1+i, data = dtf),
          error = function(e) { NULL })
      } else {
        fit <- tryCatch(fit_func(c~0+i, data = dtf),
          error = function(e) { NULL })
      }
    }
  }

  if (is.null(fit)) return(NULL)

  if (simple_output) {
    return(fit)
  }

  if (is.null(fit)) {
    return(list('message' = 'model_fit_failed'))
  }

  coef_mat <- summary(fit)$coefficients

  out <- list(
    'message' = 'OK',
    'n_patients' = dtf[, .N],
    ## High should mean good correspondence to linear model
    'deviance_red' = 1 - fit$deviance / fit$null.deviance,
    'converged' = fit$converged,
    'median_TMB' = median(fit$data$i)
  ) %>% { . }

  ## This is rather fugly, but pragmatic..?
  p_index = 4
  t_index = 3
  rwNULL <- function(x) tryCatch(x, error = function(e) { NULL })
  if (fit_PS) {
    if (fit_CYT) {
      if (fit_offset) {
        out <- c(out, list(
          'offset' = coef_mat[1, 1],
          'intercept' = coef_mat[2, 1],
          'p_val_intercept' = coef_mat[2, p_index],
          'rc' = coef_mat[3, 1],
          'rc_CYT' = coef_mat[4, 1],
          ## Should be in range [-1, 0]
          'delta' = coef_mat[3, 1] / coef_mat[2, 1],
          'p_val' = coef_mat[3, p_index],
          'p_val_CYT' = coef_mat[4, p_index],
          't_val' = rwNULL(coef_mat[3, t_index]))
        )
      } else {
        out <- c(out, list(
          'intercept' = coef_mat[1, 1],
          'p_val_intercept' = coef_mat[1, p_index],
          'rc' = coef_mat[2, 1],
          'rc_CYT' = coef_mat[3, 1],
          ## Should be in range [-1, 0]
          'delta' = coef_mat[2, 1] / coef_mat[1, 1],
          'p_val' = coef_mat[2, p_index],
          'p_val_CYT' = coef_mat[3, p_index],
          't_val' = rwNULL(coef_mat[2, t_index]))
        )
      }
    } else {
      if (fit_offset) {
        out <- c(out, list(
          'offset' = coef_mat[1, 1],
          'intercept' = coef_mat[2, 1],
          'p_val_intercept' = coef_mat[2, p_index],
          'rc' = coef_mat[3, 1],
          ## Should be in range [-1, 0]
          'delta' = coef_mat[3, 1] / coef_mat[2, 1],
          'p_val' = coef_mat[3, p_index],
          't_val' = rwNULL(coef_mat[3, t_index]))
        )
      } else {
        out <- c(out, list(
          'intercept' = coef_mat[1, 1],
          'p_val_intercept' = coef_mat[1, p_index],
          'rc' = coef_mat[2, 1],
          ## Should be in range [-1, 0]
          'delta' = coef_mat[2, 1] / coef_mat[1, 1],
          'p_val' = coef_mat[2, p_index],
          't_val' = rwNULL(coef_mat[2, t_index]))
        )
      }
    }
  } else {
    if (fit_CYT) {
      if (fit_offset) {
        out <- c(out, list(
          'offset' = coef_mat[1, 1],
          'intercept' = coef_mat[2, 1],
          'p_val_intercept' = coef_mat[2, p_index],
          'rc_CYT' = coef_mat[3, 1],
          'p_val_CYT' = coef_mat[3, p_index]
        ))
      } else {
        out <- c(out, list(
          'intercept' = coef_mat[1, 1],
          'p_val_intercept' = coef_mat[1, p_index],
          'rc_CYT' = coef_mat[2, 1],
          'p_val_CYT' = coef_mat[2, p_index]
        ))
      }
    } else {
      if (fit_offset) {
        out <- c(out, list(
          'offset' = coef_mat[1, 1],
          'intercept' = coef_mat[2, 1],
          'p_val_intercept' = coef_mat[2, p_index]
        ))
      } else {
        out <- c(out, list(
          'intercept' = coef_mat[1, 1],
          'p_val_intercept' = coef_mat[1, p_index]
        ))
      }
    }
  }

  if (add_resid_ac) {
    acf_res <- acf(residuals(fit), pl=FALSE, lag.max = 5)
    for (idx in 1:5) {
      ## Plus 1 to idx as the acf result is zero-indexed and the zero
      ## index is always 1, i.e. uninformative
      out[[glue::glue('resid_ac_{idx}')]] <- acf_res$acf[, , 1][idx+1]
    }

    dw_test <- car::durbinWatsonTest(fit, max.lag=1) %>%
      as.list() %>%
      { setNames(., paste0('dw_', names(.))) } %>%
      { . }
    out <- c(out, dw_test)
  }

  out$baseline_yr <- coef(fit)[['(Intercept)']] + 
        mean(fit$data$i) * coef(fit)[['i']]

  if (diagnose) {
    leverage <- hatvalues(fit)
    cooksd <- cooks.distance(fit)
    out$leverage_evenness <- maartenutils::compute_evenness(abs(leverage))
    out$cooksd_evenness <- maartenutils::compute_evenness(abs(cooksd))

    ps_fit <- lm(i ~ ol, data = dtf, weights = leverage)
    # if (coef(ps_fit)['(Intercept)'] < 0) browser()
    out$ol_vs_i_rc <- coef(ps_fit)['ol']
    out$ol_vs_i_rc_n <- coef(ps_fit)['ol'] / coef(ps_fit)['(Intercept)']
  }

  if (fit_offset) {
    out$baseline_NAYR <- coef(fit)[['(Intercept)']] / out$median_TMB + 
      coef(fit)[['i']]
    if (fit_PS) {
      out$slope_NAYR <- 1/out$median_TMB * coef(fit)[['i:ol']]
    }
  }

  return(out)
}


fit_glm_log_ <- function(dtf, simple_output = FALSE,
  add_resid_ac = FALSE,
  fit_offset = TRUE, fit_PS = TRUE, fit_CYT = FALSE,
  do_test_data_sufficiency = TRUE) {

  if (maartenutils::null_dat(dtf)) return(NULL)
  if (is.null(dtf$c) || all(is.na(dtf$c))) return(NULL)
  if (is.null(dtf$i) || all(is.na(dtf$i))) return(NULL)

  dtf %>%
    dplyr::mutate(c = log10(c + 1)) %>%
    dplyr::mutate(i = log10(i + 1)) %>%
    fit_glm_(fit_offset = fit_offset, fit_PS = fit_PS,
      add_resid_ac = add_resid_ac,
      fit_CYT = fit_CYT, 
      do_test_data_sufficiency = do_test_data_sufficiency,
      simple_output = simple_output)
}


gen_sample_data <- function(N_p = 500, sd = 1) {
  dtf <- tibble(
    i = rpois(N_p, 20),
    c = rnorm(N_p, i * .15, sd),
    ol = runif(N_p)
  )
  return(dtf)
}
# print(fit_glm_(gen_sample_data(sd = 10)))
# print(fit_glm_(gen_sample_data()))


fit_nb_mod = function(x, formula) {
  tryCatch(suppressWarnings(
    MASS::glm.nb(formula, data = x,
    control = glm.control(maxit = 1000))
  ), error = function(e) { NULL })
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


fit_dropout_aware_nb <- function(
  dtf,
  remove_low_TMB_patients = T,
  remove_outliers = F,
  N_draws = 1e3,
  draw_res = 1e2,
  verbose = F) {

  time_start <- Sys.time()

  if (remove_low_TMB_patients) {
    l_dtf <- filter_low_TMB_patients(dtf)
    nz_report <- attr(l_dtf, 'nz_report')

    if (is.null(l_dtf)) return(tibble_row(message = 'NZ_failed'))
    if (!test_data_sufficiency(l_dtf)) {
      return(
        tibble_row(message = 'data_insufficient_after_NZ_filtering')
      )
    }
  } else {
    l_dtf <- dtf
    nz_report <- list()
  }

  ## PHASE 1.5 remove outlying patients. Not all that necessary
  ## because of the IWLS in glm.nb
  # test_plot({ par(mfrow = c(2, 2)); plot(mod) }, w = 17.4, h = 20)
  if (remove_outliers) {
    start_n <- nrow(l_dtf)
    prev_n <- NULL
    i <- 0
    while ((is.null(prev_n) || prev_n != nrow(l_dtf)) && i <= 10) {
      prev_n <- nrow(l_dtf)
      m3 <- fit_nb_mod(l_dtf, c ~ 1 + i + ol + i:ol)
      if (is.null(m3)) return(tibble_row(
          message = 'remove_outlying_patients_failed'))
      l_dtf <- l_dtf[abs(studres(m3)) < 3, ]
      i <- i + 1
    }
    outlying_pts <- start_n - prev_n
    # cooks <- NULL
    # i <- 0
    # while (nrow(l_dtf) > 0 && is.null(cooks)) {
    #   i <- i + 1
    #   cooks <- cooks.distance(mod)
    #   if (!is.null(cooks)) {
    #     cooks_thresh <- quantile(cooks, .975)
    #     allowed <- which(cooks.distance(mod) < cook_thresh)
    #     l_dtf <- l_dtf[allowed, ]
    #   }
    #   mod <- fit_nb_mod(l_dtf)
    #   # test_plot(hist(cooks, breaks = 100))
    # }
    # if (!test_data_sufficiency(l_dtf)) {
    #   return(NULL)
    # }
  } else {
    i <- 0
    outlying_pts <- 0
  }

  ## PHASE 2: test correlation between i and i:ol, strong enough for
  ## inclusion?
  mods <-
    list(
      c ~ 1 + i + ol,
      # c ~ 1 + i + i:ol,
      c ~ 1 + i + ol + i:ol
    ) %>% map(fit_nb_mod, x = l_dtf)
  mods <- mods[!sapply(mods, is.null)]
  if (length(mods) == 0) return(tibble_row(
      message = 'fit_full_models_failed'))

  if (F) {
    test_plot(plot(l_dtf$c, predict(m3)), w = 20, h = 20)
    # test_plot({ par(mfrow = c(2, 2)); plot(m3) }, w = 17.4, h = 20)
    summary(residuals(mods[[1]]))
    summary(residuals(mods[[2]]))
    summary(residuals(m3))
    cor(studres(mods[[1]]), studres(mods[[2]]))
    cor(studres(mods[[2]]), studres(m3))
    cor(studres(mods[[1]]), studres(m3))
    test_plot(hist(studres(m3)), w = 20, h = 20)
  }

  AICs <- map_dbl(mods, AIC)
  pref_mod <- which.min(AICs)
  ## Compute relative likelihood of models
  ## https://en.wikipedia.org/wiki/Akaike_information_criterion
  rel_lik <- exp((AICs[pref_mod] - AICs) / 2) %>% { . / sum(.) }
  # corM <- cov_to_cor(vcov(mods[[2]]))

  N_mut = median(dtf$i)

  ## PHASE 3: fit model to all sufficiently mutated patients
  ## Transform random draws to quantity of interest, delta, using
  ## random sampling from parameter distribution
  ## Compare h = 1 and h = 0 for a median number of mutations
  ## h = ol in code

  ## Another, exact, option would be 'ray-scanning'
  ## https://nl.mathworks.com/matlabcentral/fileexchange/84973
  ## -integrate-and-classify-normal-distributions

  delta_1 <- mvtnorm::rmvnorm(
      n = min(draw_res^3, N_draws), mean = coef(mods[[1]]),
      sigma = vcov(mods[[1]])
    ) %>% {
      base <- .[,'(Intercept)'] + .[, 'i'] * N_mut; .[, 'ol'] / base
    } %>% {
      c('delta_median' = median(.), 'delta_mad' = mad(.))
    }

  # delta_2 <- mvtnorm::rmvnorm(
  #     n = min(draw_res^3, N_draws),
  #     mean = coef(mods[[2]]), sigma = vcov(mods[[2]])
  #   ) %>% {
  #     base <- .[, '(Intercept)'] + .[, 'i'] * N_mut;
  #     (.[, 'i:ol'] * N_mut) / base
  #   } %>% {
  #     c('delta_median' = median(.), 'delta_mad' = mad(.))
  #   }

  delta_3 <- mvtnorm::rmvnorm(
      n = min(draw_res^4, N_draws), mean = coef(mods[[2]]),
      sigma = vcov(mods[[2]])
    ) %>%
    {
      base <- .[, '(Intercept)'] + .[, 'i'] * N_mut;
      (.[, 'ol'] + .[, 'i:ol'] * N_mut) / base
    } %>%
    { c('delta_median' = median(.), 'delta_mad' = mad(.)) }

  delta <- weighted.mean(c(delta_1[1], delta_3[1]), rel_lik)
  delta_mad <- weighted.mean(c(delta_1[2], delta_3[2]), rel_lik)

  rlm_delta <- fit_rlm_model(l_dtf, include_lm = include_lm)

  out <- tibble_row(
    message = 'OK',
    delta = delta,
    delta_mad = delta_mad,
    mod1_converged = mods[[1]]$converged,
    # mod2_converged = mods[[2]]$converged,
    mod3_converged = mods[[2]]$converged,
    aic1 = AICs[1],
    # aic2 = AICs[2],
    aic3 = AICs[2],
    aic_confidence1 = rel_lik[1],
    # aic_confidence2 = rel_lik[2],
    aic_confidence3 = rel_lik[2],
    cleaning_i = i,
    dropped_pts = nrow(dtf) - nrow(l_dtf),
    low_tmb_pts = nrow(dtf) - nrow(l_dtf) + outlying_pts,
    outlying_pts = outlying_pts,
    theta1 = mods[[1]]$theta,
    # theta2 = mods[[2]]$theta,
    theta3 = mods[[2]]$theta
  ) %>% append(nz_report) %>% append(rlm_delta)

  stopifnot(!any(duplicated(names(out))))

  return(out)
}


fit_lm_model <- function(dtf, include_lm = F) {
  dtf <- dtf[is.finite(y_var) & is.finite(ol)]

  lc <- tryCatch(lm(y_var ~ ol, data = dtf),
    error = function(e) { print(e); NULL })
  coef_mat <- summary(lc)$coefficients
  if (is.null(lc) || dim(coef_mat)[1] == 1) return(NULL)

  res <- list(
    # 'lm' = lc,
    'intercept' = coef_mat[1, 1],
    'p_val_intercept' = coef_mat[1, 4],
    'rc' = tryCatch(coef_mat[2, 1], error = function(e) { NA_real_ }),
    'p_val' = coef_mat[2, 4],
    't_val' = coef_mat[2, 3],
    'n_patients' = dtf[, .N]
  )

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


glm_fit_model <- function(
  dtf,
  remove_low_TMB_patients = FALSE,
  log_transform = FALSE,
  use_lmPerm = FALSE,
  include_perms = FALSE,
  # N_perms = 1000) {
  N_perms = N_perms_global) {

  if (null_dat(dtf)) return(NULL)
  dtf <- dtf[is.finite(y_var) & is.finite(ol)]
  if (null_dat(dtf)) return(NULL)

  if (remove_low_TMB_patients) {
    l_dtf <- filter_low_TMB_patients(dtf)
    nz_report <- attr(l_dtf, 'nz_report')
    if (!test_data_sufficiency(l_dtf)) {
      return(
        tibble_row(message = 'data_insufficient_after_NZ_filtering'))
    }
  } else {
    l_dtf <- dtf
    nz_report <- list()
  }

  if (log_transform) {
    fit <- fit_glm_log_(l_dtf, add_resid_ac = TRUE)
  } else {
    fit <- fit_glm_(l_dtf, add_resid_ac = TRUE)
  }

  if (is.null(fit))
    return(tibble_row(message = 'glm_model_fit_failed'))

  if (!is.null(N_perms) && N_perms > 0) {
    if (!use_lmPerm) {
      perm_res <- purrr::map_dfr(1:N_perms, function(n) {
        ord <- sample(1:nrow(l_dtf))
        p_dtf <- l_dtf[, .(ol = ol[ord], c, i)]
        p_fit <- fit_glm_(p_dtf, simple_output = FALSE)
        return(p_fit)
      })

      num_cols <- map_lgl(perm_res, is.numeric) %>%
        which %>% names %>% { . }

      ## Compute stats over all numerical values of the permutation
      ## distribution
      perm_stats <-
        perm_res %>%
        dplyr::select(any_of(num_cols)) %>%
        summarise(across(everything(),
            c('median' = median, 'mad' = mad, 'mean' = mean),
            na.rm = T)) %>%
        { . }

      # delta_mean_c_mean <-
      #   fit$delta_mean - perm_stats[['delta_mean_mean']]
      # perm_var <- 1/N_perms * sum(perm_res$delta_SE^2)
      # ## SE of 'corrected' distribution is sum of original and
      # ## permutation dist

      ## Compute what fraction of the permutations are more extreme in
      ## absolute sense. Do f(x) 1-x on this to get permutation
      ## p-values
      for (cn in num_cols) {
        acn <- glue('{cn}_apq')
        perm_stats[[acn]] <-
          mean(abs(fit[[cn]]) >= abs(perm_res[[cn]]))
      }

      for (cn in num_cols) {
        acn <- glue('{cn}_pq')
        perm_stats[[acn]] <-
          mean(fit[[cn]] >= perm_res[[cn]])
      }

      perm_stats <- perm_stats %>%
        set_names(paste('perm', names(.), sep = '_')) %>%
        { . }
    } else {
    }
  } else {
    perm_stats <- list()
  }

  out <-
    list(message = 'OK') %>%
    c(fit) %>%
    c(perm_stats)

  if (include_perms) {
    attr(out, 'perms') <- perm_res
  }

  return(out)
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


lmp_fit_model <- function(
  dtf,
  remove_low_TMB_patients = FALSE,
  log_transform = TRUE,
  N_perms = N_perms_global) {

  if (null_dat(dtf)) return(NULL)
  dtf <- dtf[is.finite(y_var) & is.finite(ol)]
  if (null_dat(dtf)) return(NULL)

  if (remove_low_TMB_patients) {
    l_dtf <- filter_low_TMB_patients(dtf)
    nz_report <- attr(l_dtf, 'nz_report')
    if (!test_data_sufficiency(l_dtf)) {
      return(
        tibble_row(message = 'data_insufficient_after_NZ_filtering'))
    }
  } else {
    l_dtf <- dtf
    nz_report <- list()
  }

  if (log_transform) {
    dtf <- dtf %>%
      dplyr::mutate(c = log10(c + 1)) %>%
      dplyr::mutate(i = log10(i + 1)) %>%
      { . }
  }

  fit_func <- pryr::partial(lmPerm::lmp, center = F, scale = F)

  # fit_no_CYT <- tryCatch(fit_func(c~1+i+i:ol, data = dtf),
  #   error = function(e) { NULL })

  fit <- tryCatch(fit_func(c~1+i+i:ol+i:CYT, data = dtf),
    error = function(e) { NULL })

  if (is.null(fit))
    return(tibble_row(message = 'model_fit_failed'))

  out <-
    list(message = 'OK') %>%
    c(broom2row(tidy(fit)))

  out$delta <- out$rc / out$intercept

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


#' Master function for continuous IE detection. Perform and plot IE
#' analyses for all samples combined and for each level of column
#' \code{project_extended} separately
#'
#' @param reg_plot Include regression plots of the combined cohort and
#' all cohorts (levels of \code{project_extended}) separately
#'
test_continuous_IE <- function(
  overlap_var = 'mean_score',
  patient_inclusion_crit = 'none',
  LOH_HLA = 'no_LOHHLA',
  analysis_name = 'twoD_sens_analysis',
  focus_allele = 'A0201',
  analysis_idx = 1,
  reg_method = 'rlm',
  project_extended = NULL,
  hla_sim_range = NULL,
  z_normalize = F,
  redo = F,
  permute_id = NULL,
  return_grob = F,
  ncores = 1,
  check_res = F,
  partition_vars = c(),
  stats_idx = 1,
  skip_missing = F,
  verbose = F,
  permute_input = c(),
  include_call = F,
  ds_frac = NULL,
  return_res = TRUE,
  debug = F) {

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
    permute_input = permute_input,
    LOH_HLA = LOH_HLA,
    patient_inclusion_crit = patient_inclusion_crit,
    overlap_var = overlap_var
  )

  ## Don't cache downsampled/permuted sub-analyses
  if (skip_missing ||
      (is.null(permute_id) && is.null(ds_frac) &&
       !IE_requires_computation(o_fn, verbose = F) &&
       !redo && !check_res)) {
    if (return_res) {
      ret_val <- tryCatch(readRDS(o_fn),
        warning = function(w) {},
        error = function(e) {
          if (!skip_missing) file.remove(o_fn); NULL
        })
      return(ret_val)
    } else {
      return(NULL)
    }
  }

  if (debug) browser()

  prep <- prep_cont_IE_analyses(
    focus_allele = focus_allele,
    redo = redo,
    analysis_name = analysis_name,
    project_extended = project_extended,
    hla_sim_range = hla_sim_range,
    z_normalize = z_normalize,
    LOH_HLA = LOH_HLA,
    permute_y = !is.null(permute_id),
    analysis_idx = analysis_idx,
    overlap_var = overlap_var,
    permute_input = permute_input,
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
        source(file.path(IE_root, 'continuous_IE_detection_init.R'))
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

  reg_plots <- plyr::llply(maartenutils::auto_name(projects),
    function(pe) {
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



#'  Prepare data for multiple types of downstream analyses. Filter
#'  patients and prepare the required regression variables.
#'
#'
prep_cont_IE_analyses <- function() {
  ## Determine donor_summary file to load in
  obj_fn <- gen_cont_IE_ds_fn(
    hla_allele = focus_allele,
    analysis_name = analysis_name,
    analysis_idx = analysis_idx,
    permute_input = permute_input,
    overlap_var = overlap_var
  )
  if (is.null(obj_fn) || is.na(obj_fn) || length(obj_fn) == 0)
    return(NULL)

  ## Load in the data and compute the ordinal variable, save
  ## attributes
  dtf <- readRDS(obj_fn)
  start_N <- dtf[, .N]

  # analysis_name = 'rooney_param_titration'
  if (analysis_name == 'rooney_param_titration') {
    cond_setnames(dtf, 'analyis_name', 'analysis_name')
    dtf <- dtf[grepl('Rooney UMF', analysis_name)]
  }

  attr(dtf, 'y_type') <- if (!grepl('rooney', analysis_name))
    'neo-antigen yield rate' else 'observed/expected neo-antigens'
  attr(dtf, 'x_label') <- sprintf('%s presentation score',
    quickMHC::ppHLA(focus_allele))

  atts <- attributes(dtf)
  dtf <- determine_continuous_IE_yval(dtf,
    replace_zeroes = replace_zeroes)

  ## Extract the  tumor projects present in this object
  dtf <- subset_project(dtf, project_extended)

  ## Downsample patients
  if (!is.null(ds_frac)) {
    dtf <- dtf[sample(1:nrow(dtf), ceiling(nrow(dtf)*ds_frac),
      replace = T), ]
  }

  if (!is.null(permute_y) && permute_y) {
    perm <- sample(1:nrow(dtf))
    dtf$y_var <- dtf$y_var[perm]
    dtf$c <- dtf$c[perm]
    dtf$i <- dtf$i[perm]
  }

  ## Annotate with required additional variables
  if ('IE_essentiality_impaired' %in% partition_vars ||
      any(patient_inclusion_crit != 'none')) {
    dtf <- setDT(dtf)
    FDR_thresh <- switch(patient_inclusion_crit,
      'strict_TR' = .01, 'TR' = .1,
      'FDR10' = .1, 'FDR1' = .01,
      'FDR0.1' = 0.001, 'FDR0.01' = 0.0001,
      'lenient' = 0.0001, 'moderate' = 0.01, 'stringent' = .1
    )
    dtf <-
      annotate_donor_by_IE_essentiality(dtf,
        FDR_thresh = FDR_thresh, verbose = verbose)
    dtf[, mean(IE_essentiality_impaired)]
  }
  if (any(patient_inclusion_crit != 'none')) {
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
  valid_h_N <- dtf[, .N]

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
  valid_y_N <- dtf[, .N]
  if (!is.null(hla_sim_range) && length(hla_sim_range) == 2) {
    dtf <- dtf[ol >= hla_sim_range[1] & ol <= hla_sim_range[2]]
  }
  if (null_dat(dtf) || all(is.na(dtf$y_var)) || all(is.na(dtf$ol))) {
    return(NULL)
  }

  ## Reduce object to crucial columns only in order to make it fit in
  ## memory more easily
  if (reduce_columns) {
    dtf <- reduce_to_essential_IE_columns(dtf)
  }

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

  h_dist <- table(cut(dtf$ol, seq(0, 1, by = .1)))

  out <- list(
    'projects' = observed_projects,
    'dtf' = dtf,
    'h_dist' = h_dist,
    'obj_fn' = obj_fn,
    'patient_counts' = c(
      'start' = start_N,
      'valid_h' = valid_h_N, 'valid_y' = valid_y_N,
      'final' = dtf[, .N])
  )
  return(out)
}
formals(prep_cont_IE_analyses) <- formals(test_continuous_IE)
formals(prep_cont_IE_analyses)$min_points_lq <-
  formals(prep_cont_IE_analyses)$min_points_uq <- 3
formals(prep_cont_IE_analyses)$reduce_columns <- T
formals(prep_cont_IE_analyses)$replace_zeroes <- T
formals(prep_cont_IE_analyses)$permute_y <- F


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
  permute_id = NULL,
  stats_idx = 1,
  return_res = F,
  overlap_var = 'mean_score',
  skip_missing = F,
  plot_ranges = NULL,
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
      permute_id = permute_id,
      analysis_name = analysis_name,
      analysis_idx = analysis_idx,
      skip_missing = skip_missing,
      patient_inclusion_crit = patient_inclusion_crit,
      overlap_var = overlap_var
    )
    if (null_dat(dtf)) return(NULL)
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
    projects <- names(dtf$stats)[good_idx]
    lol <- lol[good_idx]
    if (is.null(lol)) {
      return(NULL)
    }

    ## `test_continuous_IE` can return bare Bayesian model fits or
    ## lists of summary stats already extracted from these fits. Both
    ## scenarios are accounted for here.
    if (length(lol) > 0 && all(purrr::map(lol, class) == 'brmsfit')) {
      if (T) {
        ## 2021-04-02 11:49 This shouldn't occur anymore
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
        stop(o_fn)
        file.remove(o_fn)
        return(NULL)
      }
      IE_tests <- purrr::map(lol, perform_bayesian_IE_test)
    } else {
      IE_tests <- purrr::map(lol, ~ .x[names(.x) != 'lm'])
    }

    ret_val <- tryCatch(rbindlist(IE_tests, fill = T) %>%
      .[, 'project_extended' := projects] %>%
      .[, 'analysis_idx' := analysis_idx],
      error = function(e) { print(e); browser(); NULL })
    return(ret_val)
  }, .parallel = (ncores > 1)) %>% rbindlist(fill = T)

  if (null_dat(dtf) || ncol(dtf) == 2) return(NULL)

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
  permute_id = NULL,
  focus_settings = NULL,
  check_data_availability = F,
  skip_missing = F,
  N_downsample_settings = NULL,
  exclude_rooney = TRUE,
  include_non_ds_tallies = (idxs_name != 'main_analysis')) {

  fa_flag <- paste0('-focus_hlas_', paste(hla_alleles,
      collapse = '_'))
  if (fa_flag == '-focus_hlas_A0201_A1101_B0702_B2705_B4001') {
    fa_flag <- ''
  }
  idxs_name <- attr(analysis_idxs, 'name') %||%
    gsub('_idxs$', '', deparse(substitute(analysis_idxs)))
  analysis_idxs_flag <- paste0('-analysis_idxs=', idxs_name)
  if (!is.null(ds_frac)) stopifnot(!is.null(iter))
  all_coef_fn <- file.path(rds_dir,
    glue::glue('all_coef_overviews\\
      {fa_flag}\\
      {make_flag(z_normalize)}\\
      {analysis_idxs_flag}\\
      {make_flag(reg_method)}\\
      {make_flag(ds_frac)}\\
      {make_flag(iter)}\\
      {make_flag(exclude_rooney)}\\
      {make_flag(include_non_ds_tallies)}\\
      {make_flag(N_downsample_settings)}\\
      {make_flag(permute_id)}\\
      .rds')
    )

  if (!IE_requires_computation(all_coef_fn) && !redo &&
      is.null(focus_settings) && !skip_missing) {
    return(
      format_overview_res(
        dtf = readRDS(all_coef_fn),
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

  if (exclude_rooney) {
    all_settings %<>% 
      dplyr::filter(analysis_name != 'rooney_param_titration')
  }

  if (!include_non_ds_tallies) {
    all_settings %<>% 
      dplyr::filter(analysis_name == 'twoD_sens_analysis')
  }

  if (!is.null(focus_settings)) {
    setDT(all_settings)
    sub_args <- focus_settings[intersect(names(focus_settings),
      colnames(all_settings))]
    setkeyv(all_settings, names(sub_args))
    all_settings <- all_settings[sub_args]
  }

  if (!is.null(N_downsample_settings)) {
    set.seed(1)
    idxs <- sample(1:nrow(all_settings), N_downsample_settings)
    all_settings <- all_settings[idxs, ]
  }

  ## Make sure donor summaries are available for all hla_alleles and
  ## analysis_idxs
  if (check_data_availability) {
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
  }

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

    ## Create single grid for all tumor types and this particular set
    ## of IE settings
    dtf <- prep_continuous_param_grid(
      focus_allele = r[['focus_allele']],
      LOH_HLA = r[['LOH_HLA']],
      analysis_name = r[['analysis_name']],
      overlap_var = r[['overlap_var']],
      patient_inclusion_crit = r[['patient_inclusion_crit']],
      analysis_idxs = analysis_idxs,
      # hla_sim_range = hla_sim_range,
      redo = redo_subanalyses,
      permute_id = permute_id,
      ds_frac = ds_frac,
      ncores = ncores,
      skip_missing = skip_missing,
      reg_method = reg_method,
      z_normalize = z_normalize
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

  if (is.null(focus_settings) && !skip_missing &&
    is.null(N_downsample_settings)) {
    saveRDS(all_coef_overviews, all_coef_fn)
    mymessage(
      msg = sprintf('Wrote all_coef_overviews to %s', all_coef_fn),
      instance = 'compile_all_coef_overview'
    )
  }

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
  print_overview_stats(dtf, stage_id = 'before any filtering')

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

  print_overview_stats(dtf, stage_id = 'after analysis type-filtering')

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
      print_overview_stats(dtf, stage_id = 'after picking representative allele')
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
      stage_id = 'After filtering for rc error'
    )
  }

  if (!is.null(delta_filter_error)) {
    dtf <- dtf[-log2_est_error >= delta_filter_error]
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After filtering for delta error'
    )
  }

  if (!single_hla_mode && pick_allele_method == 'est_error') {
    s_dtf <- copy(dtf[order(log2_est_error), .SD[1],
      by = analysis_grp_vars])
    dtf <- s_dtf
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After picking representative alleles'
    )
  }

  if ('intercept' %in% colnames(dtf) &&
      filter_intercept && dtf[intercept < 0, .N] > 0) {
    dtf <- dtf[intercept >= 0]
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'After removing negative intercept sub-analyses'
    )
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
