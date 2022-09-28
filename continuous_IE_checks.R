
check_perm_q <- function(dtf) {
  ## Check whether the assumption of all tumor types having delta_n in
  ## all three boxes is true
  dtf[, .N, 
    keyby = list(
      project_extended, 
      cut(perm_delta_pq, c(0, 0.05, .95, 1), 
        include.lowest = TRUE))] %>%
  .[, .N == 3, by = project_extended] %>%
  .[, all(V1)]

  dtf[, .N, 
    keyby = list(project_extended, 
      cut(perm_delta_pq, c(0, 0.05, .95, 1), 
        include.lowest = TRUE))] %>%
  .[, N[1] / (N[1] + N[3]), by = project_extended] %>%
  .[order(V1)] %T>%
  .[, cat('Fraction favoring editing: ', mean(V1 > .5), '\n\n')] %>%
  { . }
}


filter_coef_overview_ <- function(dtf, code,
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

  cat('Overall decrease: ', 1-mean(bools), '\n\n')

  stats <- 
    purrr::map_dfr(names(comp_before), function(n) {
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
    dplyr::mutate(fc = 2^log2_fc) %>%
    dplyr::select(-log2_fc) %>%
    dplyr::rename(N_before = N.x, N_after = N.y)

  # return(list('dtf' = dtf, 'stats' = stats))
  attr(dtf, 'latest_stats') <- stats
  return(dtf)
}


filter_coef_overview <- function(
  dtf,
  print_messages = F,
  force_positive_intercept = F,
  force_intercept_ge_rc = F,
  min_patients = NULL,
  min_nz_patients = NULL,
  min_project_size = NULL,
  max_delta_n_cov = NULL,
  intercept_filter_p_val = NULL,
  intercept_filter_magnitude = NULL,
  force_sensible_perm_delta_CI = F,
  norm_scale_filter = NULL,
  scale_avg_norm_filter = NULL,
  delta_SE_filter = NULL,
  delta_c_SE_filter = NULL,
  AFDP_50_filter = NULL,
  AFDP_75_filter = NULL,
  AFDP_90_filter = NULL,
  NMADR_q50_filter = NULL,
  NMADR_q75_filter = NULL,
  adaptive_scale_filter = FALSE,
  adaptive_NMADR_q50_filter = FALSE,
  adaptive_NMADR_q75_filter = FALSE,
  adaptive_delta_SE_filter = FALSE,
  adaptive_delta_mean_c_SE_filter = FALSE,
  adaptive_norm_scale_filter = FALSE) {

  stopifnot(!any(duplicated(colnames(dtf))))

  ## VE-threshold of 5 leads to ~0 intercept
  dtf <- dtf[VE_threshold != '5']
  dtf[, VE_threshold := droplevels(VE_threshold)]

  # dtf <- dtf[, !duplicated(colnames(dtf)), with = F]
  if ('delta_n' %in% colnames(dtf)) {
    dtf <- dtf[order(delta_n), ]
  }

  if (null_dat(dtf)) {
    message('No rows before filtering')
    return(NULL)
  }

  dtf <- filter_coef_overview_(
    dtf = dtf,
    print_messages = print_messages,
    code = !is.na(rc) & !is.na(intercept) & converged
  )
  if (print_messages)
    print_overview_stats(dtf,
      plot_fishtails = plot_fishtails,
      stage_id = 'Filtering non-converged sub-analyses'
    )

  if (!is.null(max_delta_n_cov)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = delta_n_cov <= max_delta_n_cov
    )
    if (print_messages)
      print_overview_stats(
        dtf = dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Delta cov filtering'
      )
  }

  if (!is.null(intercept_filter_p_val)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = p_val_intercept <= intercept_filter_p_val
    )
    if (print_messages)
      print_overview_stats(
        dtf = dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Intercept p-val filtering'
      )
  }

  if (!is.null(intercept_filter_magnitude)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(intercept) >= intercept_filter_magnitude
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Intercept magnitude filtering'
      )
  }

  if (force_positive_intercept) {
    dtf <- filter_coef_overview_(
      print_messages = print_messages,
      dtf = dtf, code = intercept > 0
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Force positive intercept'
      )
  }

  if (force_intercept_ge_rc) {
    dtf <- filter_coef_overview_(
      print_messages = print_messages,
      dtf = dtf, code = abs(intercept) >= abs(rc)
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Force intercept >= rc'
      )
  }

  any_bool_filter <- ls(pattern = 'scale.*filter') %>%
    setdiff('norm_scale_filter') %>%
    map_lgl(~get(.x) %||% F) %>%
    any()
  if (any_bool_filter) {
    tmp <- dtf[
      maartenutils::eps(scale, 0, 1e-16) &
      intercept > 1e-2]
    if (nrow(tmp) > 0) {
      tmp[, .(delta_mean, n_patients, intercept, scale)]
      args <- tmp[1] %>% pan_IE_res_to_call_args
      print_plot(plot_rlm_pick(tmp[1]))
      stop('Implement zero error filtering!')
    }
  }

  if (force_sensible_perm_delta_CI) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = perm_delta_CI_l_median < 0 & perm_delta_CI_h_median > 0
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'Restrict to sensible permutation CI'
      )
  }

  if (!is.null(delta_SE_filter)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = delta_SE <= delta_SE_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After delta_SE filtering'
      )
  }

  if ('perm_delta_mean_c_SE' %in% colnames(dtf) &&
      !is.null(delta_c_SE_filter)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = perm_delta_mean_c_SE <= delta_c_SE_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After delta_c_SE filtering'
      )
  }

  if (!is.null(norm_scale_filter)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      # code = abs(norm_scale) <= norm_scale_filter
      code = abs(scale) <= norm_scale_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After norm_scale filtering'
      )
  }

  if (!is.null(scale_avg_norm_filter)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(scale_avg_norm) <= scale_avg_norm_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After scale_avg_norm filtering'
      )
  }

  if (!is.null(NMADR_q75_filter)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = NMADR_q75 <= NMADR_q75_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After NMADR_q75 filtering'
      )
  }

  if (!is.null(AFDP_50_filter)) {
    dtf <- filter_coef_overview_(
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
    dtf <- filter_coef_overview_(
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
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = AFDP_90 <= AFDP_90_filter
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After AFDP_90 filtering'
      )
  }

  if (!is.null(adaptive_norm_scale_filter) &&
      adaptive_norm_scale_filter &&
      dtf[, any(intercept < 0)]) {
    stop('Implement me')
    ## Determine what scale is not tolerable
    dtf[is.finite(scale_avg_norm) & intercept < 0]
    thresh <- dtf[is.finite(scale_avg_norm) & intercept < 0,
      min(abs(scale_avg_norm))]
    print(thresh)
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(norm_scale) < thresh
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After adaptive norm_scale filtering'
      )
  }

  if (!is.null(adaptive_scale_filter) &&
      adaptive_scale_filter &&
      dtf[, any(intercept < 0)]) {
    ## Determine what scale is not tolerable
    thresh <- dtf[intercept < 0, min(abs(norm_scale))]
    print(thresh)
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(norm_scale) < thresh
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After adaptive norm_scale filtering'
      )
  }

  if (!is.null(adaptive_NMADR_q50_filter) &&
      adaptive_NMADR_q50_filter &&
      dtf[, any(intercept < 0)]) {
    ## Determine what NMADR_q50 is not tolerable
    thresh <- dtf[intercept < 0, min(NMADR_q50)]
    cat('Required threshold inferred to be: ', thresh, '\n')
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(NMADR_q50) < thresh
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After adaptive NMADR_q50 filtering'
      )
  }

  if (!is.null(adaptive_NMADR_q75_filter) &&
      adaptive_NMADR_q75_filter &&
      dtf[, any(intercept < 0)]) {
    ## Determine what NMADR_q75 is not tolerable
    thresh <- dtf[intercept < 0, min(NMADR_q75)]
    cat('Required threshold inferred to be: ', thresh, '\n')
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = abs(NMADR_q75) < thresh
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After adaptive NMADR_q75 filtering'
      )
  }

  if (!is.null(adaptive_delta_SE_filter) &&
      adaptive_delta_SE_filter &&
      dtf[, any(intercept < 0)]) {
    ## Determine what scale is not tolerable
    thresh <- dtf[intercept < 0, min(delta_SE)]
    print(thresh)
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = delta_SE < thresh
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After adaptive delta_SE filtering'
      )
  }

  if (!is.null(adaptive_delta_mean_c_SE_filter) &&
      'perm_delta_mean_c_SE' %in% colnames(dtf) &&
      adaptive_delta_mean_c_SE_filter) {
    ## Determine what scale is not tolerable
    thresh <- dtf[abs(perm_delta_mean_c_mean) > 1, 
      min(abs(perm_delta_mean_c_SE))]
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = perm_delta_mean_c_SE < thresh
    )
    if (print_messages)
      print_overview_stats(dtf,
        plot_fishtails = plot_fishtails,
        stage_id = 'After adaptive perm_delta_mean_c_SE filtering'
      )
  }

  ## Apply these filters after all others to ensure we're counting
  ## 'valid' analyses rather than 'all' analyses
  if (!is.null(min_patients) &&
      dtf[, any(n_patients < min_patients)]) {
    dtf <- filter_coef_overview_(
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

  ## Apply these filters after all others to ensure we're counting
  ## 'valid' analyses rather than 'all' analyses
  # summary(dtf$n_nz_patients)
  if (!is.null(min_nz_patients) &&
      any(min_nz_patients > dtf$n_nz_patients)) {
    dtf <- filter_coef_overview_(
      dtf = dtf,
      print_messages = print_messages,
      code = n_nz_patients >= min_nz_patients
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
  ## Order projects by mean perm_delta_mean_pc
  if (F) {
    pe <- dtf[,
      wa_var(.SD, varname = 'perm_delta_mean_pc',
        error_var = 'norm_scale'),
      by = .(tumor_type, project_extended)] %>%
      .[order(V1), .(project_extended, tumor_type, V1)]
  } else if (F && 'perm_delta_mean_pc' %in% colnames(dtf)) {
    pe <- dtf[, mean(perm_delta_mean_pc),
      by = .(tumor_type, project_extended)] %>%
      .[order(V1), .(project_extended, tumor_type, V1)]
  } else if (F && 'delta' %in% colnames(dtf)) {
    pe <- dtf[, mean(delta),
      by = .(tumor_type, project_extended)] %>%
      .[order(V1), .(project_extended, tumor_type, V1)]
  } else if ('delta_n' %in% colnames(dtf)) {
    pe <- dtf[, median(delta_n),
      by = .(tumor_type, project_extended)] %>%
      .[order(V1), .(project_extended, tumor_type, V1)]
  }
  if (print_messages) print(pe)
  dtf[, project_extended := factor(project_extended,
    levels = pe$project_extended)]
  dtf[, tumor_type := factor(tumor_type, levels = pe$tumor_type)]

  ## Some sanity checks
  # stopifnot(dtf[, all(intercept > 0)])
  # browser(expr = dtf[, all(intercept > 0)])
  stopifnot(!dtf[, all(maartenutils::eps(intercept, 1, 1e-1))])
  dtf[, log10_n_patients := log10(n_patients)]

  return(dtf)
}



