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

  # out$baseline_yr <- coef(fit)[['(Intercept)']] + 
  #       mean(fit$data$i) * coef(fit)[['i']]

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


