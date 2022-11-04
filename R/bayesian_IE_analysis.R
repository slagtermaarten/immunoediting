fit_bayesian_regression <- function(
  dtf, orig_fit = NULL, N_iter = 2000,
  intercept_resolution = .01,
  intercept_scale = NULL,
  compare_unbiased = F, model_name = 'biased',
  rlm_diff_threshold = NULL) {

  source(file.path(ma_dir, 'immune_editing', 'load_brms_models.R'))
  ## Save example data for use in generating scaffold models
  if (F) saveRDS(dtf, file.path(rds_dir, 'scaffold_rstan_dat.rds'))

  ## If orig_fit not supplied to function, try to look it up from
  ## scaffold
  if (is.null(orig_fit) || is.character(orig_fit)) {
    if (!grepl('_fit$', model_name)) {
      lu_model_name <- glue('{model_name}_fit') %>% as.character()
      orig_fit <- brm_scaffold[[lu_model_name]]
    } else if (model_name == 'intercept_matched') {
      orig_fit <- load_intercept_matched_brm(
        dtf = dtf,
        intercept_resolution = intercept_resolution
      )
    } else {
      lu_model_name <- model_name
      model_name <- gsub('_fit$', '', model_name)
      orig_fit <- brm_scaffold[[lu_model_name]]
    }
  }

  dtf <- apply_intercept_scale(
    dtf = dtf,
    intercept_scale = intercept_scale
  )
  if (is.null(dtf)) return(NULL)

  fit <- tryCatch({
    update(orig_fit, newdata = dtf, refresh = 0, iter = N_iter)
  }, error = function(e) { print(e); NULL })
  if (is.null(fit)) return(NULL)

  if (!is.null(rlm_diff_threshold)) {
    coef_frac_change <- compare_bayes_to_rlm(fit, prep)
    ## (1 + f_1) / (1 + f_0) = factor difference between first NAYR
    ## and second NAYR, wherein f_1 is fractional change in rc and f_0
    ## is fractional change in intercept.
    total_change_in_nayr <-
      (1 + coef_frac_change[2]) / (1 + coef_frac_change[1])

    if (abs(total_change_in_nayr - 1) > rlm_diff_threshold) {
      warning('Detected a large deviation from rlm')
      if (interactive()) {
        browser()
      }
    }
  }

  fit$model_name <- model_name
  fit$orig_b0 <- (attr(dtf, 'orig_b0') %||% NULL)
  return(fit)
}


#' Compute Bayes factors
#'
#'
perform_bayesian_IE_test <- function(fit) {
  ## Change in NAYR between PS == 0 and PS == 1 ->
  ## ((b0 + b1) - b0) / b0 = b1 / b0 < 0 is consistent with editing
  if (is.null(fit)) return(NULL)
  editing_test <- tryCatch(suppressWarnings({
    brms::hypothesis(x = fit, hypothesis = '(ol / Intercept) < 0')
  }), error = function(e) { print(e); browser() })

  out <- editing_test$hypothesis

  dtf <- fit$data %>% dplyr::filter(!is.na(y_var) & !is.na(ol))
  out$N_patients <-
    tryCatch(nrow(dtf), error = function(e) { NA })

  out$sampling_convergence <- check_rhat(fit$fit)
  out <- c(out, extract_brms_pars(fit))
  out$orig_b0 <- fit$orig_b0
  out$model_name <- fit$model_name
  names(out) <- make_names_df_friendly(names(out))
  return(out)
}



#' Simulate test data and fit a Bayesian model
#'
#' @param scale_intercept Passed to sim_data. Simulate data with
#' observed intercept multipled by \code{scale_intercept}
#' @param beta_x Passed to sim_data
#'
sim_and_test <- function(
  prep, rc = 0, noise_bound = .2,
  intercept_scale = NULL,
  scale_intercept = NULL,
  orig_fit = NULL, model_name = 'unbiased', ...) {

  dots <- list(...)

  s_dtf <- purrr::exec(sim_data,
    prep = prep, rc = rc,
    beta_x = beta_x,
    scale_intercept = scale_intercept,
    noise_bound = noise_bound)

  fit <- fit_bayesian_regression(
    dtf = s_dtf,
    orig_fit = orig_fit,
    intercept_scale = intercept_scale,
    # intercept_scale = NULL,
    model_name = model_name
  )
  if (is.null(fit)) return(NULL)

  test <- perform_bayesian_IE_test(fit)
  test$rc <- rc
  test$noise_bound <- noise_bound
  test$orig_b0 <- fit$orig_b0

  rc_test <- tryCatch(suppressWarnings({
    brms::hypothesis(x = fit, hypothesis = 'ol < 0')
  }), error = function(e) { print(e); browser() })

  # qplot(x = b_Intercept, y = b_ol, data = posterior_samples(fit))
  # qplot(x = b_ol / b_Intercept, data = posterior_samples(fit))
  # post_samples <- posterior_samples(fit, add_chain = TRUE)
  # sd(with(post_samples, b_ol / b_Intercept))
  # mean(with(post_samples, b_ol / b_Intercept))
  # test$intercept_estimate

  test$rc_error <- rc_test$hypothesis$Est.Error

  return(test)
}

plot_bayesian_regression_scatter <- function(dtf, fit, N_coefs = 200) {
  abline_coefs <- dplyr::sample_n(posterior_samples(fit), N_coefs)
  p <- qplot(data = dtf, x = ol, y = y_var) +
    geom_abline(
      data = abline_coefs,
      mapping = aes(intercept = b_Intercept, slope =  b_ol),
      colour = 'indianred3', alpha = .1)
}

sim_and_plot <- function(prep, rc = 0, noise_bound = .3) {
  s_dtf <- sim_data(
    prep = prep,
    rc = rc,
    noise_bound = noise_bound
  )
  fit <- fit_bayesian_regression(s_dtf, model_name = 'unbiased')
  # perform_bayesian_IE_test(fit)
  editing_test <- tryCatch(suppressWarnings({
    brms::hypothesis(x = fit, hypothesis = '(ol / Intercept) < 0')
  }), error = function(e) { print(e); browser() })
  p_hypo <- plot_hypothesis(editing_test)
  return(list(p, p_hypo))
}


#' Check R-hat
#'
#' url <- 'https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/divergences_and_bias/stan_utility.R'
#' devtools::source_url(url)
#'
check_rhat <- function(fit) {
  if (class(fit) == 'brmsfit') {
    fit <- fit$fit
  }
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  valid_chains = T
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      valid_chains <- F
    }
  }
  return(valid_chains)
}


extract_brms_pars <- function(fit, pars = c('Intercept', 'ol')) {
  purrr::map(pars, function(par) {
    suppressWarnings(summary(fit))$fixed[par, ] %>%
      { set_names(., paste0(tolower(par), '_', names(.))) } %>%
      { set_names(., make_names_df_friendly(names(.))) }
  }) %>% unlist
}


#' Permute data and do Bayesian parameter estimation for each permutation
#'
#'
bayes_negative_control <- result_cacher(
  f = function(args, N_repeats = 1000) {
    prep <- call_func(prep_cont_IE_analyses, args)

    neg_c <- furrr::future_map_dfr(1:(N_repeats+1), function(i) {
      l_dtf <- prep$dtf
      if (i != 1) {
        l_dtf <- l_dtf %>% mutate(y_var = sample(y_var))
      }
      b_lm <- fit_bayesian_regression(l_dtf)
      perform_bayesian_IE_test(b_lm)
    }, .id = 'i')

    ## Permutation p-value
    valid_perms <- setDT(neg_c) %>%
      .[sampling_convergence == T] %>%
      .[est_error <= 1] %>%
      { . }

    args %>%
      append(neg_c[1]) %>%
      append(list(
        'perm_p' = valid_perms[, mean(.SD[1, estimate] >= .SD[-1, estimate])],
        'N_perms' = as.integer(nrow(valid_perms))
      ))
  },
  filename = function() {
    file.path(rds_dir, 'bayes_negative_control',
      paste0(paste(args, collapse = '-'), make_flag(N_repeats), '.rds'))
  },
  min_mod_time = '2021-03-04 11:00'
)


