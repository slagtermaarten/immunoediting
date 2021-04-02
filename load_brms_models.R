pacman::p_load('brms')


basic_model <- function(dtf, prior = NULL, family = stats::gaussian()) {
  brms::brm(
    y_var ~ 0 + Intercept + ol,
    data = dtf,
    family = family,
    prior = prior,
    sample_prior = T,
    control = list(adapt_delta = 0.95)
  )
}


#' Abandoned function
#'
#'
gen_student_brms <- function(mean_y_var, mean_residual, sd_residual) {
  prior <- brms::get_prior(y_var ~ 0 + Intercept + ol,
    data = dtf, family = student(), control = list(adapt_delta = 0.95))
  # ## Alpha, beta (shape, rate) parameterization
  # ## Mean = alpha/beta, mode = (alpha - 1)/beta, variance = alpha/beta^2
  # ## Put average of intercept at mean_y_var
  # prior$prior[2] <- glue::glue('gamma({mean_y_var * conc}, {conc})')
  # ## Put average of nu (df) on 2
  # prior$prior[4] <- glue::glue('gamma({2 * conc}, {conc})')
  # fit <- basic_model(prep$dtf, prior, family = student())
}


gen_scaffold <- function() {
  ## Only create new env when it doesn't exist already in this session
  brm_scaffold <<- new.env(parent = baseenv())
  dtf <- readRDS(file.path(rds_dir, 'scaffold_rstan_dat.rds'))

  ## Explicitly model the intercept such that we define our own prior
  ## for it
  prior <- brms::get_prior(y_var ~ 0 + Intercept + ol,
    data = dtf, family = gaussian(),
    control = list(adapt_delta = 0.95))
  ## Intercept
  prior$prior[2] <- 'student_t(3, 0, 2.5)'
  ## Slope
  prior$prior[3] <- 'student_t(3, 0, 2.5)'
  ## sigma
  prior$prior[4] <- 'student_t(3, 0, 2.5)'
  assign('unbiased_fit',
    basic_model(dtf = dtf, prior = prior),
    envir = brm_scaffold)

  ## Unbiased, uninformative prior on intercept
  ## Intercept
  prior$prior[2] <- 'student_t(20, 0, 2.5)'
  assign('unbiased_ui_fit',
    basic_model(dtf = dtf, prior = prior),
    envir = brm_scaffold)

  biased_prior <- prior
  ## Intercept
  biased_prior$prior[2] <- 'exponential(.5)'
  assign('biased_fit',
    basic_model(dtf = dtf, prior = biased_prior),
    envir = brm_scaffold)

  t_prior <- brms::get_prior(y_var ~ 0 + Intercept + ol,
    data = dtf, family = student(),
    control = list(adapt_delta = 0.95))
  ## Intercept
  # t_prior$prior[2] <- 'student_t(3, 0, 2.5)'
  t_prior$prior[2] <- 'exponential(.5)'
  ## Nu
  t_prior$prior[4] <- 'gamma(20, .1)'
  assign('biased_robust_fit',
    basic_model(dtf = dtf, prior = t_prior, family = student()),
    envir = brm_scaffold)
}


peek_into_scaffold <- function() {
  for (o in ls(brm_scaffold)) {
    cat('\n\nObject :', o, '\n\n')
    print(brm_scaffold[[o]])
  }
}
# peek_into_scaffold()
# plot(brm_scaffold$neg_biased_fit)
# plot(brm_scaffold$neg_biased_fit)
# plot(brm_scaffold$unbiased_fit)
# stancode(brm_scaffold$biased_fit)


clear_scaffold <- function() {
  tryCatch(rm(list = ls(brm_scaffold), envir = brm_scaffold),
    error = function(e) { })
}


rds_scaffold <- function() {
  clear_scaffold()
  gen_scaffold()
  saveRDS(brm_scaffold, file.path(rds_dir, 'brm_scaffold.rds'))
  peek_into_scaffold()
}
# rds_scaffold()


if (!exists('brm_scaffold')) {
  brm_scaffold <<- readRDS(file.path(rds_dir, 'brm_scaffold.rds'))
  # stancode(brm_scaffold$unbiased_fit)
  # stancode(brm_scaffold$biased_fit)
  # stancode(brm_scaffold$biased_robust_fit)
}


find_fitting_intercept <- function(dtf, intercept_resolution = .01) {
  b0 <- coef(lm(y_var ~ ol, prep$dtf))[1]
  f <- 10^(-log10(intercept_resolution))
  unname(round(f * b0) / f)
}


#' Run Bayesian regression with intercept prior centered on it's expected
#' location. 
#'
#' @details The returned model is not necessarily fitted to the data in dtf but
#' does have the intercept centered at theќ
#'
#'
load_intercept_matched_brm <- result_cacher(
  f = function(dtf, model_name = 'unbiased', 
    intercept_resolution = .01) {
    intercept <- find_fitting_intercept(
      dtf = dtf,
      intercept_resolution = intercept_resolution
    )
    if (model_name == 'unbiased') {
      prior <- brms::get_prior(y_var ~ 0 + Intercept + ol,
        data = dtf, family = gaussian(),
        control = list(adapt_delta = 0.95))
      ## Intercept
      prior$prior[2] <- glue::glue('student_t(3, {intercept}, 2.5)')
      ## Slope
      prior$prior[3] <- 'student_t(3, 0, 2.5)'
      ## sigma
      prior$prior[4] <- 'student_t(3, 0, 2.5)'
    }
    mod <- basic_model(dtf = dtf, prior = prior)
    return(mod)
  },
  filename = function() {
    intercept <- find_fitting_intercept(
      dtf = dtf,
      intercept_resolution = intercept_resolution
    )
    file.path(rds_dir, 'brms_lm_intercept_titration',
      glue('unbiased\\
        {make_flag(model_name)}\\
        {make_flag(intercept)}\\
        {make_flag(intercept_resolution)}.rds'))
  },
  min_mod_time = '2021-03-14 20:13'
)
