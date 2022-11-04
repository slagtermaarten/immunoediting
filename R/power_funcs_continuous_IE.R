#' Also compute power the 'classical way', i.e. compute the probability of
#' obtaining a significant result given all other relevant parameters
#'
compute_power_continuous_IE_formal <- function() {
  prep <- prep_cont_IE_analyses(focus_allele = as.character(focus_allele),
                                redo = redo,
                                analysis_name = as.character(analysis_name),
                                tumor_type = tumor_type,
                                reg_method = reg_method,
                                # hla_sim_range = hla_sim_range,
                                repertoire_overlap_dat = repertoire_overlap_dat,
                                LOH_HLA = LOH_HLA,
                                idx = idx,
                                overlap_var = overlap_var,
                                patient_inclusion_crit = patient_inclusion_crit)
  if (null_dat(prep$obj)) return(NULL)

  repertoire_overlap_dat <- get_repertoire_overlap(ro = repertoire_overlap_dat,
                                                   hla_allele = focus_allele)

  IE_tests <- test_continuous_IE(focus_allele = as.character(focus_allele),
                                 tumor_type = tumor_type,
                                 redo = redo,
                                 repertoire_overlap_dat = repertoire_overlap_dat,
                                 analysis_name = analysis_name,
                                 reg_method = reg_method,
                                 LOH_HLA = as.character(LOH_HLA),
                                 # hla_sim_range = hla_sim_range,
                                 idx = idx,
                                 return_res = T,
                                 fn_suffix = fn_suffix,
                                 overlap_var = overlap_var,
                                 patient_inclusion_crit = patient_inclusion_crit)
  stats <- IE_tests$stats
  stats <- stats[which(!sapply(stats, is.null))]
  projects <- intersect(prep$projects, names(stats))

  res <- plyr::llply(prep$projects, function(pe) {
    l_stat <- stats[[pe]][[stats_idx]]
    if (is.null(l_stat) || 'error' %in% class(l_stat)) return(NULL)

    l_dtf <- subset_project(prep$obj, pe)
    l_dtf <- normalize_columns(l_dtf)
    library(pwr)
    if (F) {
      ## Fix effect size and retrieve power
      effect_sizes <- as.vector(outer(c(1, 2.5, 5, 7.5, 9), 10^-c(0:5))) %>%
        sort %>% { .[. <= 1] }
      conv_p_test <- pwr.f2.test(f2 = effect_sizes,
                                 u = u, v = l_stat$n_patients - u - 1,
                                 sig.level = p_val_bound)
    } else {
      ## Fix power and retrieve effect size f2
      power_vals <- seq(0.1, 1, by = .05)

      ## Cohen's f-squared to correlation coefficient for MCA
      ## Cohen, J. (1988). Statistical power analysis for the behavioral
      ## sciences (2nd ed.). Hillsdale,NJ: Lawrence Erlbaum.
      f2_to_r <- function(f2 = .35) {
        #' f2 <- r2 / (1 - r2)
        r2 <- f2 / (1 + f2)
        sqrt(r2)
      }
      y_var_sd <- l_dtf[, unique(y_var_sd)]
      intercept <- IE_tests$stats[[pe]][[stats_idx]]$intercept
      rc <- IE_tests$stats[[pe]][[stats_idx]]$rc

      conv_p_test <- purrr::map(auto_name(power_vals), function(p) {
        pwr.f2.test(power = p,
                    u = u,
                    v = l_stat$n_patients - u - 1,
                    sig.level = p_val_bound)$f2 %>% f2_to_r
        })

      conv_p_test <- purrr::map(conv_p_test, function(x) {
        if (!is.null(x) & !is.na(x)) {
          ## (intercept + 1 * x * y_var_d) / intercept (in which x <= 0)
          ## ==> 1 + x * y_var_d / intercept
          return(1 + x * y_var_sd / abs(intercept))
        } else {
          return(NULL)
        }})

      return(c(list('project_extended' = pe), conv_p_test,
               list('y_var_sd' = y_var_sd, 'intercept' = intercept)))
    }
    return(conv_p_test)
  }, .parallel = (ncores > 1)) %>% rbindlist(fill = T)
  return(res)
}
formals(compute_power_continuous_IE_formal) <- formals(test_continuous_IE)
formals(compute_power_continuous_IE_formal)$p_val_bound <- 0.05
formals(compute_power_continuous_IE_formal)$u <- 1




#' Compute quantity to be minimized for power computation
#'
#'
compute_diff_from_desired_p <- function(es,
                                        dtf, stat, focus_allele, 
                                        fill_var = 'rc',
                                        stats_idx = 1, p_val_bound = .05,
                                        PPV = .45,
                                        PS_accuracy = 1,
                                        reg_method = 'rlm_one_group') {
  # dtf$y_var <- dtf$y_var_bak
  # dtf$y_var_bak <- dtf$y_var
  ## TODO consider deducting both CYT and PS effects to arrive at cleaner
  ## 'baseline' y_var
  if (fill_var == 'rc') {
    es_PS <- es
    es_CYT <- 0
  } else {
    es_PS <- 0
    es_CYT <- es
  }
  # dtf$y_var <- dtf$y_var_bak
  dtf$y_var_bak <- dtf$y_var
  # mean(dtf$y_var, na.rm = T)
  dtf$y_pred <- predict(stat$lm, dtf)

  dtf$y_error <- dtf$y_pred - dtf$y_var
  # median(log2(dtf$y_pred) - log2(dtf$y_var), na.rm = T)
  # median(dtf$y_error, na.rm = T)

  # plot(dtf$y_var, dtf$y_pred)
  # dtf %>% dplyr::filter(error > 0)
  # dtf %>% dplyr::filter(error < 0)
  # if (T) {
  #   for (i in 1) {
      # print(i)
      # dtf$y_var <- dtf$y_var - (stat[['rc']] + es_PS) * dtf$o -
      #   (stat[['rc_CYT']] + es_CYT) * dtf$CYT
      # dtf$y_var <-  coef(stat$lm)[1] + dtf$y_error -
      #   (stat[['rc']] + es_PS) * dtf$o -
      #   (stat[['rc_CYT']] + es_CYT) * dtf$CYT
      ## 2019-12-17 07:46 New way of computing new y_var
      # dtf$y_var <-  coef(stat$lm)[1] + dtf$y_error -
      #   PPV * (es_PS * dtf$o - es_CYT * dtf$CYT)
      ## 2020-07-07 17:30 New way of computing new y_var
      dtf$y_var <- coef(stat$lm)[1] + dtf$y_error -
        PS_accuracy * PPV * (es_PS * dtf$o + es_CYT * dtf$CYT)
      stat_test <- tryCatch(compute_continuous_IE_statistics(
          obj = dtf, reg_method = reg_method, focus_allele = focus_allele,
          debug = F, partition_vars = c()),
        error = function(e) { print(e); NULL })
      if (is.null(stat_test)) return(1)
      # print(stat_test[[stats_idx]]$lm)
    # }
  # }
  p_var <- switch(fill_var, 'rc' = 'p_val', 'rc_CYT' = 'p_val_CYT')
  ## Squared loss function
  return((stat_test[[stats_idx]][[p_var]] - p_val_bound)^2)
}


optim_compute_power_continuous_IE <- function() {
  if (eps(p_val_bound, 0, 1e-12)) {
    return(NULL)
  }
  ## Determine name of outputfile
  o_fn <- gen_cont_IE_fn(base_name = 'cont-IE',
                         focus_allele = focus_allele,
                         analysis_name = analysis_name,
                         sub_dir = 'cont-IE-power_analysis',
                         # hla_sim_range = hla_sim_range,
                         tumor_type = tumor_type,
                         reg_method = reg_method,
                         LOH_HLA = LOH_HLA,
                         fill_var = fill_var,
                         idx = idx,
                         fn_suffix = glue::glue('power_analysis-\\
                                                ppv_{PPV}{fn_suffix}'),
                         p_val_bound = p_val_bound,
                         overlap_var = overlap_var,
                         patient_inclusion_crit = patient_inclusion_crit)

  if (!optim_IE_requires_computation(o_fn) && !redo) {
    if (return_res %||% T) {
      return(readRDS(o_fn))
    } else {
      return(NULL)
    }
  }

  repertoire_overlap_dat <- get_repertoire_overlap(ro = repertoire_overlap_dat,
                                                   hla_allele = focus_allele)

  prep <- prep_cont_IE_analyses(focus_allele = as.character(focus_allele),
                                redo = redo,
                                analysis_name = as.character(analysis_name),
                                tumor_type = tumor_type,
                                reg_method = reg_method,
                                # hla_sim_range = hla_sim_range,
                                repertoire_overlap_dat = repertoire_overlap_dat,
                                LOH_HLA = LOH_HLA,
                                idx = idx,
                                overlap_var = overlap_var,
                                patient_inclusion_crit = patient_inclusion_crit)

  if (null_dat(prep$obj)) return(NULL)

  IE_tests <- test_continuous_IE(focus_allele = as.character(focus_allele),
                                 tumor_type = tumor_type,
                                 redo = redo,
                                 repertoire_overlap_dat = repertoire_overlap_dat,
                                 analysis_name = analysis_name,
                                 reg_method = reg_method,
                                 LOH_HLA = as.character(LOH_HLA),
                                 # hla_sim_range = hla_sim_range,
                                 z_normalize = F,
                                 idx = idx,
                                 return_res = T,
                                 fn_suffix = fn_suffix,
                                 overlap_var = overlap_var,
                                 patient_inclusion_crit = patient_inclusion_crit)

  stats <- IE_tests$stats
  stats <- stats[which(!sapply(stats, is.null))]
  projects <- intersect(prep$projects, names(stats))

  res <- plyr::llply(prep$projects, function(pe) {
    # pe <- 'Pan*\'-\'*cancer'
    l_stat <- stats[[pe]][[stats_idx]]
    # vcov(l_stat$lm)
    if (is.null(l_stat)) return(NULL)

    optim_grid_size <- 1e-2
    optim_grid_size <- 1e-3
    x <- seq(0, 3, by = optim_grid_size)
    # x <- seq(0, 3, by = 1e-1)
    l_dtf <- subset_project(prep$obj, pe)
    # l_dtf <- normalize_columns(l_dtf)

    l_dtf <- l_dtf[is.finite(y_var) & is.finite(ol) & is.finite(weight)]
    if (nrow(l_dtf[y_var > 0]) <= 10) {
      return(NULL)
    }

    ## Test if there are any non-zero values
    if (all(l_dtf$y_var == 0 | is.na(l_dtf$y_var))) {
      return(NULL)
    }

    y_var_sd <- unique(prep$obj$y_var_sd)

    ## Determine initial solution by minimizing over a 1D grid of values
    y <- tryCatch(vapply(x, function(xv)
      compute_diff_from_desired_p(xv, dtf = l_dtf, stat = l_stat,
                   fill_var = fill_var,
                   focus_allele = focus_allele,
                   PPV = PPV,
                   reg_method = reg_method,
                   stats_idx = stats_idx,
                   p_val_bound = p_val_bound), numeric(1)),
                  error = function(e) { print(e); browser() })

    # compute_diff_from_desired_p(x[2496], dtf = l_dtf, stat = l_stat,
    #              fill_var = fill_var,
    #              focus_allele = focus_allele,
    #              PPV = PPV,
    #              reg_method = reg_method,
    #              stats_idx = stats_idx,
    #              p_val_bound = p_val_bound)
    # plot(l_dtf$ol, l_dtf$y_var)
    # plot(x, y, type = 'l')

    if (all(is.na(y))) {
      mywarning('Optimization problem, all ys NA',
          glue::glue('{o_fn}-{pe}'))
      perm_browser()
    }

    ## 2019-11-08 17:15 Inform when optimization landscape is flat and forego of
    ## optimizing further.
    if (eps(var(y, na.rm = T), 0, eps = 1e-6) &&
        length(setdiff(unique(y), NA)) == 1) {
      ## This is probably due to applying the fitted regression coefficient
      ## already causes a significant regression (negative but close to zero, in
      ## combination with large patient numbers)
      if (!is.na(y[1]) && maartenutils::eps(y[1], p_val_bound^2)) {
        ## Zero editing strength already suffices
        ret_val <- l_stat %>%
          { .[2:length(.)] } %>%
          c('es' = 0,
            'message' = 'flat_optimization_landscape',
            'par' = 0,
            'project_extended' = pe)
        return(ret_val)
      } else if (!is.na(y[1])) {
        ## No tested value suffices
        ret_val <- l_stat %>%
          { .[2:length(.)] } %>%
          c('es' = 1,
            'message' = 'flat_optimization_landscape',
            'par' = 1,
            'project_extended' = pe)
        return(ret_val)
      } else {
        erm_browser()
        mywarning('Optimization problem, unknown cause',
            glue::glue('{o_fn}-{pe}'))
        mail_notify(subject = 'Unknown problem with IE optimization',
          msg = glue::glue('{o_fn}-{pe}'))
      }
    }
    # max_y <- max(y, na.rm = T)
    min_y <- min(y, na.rm = T)
    ## 2019-11-03 11:23 Return the highest x that gives the lowest y
    ## 2019-11-08 16:35 Bug fix herein
    init <- x[y == min_y] %>% first

    if (fill_var == 'rc_CYT') {
      ub <- 10
    } else {
      # ub <- max(2, init * 2, na.rm = T)
      ub <- max(init * 8, .1)
    }

    ## Optimize starting from grid-identified initial solution
    i <- 0
    maxi <- 10
    res <- NULL
    while ((is.null(res) || eps(res$par, ub, 1e-5)) && i <= maxi) {
      i <- i + 1
      pacman::p_load('pso')
      elapsed_time <- system.time(res <- tryCatch({
        pso::psoptim(init, compute_diff_from_desired_p,
                     dtf = l_dtf,
                     stat = l_stat,
                     focus_allele = focus_allele,
                     stats_idx = stats_idx,
                     fill_var = fill_var,
                     PPV = PPV,
                     p_val_bound = p_val_bound,
                     reg_method = reg_method,
                     control = list(abstol = 1e-6),
                     # method = 'Brent',
                     # method = 'Nelder-Mead',
                     # control = list(maxit = 1000, abstol = 1e-6),
                     lower = c(0),
                     ## 2019-09-20 12:39 Adapted upper bound of search space,
                     ## minimum of 0.05
                     upper = i * ub)
      }, error = function(e) { print(e) }))
    }

    if (elapsed_time['elapsed'] > 10 && debug_optim) {
      et <- elapsed_time['elapsed']
      mail_notify(
        subject = glue::glue('Optimization is taking too long ({et})'),
        msg = glue::glue('{o_fn}-{pe}'))
    }

    if (eps(res$par, maxi * ub, 1e-5)) {
      mywarning(glue::glue('Optimization problem'), glue::glue('{o_fn}-{pe}'))
      mail_notify(subject = 'Optimization solution converged to upper boundary',
        msg = glue::glue('{o_fn}-{pe}'))
      # browser()
    }

    ## If solution lies further than grid size away from the initial solution,
    ## something is wrong. Optimization space is not necessarily convex
    ## (determined emperically) so this test does not make sense
    if (F && !eps(res$par, init, optim_grid_size)) {
      mywarning(glue::glue('Optimization problem'), glue::glue('{o_fn}-{pe}'))
      mail_notify(subject = 'Unexpected optimization solution',
        msg = glue::glue('{o_fn}-{pe}'))
      perm_browser()
    }
    # return(c(list('project_extended' = pe), conv_p_test,
    #          list('y_var_sd' = y_var_sd, 'intercept' = intercept)))
      # conv_p_test <- purrr::map(conv_p_test, function(x) {
      #   if (!is.null(x) & !is.na(x)) {
      #     ## (intercept + 1 * x * y_var_d) / intercept (in which x <= 0)
      #     ## ==> 1 + x * y_var_d / intercept
      #     return(1 + x * y_var_sd / abs(intercept))
      #   } else {
      #     return(NULL)
      #   }})

    ret_val <- l_stat %>%
      { .[2:length(.)] } %>% ## Exclude the 'lm' call
      ## 2019-09-20 12:39 Check convergence of optimizer and set to NA if
      ## convergence wasn't reached (success := 0)
      c('es' = ifelse(res$convergence == 0, res$par, NA),
        'par' = res$par, 'project_extended' = pe)
      # c('r' = setNames(conv_p_test[['power']], conv_p_test[['r']]))
    return(ret_val)
  }, .parallel = (ncores > 1)) %>% rbindlist(fill = T)

  saveRDS(res, o_fn)
  if (verbose) {
    mymessage(glue::glue('Writing {o_fn}'), 'optim_compute_power_continuous_IE')
  }
  return(res)
}
formals(optim_compute_power_continuous_IE) <- formals(test_continuous_IE)
formals(optim_compute_power_continuous_IE)$debug_optim <- F
formals(optim_compute_power_continuous_IE)$verbose <- T
formals(optim_compute_power_continuous_IE)$p_val_bound <- 0.05
formals(optim_compute_power_continuous_IE)$fill_var <- 'rc'
formals(optim_compute_power_continuous_IE)$PPV <- .45
# formals(optim_compute_power_continuous_IE)$adaptive_p_val_bound <- T


#' Block of code moved over from prep_pan_IE_heatmap
#'
#'
compute_power_block <- rlang::expr({
  if (include_power) {
    ## Adjust the p_val_bound
    if (adaptive_p_val_bound) {
      ## Overwrite p_val_bound with the maximum unadjusted p-value that is still
      ## 'significant' after MTC.
      ## This will be the new threshold for individual tests to be significant
      ## in the following power analysis
      signif_subset <- dplyr::filter(dtf, '{p_var.adj}' <= p_val_bound)
      if (null_dat(signif_subset)) {
        message('None of the analyses are significant after multiple ',
                'testing correction')
        p_val_bound <- .05
      } else {
        p_val_bound <- signif_subset %>% pull(p_val) %>% max(na.rm = T)
      }
    }

    ## Annotate the analyses with statistical power
    if (!is.null(p_val_bound) && is.finite(p_val_bound)) {
      if (F) {
        power_fun <- compute_power_continuous_IE_formal
        merge_cols <- analysis_id_cols
      } else {
        power_fun <- optim_compute_power_continuous_IE
        merge_cols <- c('focus_allele', analysis_id_cols)
      }
      power_analyses <-
        dtf[, c('focus_allele', analysis_id_cols), with = F] %>%
        cbind('tumor_type' = tumor_type) %>%
        unique %>%
        { cbind(., 'i' = 1:nrow(.)) } %>%
        { cbind(., 'total' = nrow(.)) } %>%
        # dplyr::filter(i == 1366) %>%
        plyr::adply(1, function(r) {
          i <- r[['i']]
          total <- r[['total']]
          mymessage(instance = 'Assembling power analyses',
                    msg = glue::glue('Processing {i}/{total}'))

          settings <- as.list(r)
          if (is.null(settings)) browser()

          if (use_PS_accuracy) {
            PS_accuracy <- PS_r2[settings[['focus_allele']], r2]
          } else {
            PS_accuracy <- 1
          }

          repertoire_overlap_dat <- get_repertoire_overlap(
            ro = repertoire_overlap_dat,
            hla_alleles = settings[['focus_allele']])

          ## The 'es' variable in this data.table will hold effective size
          ## required for the corresponding fill_var
          power_dat <- tryCatch(power_fun(
            focus_allele = settings[['focus_allele']],
            overlap_var = settings[['overlap_var']],
            tumor_type = tumor_type,
            patient_inclusion_crit = settings[['patient_inclusion_crit']],
            p_val_bound = p_val_bound,
            LOH_HLA = settings[['LOH_HLA']],
            PS_accuracy = PS_accuracy,
            PPV = PPV,
            repertoire_overlap_dat = repertoire_overlap_dat,
            # hla_sim_range = settings[['hla_sim_range']],
            analysis_name = settings[['analysis_name']],
            idx = settings[['analysis_idx']],
            redo = F,
            ncores = 1,
            fn_suffix = fn_suffix), error = function(e) { print(e); NULL })
          return(power_dat)
        }, .parallel = (ncores > 1))
      dtf <- controlled_merge(dtf, power_analyses, by_cols = merge_cols)
    }
  }
})
