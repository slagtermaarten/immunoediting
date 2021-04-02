theme_set(
  theme_fas(
    legend.box.background = element_blank(),
    legend.key.width = unit(5, 'mm'),
    panel.spacing = unit(.4, 'lines'),
    # base_size = font_size,
    legend.direction = 'vertical',
    legend.position = c(.95, .95),
    legend.justification = c(1, 1)
  )
)

tts <- names(tumor_types) %>% { .[!grepl('Cutaneous', .)] }
if (F) {
  tts <- grep('Stomach~EBV', tts, value = T)[2]
}
suppressWarnings(rm(all_peps))


if (F && !exists('ov')) {
  source(file.path(ma_dir, 'immune_editing',
      'continuous_IE_detection_helpers.R'))
  ov <- compile_all_coef_overview(
    redo = F,
    ncores = ifelse(T, 1, 40),
    z_normalize = F,
    idxs = main_analysis_idxs,
    repertoire_overlap_dat = repertoire_overlap_dat) %>%
    { cbind(., ds_param_grid[.$analysis_idx, ]) }
  ov_ao <- compile_all_coef_overview(
    redo = F,
    ncores = ifelse(T, 1, 40),
    z_normalize = F,
    idxs = affinity_idxs,
    repertoire_overlap_dat = repertoire_overlap_dat) %>%
    { cbind(., ds_param_grid[.$analysis_idx, ]) }
  all_preps <- get_all_preps(idxs = main_analysis_idxs) %>%
    { cbind(., ds_param_grid[.$analysis_idx, ]) }
  all_preps_ao <- get_all_preps(idxs = affinity_only_idxs) %>%
    { cbind(., ds_param_grid[.$analysis_idx, ]) }
}


sub_repetitions <- function(v, sub_val = '') {
  v[duplicated(v)] <- sub_val
  v
}
stopifnot(sub_repetitions(c('a', 'a', 'a', 'b', 'b')) ==
          c('a', '', '', 'b', ''))


display_settings <- list(
  analysis_name = list(
    breaks = analysis_names[c(1, 3, 4, 5)],
    labels = c('All variants', 'Clonal',
               'No drivers/essentials', 'Marty oncogenic')),
  LOH_HLA = list(
    breaks = c('no_LOHHLA', 'LOHHLA', 'strict_LOHHLA'),
    labels = c('None', 'LOHHLA', 'High confidence LOHHLA')),
  patient_inclusion_crit = list(
    breaks = c('', 'TR'),
    labels = c('Included', 'Excluded')),
  overlap_var = list(
    breaks = c('mean_score', 'mean_score_AB'),
    labels = c('A, B & C alleles', 'A & B alleles')),
  sts_filtering = list(
    breaks = c('TRUE', 'FALSE'),
    labels = c('STS', 'no STS')
  )
)


var_names <-
  c('Variant\nselection' = 'analysis_name',
    'Presentation\nscore' = 'overlap_var',
    'IT-resistant\npatients' = 'patient_inclusion_crit',
    'Adapting PS\nto LOH in HLA' = 'LOH_HLA',
    'Similarity\nto self-filter' = 'sts_filtering',
    'Affinity rank\npercentile threshold' = 'percentile_rank',
    'Expression\nfiltering' = 'exp_filter') %>%
  { set_names(names(.), .) }

# var_names <-
#   c('Variant selection' = 'analysis_name',
#     'Presentation score' = 'overlap_var',
#     'IT-resistant patients' = 'patient_inclusion_crit',
#     'LOH in HLA' = 'LOH_HLA',
#     'Similarity to self filter' = 'sts_filtering',
#     'Peptide binding affinity\nrank percentile threshold' = 'percentile_rank',
#     'Expression filtering' = 'exp_filter') %>%
#   { set_names(names(.), .) }


#' Take a row from the output created by prep_pan_IE_heatmap and
#' convert it to be arguments to be backfed to to test_continuous_IE
#'
pan_IE_res_to_call_args <- function(res_oi) {
  ## Complement with all other required arguments to
  ## test_continuous_IE
  call_args <- res_oi %>%
    as.list %>%
    { .[intersect(names(.), names(formals(test_continuous_IE)))] } %>%
    c(list('include_call' = T, 'redo' = F)) %>%
    {
      modifyList(.,
        list(
          hla_allele = as.character(.$focus_allele %||% .$hla_allele),
          project_extended = as.character(.$project_extended),
          tumor_type = as.character(.$project_extended),
          analysis_name <- as.character(.$analysis_name)
        )
      )
    }
  call_args$analysis_name <- as.character(call_args$analysis_name)
  return(call_args)
}


#' The fundamental and elementary visualisation of our 'continuous'
#' IE-detection method
#'
plot_PS_vs_yr <- function(
  dtf, lm,
  reg_method = 'rlm_one_group',
  # colour_vars = 'anchor_present_f',
  colour_vars = 'hla_allele_status',
  # shape_var = 'anchor_present_f',
  shape_var = 'hla_allele_status',
  include_legend = T,
  plot_density = dtf[, .N >= 1000],
  return_grob = F,
  point_alpha = .5,
  N_coefs = 200,
  hla_allele = attr(dtf, 'ds_opts')$hla_allele,
  font_size = 6,
  y_label = attr(dtf, 'y_label') %||% '',
  x_label = attr(dtf, 'x_label') %||%
    sprintf('%s presentation score', quickMHC::ppHLA(hla_allele))) {

  bayesian_mode <- grepl('bayesian', reg_method)

  if (null_dat(dtf) || is.null(lm)) return(NULL)
  if (plot_density) {
    colour_vars <- NULL
    shape_var <- NULL
  }
  x_min <- dtf[, floor(10 * min(ol, na.rm = T)) / 10]
  x_max <- dtf[, ceiling(10 * max(ol, na.rm = T)) / 10]
  x_range <- seq(x_min, x_max, by = .01)
  x_breaks <- seq(x_min, x_max, by = .3)
  y_min <- min(setdiff(dtf[, round(y_var, 1)], c(NA, NaN)), 0)
  # y_max <- min(max(setdiff(dtf[, round(y_var, 1)], c(NA, NaN))), 1)
  y_max <- coef(lm) %>% { c(.[1], .[1] + .[2]) } %>%
    max %>% { . * 1.05 }
  y_max <- max(dtf$y_var, na.rm = T)

  if (reg_method == 'lm') {
    scatter_dat <- dtf[hla_allele_status_b == F |
      is.na(hla_allele_status_b)]
  } else {
    scatter_dat <- dtf
  }

  if (!bayesian_mode) {
    ## Take the median CYT as the CYT over the entire range
    CYT <- if ('CYT' %in% colnames(scatter_dat))
      scatter_dat[, median(CYT, na.rm = T)]
    else
      NA
    predicted_dat <- data.frame(ol = x_range, CYT = CYT) %>%
      { cbind(., 'y_var' = predict(lm, ., se.fit = T)) } %>%
      dplyr::rename(y_var = y_var.fit) %>%
      dplyr::mutate(y_var_CI_l = y_var - 1.96 * y_var.se.fit) %>%
      dplyr::mutate(y_var_CI_h = y_var + 1.96 * y_var.se.fit) %>%
      dplyr::mutate(y_var_CI_l = pmax(0, y_var_CI_l)) %>%
      dplyr::mutate(weight = 1)
  } else {
    if (is.null(N_coefs)) {
      predicted_dat <- data.frame(ol = rep(x_range, each = 10)) %>%
        { cbind(., predict(lm, newdata = .)) }
      colnames(predicted_dat) <-
        make_names_df_friendly(colnames(predicted_dat))
      predicted_dat <- predicted_dat %>%
        dplyr::rename(y_var = estimate) %>%
        dplyr::rename(y_var_CI_l = q2_5) %>%
        dplyr::rename(y_var_CI_h = q97_5) %>%
        dplyr::mutate(weight = 1) %>%
        group_by(ol) %>%
        summarize(across(everything(), mean)) %>%
        dplyr::mutate(y_var_CI_l = pmax(0, y_var_CI_l)) %>%
        { . }
    }
  }

  # scatter_dat <- scatter_dat %>%
  #   dplyr::mutate(anchor_present_f = factor(ifelse(anchor_present == T,
  #       sprintf('%s+', quickMHC::ppHLA(focus_allele)),
  #       sprintf('%s-', quickMHC::ppHLA(focus_allele))))) %>%
  #   dplyr::filter(!is.na(weight) & weight > 0)

  if (is.null(colour_vars)) {
    col_title <- ''
  } else {
    if (colour_vars == 'IE_essentiality_impaired') {
      col_title <- 'Impaired T-cell\nsensitivity'
    } else if (colour_vars == 'CYT') {
      col_title <- 'Cytolytic\nscore'
    } else if (colour_vars %in% c('anchor_present_f', 'hla_allele_status')) {
      col_title <- 'Patient carries\nfocus allele'
    }
  }

  p <- ggplot(scatter_dat, aes_string(x = 'ol', y = 'y_var',
      colour = colour_vars, shape = shape_var))

  if (!plot_density) {
    p <- p + geom_point(alpha = point_alpha)
    if (!is.null(colour_vars)) {
      l_cols <- pipeline_evaluation_cols[1:(length(unique(dtf[[colour_vars]])))] %>%
        setNames(unique(dtf[[colour_vars]]))
      l_shapes <- setNames(c(15:17), as.character(unique(dtf[[colour_vars]])))
      p <- p +
        scale_colour_manual(name = col_title, values = l_cols) +
        scale_shape_manual(name = 'Patient carries\nfocus allele',
                           values = l_shapes)
    }
  } else {
    p <- p + geom_hex(bins = min(100, dtf[, .N] / 3)) +
      scale_fill_viridis_c(name = '# patients')
  }

  if (bayesian_mode && !is.null(N_coefs)) {
    abline_coefs <- dplyr::sample_n(posterior_samples(lm), N_coefs)
    p <- p + geom_abline(
      data = abline_coefs,
      mapping = aes(intercept = b_Intercept, slope =  b_ol),
      colour = 'indianred3', alpha = .1)
  } else {
    p <- p +
      geom_ribbon(data = predicted_dat,
        mapping = aes_string(ymin = 'y_var_CI_l', ymax = 'y_var_CI_h',
                             size = NULL, colour = NULL, shape = NULL),
        fill = 'firebrick3', show.legend = F, alpha = .2) +
      geom_line(data = predicted_dat, aes(shape = NULL),
        colour = 'firebrick3', show.legend = F)
  }

  p <- p +
    scale_x_continuous(x_label,
      breaks = x_breaks, expand = c(0.01, 0)) +
    # scale_y_continuous(name = attr(dtf, 'y_label'),
    #                    limits = c(y_min, y_max),
    #                    expand = c(0.01, 0)) +
    ylab(y_label) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    guides(shape = guide_legend(title = col_title),
           colour = guide_legend(title = col_title)) +
    maartenutils::gg_legend_alpha_cancel

  if (!include_legend) {
    p <- p + theme(legend.position = 'none')
  }

  if (reg_method != 'lm') {
    if (return_grob) {
      return(ggplotGrob(p))
    } else {
      return(p)
    }
  } else {
    marginal_boxplot <- ggplot(dtf[hla_allele_status_b == T & y_var <= 1],
      aes_string(x = '1', y = 'y_var')) +
      geom_boxplot(alpha = .5, width = .5, outlier.size = .5,
        fill = 'firebrick3') +
      scale_x_continuous(name = 'Anchor carriers', breaks = c(1)) +
      scale_y_continuous(name = '', limits = c(0, y_max)) +
      theme_fas(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

    ## Combine panels
    right_mar_s <- theme_fas()$plot.margin
    right_mar_s[[2]] <- unit(0, 'cm')
    g1 <- ggplotGrob(p + ggplot2::theme(plot.margin = right_mar_s))
    left_mar_s <- theme_fas()[['plot.margin']]
    left_mar_s[[4]] <- unit(0, 'cm')
    g2 <- ggplotGrob(marginal_boxplot +
      ggplot2::theme(plot.margin = left_mar_s))
    comb <- gridExtra::arrangeGrob(g1, g2, nrow = 1, widths = c(5, 1))
    return(comb)
  }
}


#' Identify some striking analyses from a prep_pan_IE_heatmap and plot scatter
#' plots for those
#'
#'
wrapper_plot_most_interesting_PS_vs_yr <-
  function(tumor_type = "Pan*'-'*cancer", 
    dtf = NULL, p_mode = 1, ...) {

  if (null_dat(dtf)) {
    dtf <- prep_pan_IE_heatmap(tumor_type = tumor_type,
                               adaptive_p_val_bound = F,
                               # redo = T,
                               redo = F,
                               # ncores = 1,
                               ncores = 40,
                               PPV = .45,
                               z_normalize = F,
                               # repertoire_overlap_dat = readRDS(ro_fn),
                               # debug_func = T,
                               # debug_func = F,
                               # p_val_bound = .05,
                               p_val_bound = .25)
  }

  if (p_mode == 1) {
    res_oi_l <- list(
      dtf %>%
        dplyr::filter(analysis_idx == 22) %>%
        dplyr::filter(analysis_name == 'twoD_sens_analysis') %>%
        dplyr::filter(LOH_HLA == 'no_LOHHLA') %>%
        dplyr::filter(overlap_var == 'mean_score') %>%
        dplyr::filter(patient_inclusion_crit == ''),
      dtf %>%
        dplyr::filter(!grepl('rooney', analysis_name)) %>%
        dplyr::filter(rc == quantile(rc, probs = .00, type = 1)),
      dtf %>%
        dplyr::filter(!grepl('rooney', analysis_name)) %>%
        dplyr::filter(p_val_signif_effect.adj == 'L') %>%
        dplyr::filter(p_val == quantile(p_val, probs = .00, type = 1)),
      dtf %>%
        dplyr::filter(!grepl('rooney', analysis_name)) %>%
        dplyr::filter(p_val_signif_effect.adj == 'L') %>%
        dplyr::filter(p_val == quantile(p_val, probs = .00, type = 1)),
      NULL
      # dtf %>%
      #   dplyr::filter(intercept == min(abs(intercept - .1), na.rm = T))
    )
  } else if (p_mode == 2) {
    ## Verify that this particular setting gives similar results for all
    ## focus alleles
    res_oi_l <- dtf %>%
      dplyr::filter(!grepl('rooney', analysis_name)) %>%
      dplyr::filter(p_val == quantile(p_val, probs = .00, type = 1))
    res_oi_l <- map(focus_hlas, ~dplyr::mutate(res_oi_l, focus_allele = .x))
  }

  plots <- purrr::imap(res_oi_l, function(res_oi, i) {
    if (is.null(res_oi)) return(NULL)
    call_args <- pan_IE_res_to_call_args(res_oi)
    test <- do.call(test_continuous_IE,
      c(list('z_normalize' = F), call_args))
    p <- wrapper_plot_PS_vs_yr(
      test_object = test,
      colour_vars = 'anchor_present_f',
      include_legend = F,
      shape_var = NULL,
      projects = as.character(res_oi$project_extended), ...)[[1]]

    if (p_mode == 1 && i == 4) {
      p <- p + scale_y_continuous(limits = c(0, 2 * res_oi$intercept))
      p <- p + scale_x_continuous()
      p <- p + ylab('')
      p <- p + xlab('')
    } else {
      p <- p + scale_y_continuous()
      p <- p + scale_x_continuous()
      p <- p + ylab(paste0(quickMHC::ppHLA(res_oi$focus_allele),
                           '\nneo-antigen yield rate'))
      p <- p + xlab(paste0(quickMHC::ppHLA(res_oi$focus_allele),
                           ' PS'))
    }
    return(p)
  })
}


relationship_between_depletion_vals <- function(dtf) {
  p1 <- dtf %>%
    ggplot(aes(x = rc, y = depletion_mean)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme_fas()
  p2 <- dtf %>%
    ggplot(aes(x = depletion_mean, y = depletion_full)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme_fas()
  patchwork::wrap_plots(list(p1, p2), tag_level = 'new')
}


get_var_summary <- function(dtf, varn = 'rc') {
  dtf %>%
    dplyr::filter(project_extended == "Pan*'-'*cancer") %>%
    pull(varn) %>%
    summary
}


merge_overviews <- function(df_A, df_B) {
  if (!is.data.frame(df_B)) return(NULL)
  analysis_id_cols <- colnames(all_cont_IE_settings) %>%
    c('project_extended') %>%
    c('analysis_idx') %>%
    { . }
  # df_A <- setDT(df_A)[, c(analysis_id_cols, 'rc'), with = F]
  # df_B <- setDT(df_B)[, c(analysis_id_cols, 'rc'), with = F]
  if ('patient_inclusion_crit' %nin% colnames(df_A) ||
      'patient_inclusion_crit' %nin% colnames(df_B)) browser()
  df_A <- setDT(df_A)
  df_B <- setDT(df_B)
  df_A$patient_inclusion_crit <- unlist(df_A$patient_inclusion_crit)
  df_B$patient_inclusion_crit <- unlist(df_B$patient_inclusion_crit)
  # merged <- merge(df_A, df_B,
  #                 by = setdiff(analysis_id_cols, 'focus_allele'),
  #                 all.x = T, all.y = T)
  # merged[, mean(focus_allele.x == focus_allele.y)]
  merged <- merge(df_A, df_B, by = analysis_id_cols, all.x = T, all.y = T)
  return(dim(merged))
  # return(cor(merged$rc.x, merged$rc.y, use = 'pairwise.complete.obs'))
  # summary(df_B$rc)
  # summary(df_A$rc)

  stopifnot(unique(df_A$analysis_idx) == unique(df_B$analysis_idx))
  merged[, .(opt_string.x, opt_string.y)]
  merged[opt_string.x != opt_string.y, .N]
  merged[opt_string.x != opt_string.y & grepl('VExp', opt_string.x), .N]
  merged[, summary(rc.y - rc.x), opt_string.x != opt_string.y]
  print(cor(merged$rc.x, merged$rc.y, use = 'pairwise.complete.obs'))
  merged[, rc.x]
  merged %>% select(analysis_id_cols) %>% unique
  missing_settings <- merged %>%
    filter(is.na(rc.y)) %>%
    select(analysis_id_cols)
  merge(Y, missing_settings[1])
  merged %>% count
  return(merged)
}


#' Print out a bunch of interesting statistics
#'
#'
summarize_pan_IE_heatmap <- function(dtf) {
  if (null_dat(dtf)) return(NULL)
  cat('RC summary:\n')
  print(summary(dtf$rc))
  print(print.mean_CI(dtf$rc))

  cat('Effective editing summary:\n')
  cat('   Mean CI_l:\n')
  print(summary(dtf$depletion_mean_ci_l))
  cat('   Mean:\n')
  print(summary(dtf$depletion_mean))
  cat('   Mean CI_h:\n')
  print(summary(dtf$depletion_mean_ci_h))
  cat('   1st_qu:\n')
  print(summary(dtf$depletion_1st_qu))
  cat('   3rd_qu:\n')
  print(summary(dtf$depletion_3rd_qu))
}


#' Plot relationship between hypothesis error and effect size for a
#' large number of semi-dependent analyses (all pertaining to one
#' particular tumor type for instance)
#'
#'
plot_bayesian_error_vs_effect <- function(dtf, x_var = '-log2_rc_error',
  return_mode = 'patchwork',
  grouping = NULL) {
  dtf[, log2_evid_ratio := log2(evid_ratio)]
  dtf[, log2_est_error := log2(est_error)]
  # point_dat <- dtf %>%
  #   dplyr::summarize(across(matches('log2_|yr_fractional'), median))

  N_bins <- min(200, 1/10 * dtf[, .N])
  dtf <- dtf[log2_est_error <= 10]

  po <- dtf[, wa_var(.SD, 'estimate'), by = project_extended] %>%
    .[order(V1), project_extended]
  dtf[, project_extended := factor(project_extended, levels = po)]

  p1 <- ggplot(mapping = aes_string(x_var, 'log2_evid_ratio'),
    data = dtf) +
    geom_hex(bins = N_bins) +
    scale_fill_viridis_c(name = '# sub-analyses') +
    ylab('log2 evidence ratio\n(>0 in agreement with immunoediting)') +
    xlab('-log2 estimation error') +
    geom_smooth(color = 'indianred3')
    # geom_point(data = point_dat, pch = 3, color = 'red', size = 3) +
    # theme_fas(
    #   legend.key.width = unit(5, 'mm'),
    #   legend.box.background = element_blank()
    # )
  if (!is.null(grouping) && grouping == 'project') {
    p1 <- p1 +
      facet_wrap(~project_extended, labeller = label_parsed)
  }

  if (return_mode == 'patchwork') {
    # dtf_s <- dtf[abs(yr_fractional_change) <= 1 & log2_est_error <= 0]
    dtf_s <- dtf
    p2 <- ggplot(mapping = aes_string(x_var, 'yr_fractional_change'),
      data = dtf_s) +
      scale_y_continuous(
        breaks = seq(-1, 1, .05),
        limits = c(-.1, .1),
        name = 'Fractional change in YR\nbetween PS = 0 and PS = 1'
      ) +
      geom_hex(bins = N_bins) +
      scale_fill_viridis_c(name = '# sub-analyses') +
      xlab('-log2 estimation error') +
      geom_smooth(color = 'indianred3') +
      # geom_smooth() +
      # geom_point(data = point_dat, pch = 3, color = 'red', size = 3) +
      theme_fas(
        legend.key.width = unit(5, 'mm'),
        legend.box.background = element_blank()
      )
    if (!is.null(grouping) && grouping == 'project') {
      p2 <- p2 +
        facet_wrap(~project_extended, labeller = label_parsed)
    }
  }

  if (return_mode == 'patchwork') {
    library(patchwork)
    return(p1 + p2)
  } else {
    return(p1)
  }
}


call_func <- function(f, args) {
  overlapping_args <- intersect(names(formals(f)), names(args))
  # ignored_args <- setdiff(names(formals(f)), names(args))
  ignored_args <- setdiff(names(args), names(formals(f)))
  if (length(ignored_args) > 0) {
    message('Following list items from call to', ' are ignored:',
      paste0(ignored_args, collapse = ', '))
  }
  purrr::exec(f, !!!args[overlapping_args])
}


pipe_row_print <- function(dtf) {
  cat('Number of rows: ', nrow(dtf), '\n')
  return(dtf)
}


friendly_factor <- function(v) {
  res <- as.character(v)
  res[is.na(res)] <- 'none'
  factor(res, levels = unique(res))
}


wa_var <- function(dtf, varname, error_var = 'est_error',
  error_trans = identity) {
  dtf[, signif(sum(get(varname)/error_trans(get(error_var))) /
    sum(1 / error_trans(get(error_var))), 3)] %>%
    setNames(varname)
}


N_comb <- function(N = 50) {
  # N = 5
  N = 9
  ## Pick 2
  N * (N-1) / 2
  ## Pick 3
  N * (N-1) * (N-2) / 6
}


#' Compare (mean) coefficient estimates of a bayesian model to that of robust
#' frequentist estimates
#'
#'
compare_bayes_to_rlm <- function(b_fit, prep) {
  if (!check_rhat(b_fit)) return(NULL)
  brm_coefs <- summary(b_fit)$fixed[, 'Estimate'] %>%
    setNames(c('intercept', 'rc'))
  rlm_fit <- unlist(rlm_fit_model(prep$dtf)[names(brm_coefs)])
  (brm_coefs - rlm_fit) / rlm_fit
}


#' Bayesian Fishtail plot
#'
#'
var_vs_error <- function(
  dtf, x = '-log2_est_error', y, binsize = NULL) {

  binsize <- binsize %||% (1/1000 * nrow(dtf))
  x_label <- switch(x,
    '-log10(scale)' = '-log10(MAD of residuals)',
    '-scale' = '-MAD of residuals',
    '-log10(norm_scale)' = '-log10(MAD of residuals)',
    '-norm_scale' = '-MAD of residuals',
    '-log2(est_error)' = '-log2(estimation error)'
  )

  p <- dtf %>%
    ggplot(mapping = aes_string(x, y)) +
    geom_hline(yintercept = 0, color = 'grey20', size = .5) +
    # geom_vline(xintercept = 2, color = 'grey20', size = .5) +
    # geom_rect(xintercept = 2, color = 'grey20', size = .5) +
    geom_hex(bins = binsize) +
    scale_fill_viridis_c() +
    xlab(x_label) +
    theme(legend.position = 'right', legend.direction = 'vertical') +
    geom_smooth(method = 'gam', col = 'indianred3', size = .5) +
    facet_wrap(~project_extended, nrow = 5,
      labeller = label_parsed)

  return(p)
}

highlight_fishtail <- function(dtf, y,
  highlight_lev = NULL, highlight_code = NULL) {
  setDT(dtf)

  if (!is.null(highlight_lev)) {
    h_dtf <- dtf
    for (i in seq_along(highlight_lev)) {
      h_dtf <- h_dtf[get(names(highlight_lev)[i]) == highlight_lev[[i]]]
    }
  }

  if (!missing(highlight_code)) {
    h_dtf <- dtf[eval(sel)]
  }

  binsize <- 1/1000 * h_dtf[, .N]

  p1 <- var_vs_error(dtf = dtf, y = y, binsize = binsize)
  p2 <- var_vs_error(dtf = h_dtf, y = y, binsize = binsize)
  library(patchwork)
  p1 + p2
}


#' Adds the single coefficient specified in row to all rows in dtf
#'
#' Computes new t and p-values for the corrected values
#'
sum_estimates <- function(dtf, row) {
  dof <- nrow(f_setting_dtf) - nrow(dtf)
  dtf %>%
    # dplyr::mutate(t_recon = estimate / std_error) %>%
    # dplyr::mutate(p_recon = 2 * pt(t_recon, dof, lower.tail = t < 0)) %>%
    dplyr::rename_with(~sapply(.x, function(x) switch(x,
        'estimate' = 'estimate_old',
        'std_error' = 'std_error_old',
        't' = 't_old',
        'p' = 'p_old',
        'conf.high' = 'CI_h_old',
        'conf.low' = 'CI_l_old',
        x
      ))
    ) %>%
    dplyr::mutate(estimate = estimate_old + row$estimate) %>%
    dplyr::mutate(std_error = std_error_old + row$std_error) %>%
    dplyr::mutate(t = estimate / std_error) %>%
    dplyr::mutate(p = 2 * ifelse(t < 0, pt(t, dof, lower.tail = T),
        1 - pt(t, dof, lower.tail = T))) %>%
    dplyr::arrange(estimate) %>%
    dplyr::mutate(term = factor(term, levels = unique(term))) %>%
    dplyr::select(term, everything()) %>%
    { . }
}


format_coefs <- format_formatted_coefs <- function(dtf) {
  setDT(dtf)
  if ('term' %in% colnames(dtf) && 'coef' %nin% colnames(dtf)) {
    dtf[, 'coef' := term]
  }

  out <- dtf %>%
    .[grepl('Intercept', coef),
      'coef_type' := 'Intercept'] %>%
    .[, coef := gsub('\\(Intercept\\)', 'Intercept', coef)] %>%
    .[grepl('tumor_type', coef),
      'coef_type' := 'Tumor type'] %>%
    .[, coef := gsub('tumor_type', '', coef)] %>%
    .[grepl('expression_threshold', coef),
      'coef_type' := 'Bulk expression'] %>%
    .[, coef := gsub('expression_threshold', '', coef)] %>%
    .[grepl('VE_threshold', coef),
      'coef_type' := 'Variant level expression'] %>%
    .[, coef := gsub('VE_threshold', '', coef)] %>%
    .[grepl('focus_allele', coef), 'coef_type' := 'Focus HLA allele'] %>%
    .[, coef := gsub('focus_allele', '', coef)] %>%
    .[grepl('sts_filtering', coef), 'coef_type' := 'STS filtering'] %>%
    .[, coef := gsub('sts_filtering_TRUE', 'STS', coef)] %>%
    # .[, coef := gsub('sts_filtering_FALSE', 'STS', coef)] %>%
    .[grepl('patient_inclusion', coef), 'coef_type' := 'Exclusion of IT resistant patients'] %>%
    .[, coef := gsub('patient_inclusion_crit', '', coef)] %>%
    .[grepl('analysis_name', coef),
     'coef_type' := 'Response variable'] %>%
    .[, coef := gsub('analysis_name', '', coef)] %>%
    .[, coef := gsub('_param_titration', '', coef)] %>%
    .[grepl('percentile_rank', coef),
     'coef_type' := 'HLA binding affinity'] %>%
    .[, coef := gsub('percentile_rank', '', coef)] %>%
    .[grepl('LOH_HLA', coef), 'coef_type' := 'LOH of HLA'] %>%
    .[, coef := gsub('LOH_HLA', '', coef)] %>%
    .[grepl('mean_score', coef), 'coef_type' := 'Presentation score computation'] %>%
    .[, coef := gsub('overlap_varmean_score_', '', coef)] %>%
    .[grepl(':', coef), 'coef_type' := 'interaction term'] %>%
    { . }

  out$coef[out$coef_type == 'Focus HLA allele'] %<>% quickMHC::ppHLA()
  out$coef[out$coef == 'sts_filteringTRUE'] <- 'STS'
  out$coef[out$coef == 'rooney'] <- 'Observed/expected (Rooney et al.)'
  out$coef[out$coef == 'marty'] <- 'Recurrent oncogenic (Marty et al.)'
  out$coef[out$coef == 'driv_ess'] <- 'No drivers/essential passengers'
  out$coef[out$coef == 'clon'] <- 'Clonal mutations only'
  out$coef[out$coef == 'strict_LOHHLA'] <- 'Fully assessable patients only'
  out$coef[out$coef == 'LOHHLA'] <- 'LOH HLA for assessable alleles'
  out$coef[out$coef == 'FDR10'] <- 'Stringent'
  out$coef[out$coef == 'FDR1'] <- 'Lenient'
  out$label <- paste(out$coef_type, '-', out$coef)
  out$label[out$label == 'Intercept - Intercept'] <- 'Intercept'
  return(out)
}


coef_volcano <- function(dtf, label_dat = dtf, mod = identity) {
  max_x <- 1.05 * setDT(dtf)[, max(abs(estimate))]

  if (!is.null(dtf$coef_type))
    color_var = 'coef_type'
  else
    color_var = NULL

  ggplot(dtf, aes_string(x = mod('estimate'), y = mod('-log10(p)'),
      label = 'coef', color = color_var)) +
    geom_point() +
    geom_vline(xintercept = 0) +
    ggrepel::geom_text_repel(data = label_dat, force = 10, show.legend = F) +
    scale_x_continuous(limits = c(-max_x, max_x)) +
    scale_colour_discrete(name = 'Coefficient type')
}


forest_plot <- function(x, ...) {
  UseMethod('forest_plot')
}


forest_plot.data.frame <- function(t_dat,
  x_var = 'project_extended', bg_var = NULL) {

  if ('V1' %in% colnames(t_dat)) {
    y <- 'V2'; ymin <- 'V1'; ymax <- 'V3'
    col_var <- NULL
  } else if ('q50' %in% colnames(t_dat)) {
    y <- '`50%`'; ymin <- '`5%`'; ymax <- '`95%`'
    y <- 'q50'; ymin <- 'q05'; ymax <- 'q95'
    col_var <- NULL
  } else if ('estimate_corr' %in% colnames(t_dat)) {
    y <- 'estimate_corr'
    ymin <- 'estimate_corr - 1.96 * std_error_corr'
    ymax <- 'estimate_corr + 1.96 * std_error_corr'
    col_var <- NULL
  } else if (all(c('estimate', 'conf.low') %in% colnames(t_dat))) {
    y <- 'estimate'
    ymin <- 'conf.low'
    ymax <- 'conf.high'
    col_var <- 'coef_type'
    # col_var <- NULL
  } else if (all(c('Q2.5', 'Est.Error') %in% colnames(t_dat))) {
    y <- 'Estimate'
    ymin <- 'Q2.5'
    ymax <- 'Q97.5'
    col_var <- 'coef_type'
    # col_var <- NULL
  }

  if ((x_var == 'project_extended' ||
      any(grepl('~', t_dat[[x_var]]))) &&
      'project_extended' %in% colnames(t_dat)) {
    t_dat[[x_var]] <- as.factor(t_dat[[x_var]])
    labels <- parse(text = levels(t_dat[[x_var]]))
  } else {
    labels <- t_dat[[x_var]]
  }

  if (!is.null(col_var)) {
    t_dat$coef_type %<>% factor(levels = unique(.))
    t_dat$coef %<>% factor(levels = unique(.))
  }

  if (!is.null(bg_var)) {
    t_dat[[x_var]] <- factor(t_dat[[x_var]], 
      levels = unique(t_dat[[x_var]]))
    t_dat[[bg_var]] <- factor(t_dat[[bg_var]], 
      levels = unique(t_dat[[bg_var]]))
  }

  p <- ggplot(t_dat,
    aes_string(x = x_var, y = y, ymin = ymin, ymax = ymax,
      color = col_var))

  if (!is.null(bg_var)) {
    fill_cols <- rev(c('gray85', 'white'))
    v <- as.integer(t_dat[[bg_var]])
    strip_dat <- suppressWarnings(
      data.table(
        'x_min' = map_int(unique(v), ~min(which(v == .x))) - .5,
        'x_max' = map_int(unique(v), ~max(which(v == .x))) + .5,
        'strip_col' = fill_cols
      )
    )
    p <- p + geom_rect(
      data = strip_dat,
      mapping = aes(xmin = x_min, xmax = x_max,
        ymin = -Inf, ymax = Inf, fill = strip_col),
      inherit.aes = F, alpha = .5
    )
    p <- p + scale_fill_manual(
      name = '', guide = F, values = auto_name(fill_cols))
  }

  p <- p +
    geom_point(pch = 20, colour = 'black') +
    geom_hline(yintercept = 0, colour = 'grey50') +
    geom_linerange() +
    scale_x_discrete(name = '', labels = labels) +
    coord_flip()

  if (F && !is.null(col_var)) {
    p <- p + scale_color_brewer(
      palette = 'Set2',
      name = 'Coefficient type'
    )
  }
  return(p)
}

delta_label <-
  expression(Delta*'  between  '*italic(h)==0~' and '*italic(h)==1)

forest_plot.setting_titration <- function(t_dat, x_var = 'project_extended') {
  if (x_var == 'project_extended' ||
      any(grepl('~', t_dat[[x_var]]))) {
    t_dat[[x_var]] <- as.factor(t_dat[[x_var]])
    labels <- parse(text = levels(t_dat[[x_var]]))
  } else {
    labels <- t_dat[[x_var]]
  }

  pos <- position_dodge(width = .5)

  ## coord_flip doesn't play nice with coord_cartesian, so cap values
  ## ourselves
  range <- 1.10 * range(t_dat$mean)
  num_cols <- map_lgl(t_dat, is.numeric) %>% { names(.)[.] }
  cap <- function(v, range) {
    v[v <= range[1]] <- range[1]
    v[v >= range[2]] <- range[2]
    return(v)
  }
  for (cn in num_cols) {
    t_dat[[cn]] <- cap(t_dat[[cn]], range)
  }

  p <- ggplot(t_dat, aes_string(x = x_var, y = 'mean',
      ymin = 'mean_l', ymax = 'mean_h',
      # fill = 'version',
      color = 'version'))

  if (T) {
    fill_cols <- rev(c('gray85', 'white'))
    strip_dat <- suppressWarnings(
      data.table(
        'x_min' = seq(0, uniqueN(t_dat[[x_var]])) - .5,
        'x_max' = seq(1, uniqueN(t_dat[[x_var]]) + 1) - .5,
        'strip_col' = fill_cols
      )
    )
    p <- p + geom_rect(
      data = strip_dat,
      mapping = aes(xmin = x_min, xmax = x_max,
        ymin = -Inf, ymax = Inf, fill = strip_col),
      inherit.aes = F, alpha = .5
    )
    p <- p + scale_fill_manual(
      name = '', guide = F, values = auto_name(fill_cols))
  }

  p <- p +
    geom_hline(yintercept = 0, colour = 'grey50', linetype = 2) +
    geom_linerange(size = 1, position = pos) +
    geom_linerange(
      mapping = aes(ymin = pred_l, ymax = pred_h),
      inherit.aes = T,
      alpha = .5, size = .5, position = pos) +
    geom_linerange(
      mapping = aes(ymin = q25, ymax = q75),
      inherit.aes = T,
      alpha = .5, size = 1, position = pos) +
    # geom_point(position = pos, size = 1,
    #   shape = 21, colour = 'black') +
    geom_point(position = pos, size = 1.05, color = 'black',
      show.guide = F) +
    geom_point(position = pos, size = 1) +
    scale_y_continuous(
      name = delta_label,
      expand = c(0, 0)
    ) +
    scale_x_discrete(name = '', labels = labels, expand = c(0, 0)) +
    # coord_cartesian(ylim = ) +
    coord_flip() +
    scale_color_brewer(
      palette = 'Set2',
      name = 'Analysis version\n(incremental)') +
    # coord_cartesian(ylim = range(t_dat$mean)) +
    theme(legend.position = 'bottom',
      legend.justification = c(.5, .5))
}

summarize_by_project <- function(dtf) {
  setDT(dtf)
  if ('ci_lower' %in% colnames(dtf)) {
    t_dat <- dtf %>%
      pipe_row_print %>%
      .[, .(
        wa_var(.SD, 'ci_lower'),
        wa_var(.SD, 'estimate'),
        wa_var(.SD, 'ci_upper')), by = project_extended] %>%
    .[order(V2)] %>%
      { . }
  } else if ('norm_scale' %in% colnames(dtf)) {
    # pacman::p_load('reldist')
    t_dat <- dtf[, {
      # ED <- spatstat::ewcdf(
      #   x = yr_fractional_change,
      #   weights = 1 - norm_scale)
      # ED(c(.05, .5, .95))
      reldist::wtd.quantile(
        x = yr_fractional_change,
        q = c(.05, .5, .95), na.rm = FALSE,
        weight = 1 - norm_scale
      ) %>% set_names(c('q05', 'q50', 'q95')) %>% as.list()
    }, by = project_extended]
  }

  t_dat %>%
    dplyr::mutate(
      tumor_type = tumor_types[as.character(project_extended)]) %>%
    dplyr::mutate(
      tumor_type = factor(tumor_type, levels = unique(tumor_type))) %>%
    { . }
}


add_frac <- function(dtf) {
  dtf[, 'frac' := N / sum(N)]
}


# prep <- call_func(prep_cont_IE_analyses, args)
# b_g_lm <- fit_bayesian_regression(prep$dtf, prior = 'biased',
#   return_val = 'lm')
# plot_id_base <- paste(args[analysis_grp_vars], collapse = '-')
# p <- plot_PS_vs_yr(dtf = prep$dtf, reg_method = 'bayesian_biased', lm = b_g_lm)


plot_hypothesis <- function(editing_test, label_size = 6) {
  x_range <- with(editing_test$hypothesis,
    c(CI.Lower, CI.Upper) + .1 * abs(CI.Lower - CI.Upper) * c(-1, 1))
  x_range <- max(abs(x_range)) * c(-1, 1)
  p_hypo <- plot(editing_test, plot = F)[[1]] +
    xlim(x_range[1], x_range[2]) +
    # ggtitle(expression('Prior and posterior of '~Delta)) +
    ggtitle('') +
    labs(title = '', subtitle = '') +
    annotate_npc(
      label = glue('D = ', signif(editing_test$hypothesis$Estimate,
          digits = 3)),
      x = .95, y = .95, hjust = 1, vjust = 1, gp = gpar(fontsize = label_size)
    ) +
    annotate_npc(
      label = glue('SE = ',
        signif(-log2(editing_test$hypothesis$Est.Error), digits = 3)),
      x = .95, y = .85, hjust = 1, vjust = 1, gp = gpar(fontsize = label_size)
    ) +
    annotate_npc(
      label = glue('Evidence ratio = ', signif(editing_test$hypothesis$Evid.Ratio, digits = 3)),
      x = .95, y = .75, hjust = 1, vjust = 1, gp = gpar(fontsize = label_size)
    ) +
    theme_fas()
}


#' Plot low-level results for a single sub-anals
#'
#' @param pick A row from a coef overview object. Irrelevant columns will be
#' discarded
#' @param idx A row index from coef_overview sample_subs
#'
inspect_coefs <- compare_methods <- function(idx = NULL, pick = NULL) {
  source(file.path(ma_dir, 'immune_editing', 'load_brms_models.R'))
  if (is.null(idx)) {
    stopifnot(!is.null(pick))
  } else {
    sample_subs <- setting_dtf[order(yr_fractional_change)] %>%
      .[log2_est_error < 0]
    pick <- sample_subs[idx]
  }
  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    {
      modifyList(., list(redo = F, tumor_type = pick$project_extended,
          'z_normalize' = F))
    }

  prep <- call_func(prep_cont_IE_analyses, args)

  av_models <- ls(brm_scaffold) %>%
    intersect('biased_fit') %>%
    { . }
  # source(file.path('~/antigenic_space', 'maarten-analyses', 'immune_editing',
  #                  'continuous_IE_detection_init.R'))
  bayesian_mods <- av_models %>%
    auto_name %>%
    purrr::map(function(mn) {
      res <- fit_bayesian_regression(prep$dtf, model_name = mn,
        N_iter = 8000)
      return(res)
    })

  p1 <- qplot(ol, y_var, data = prep$dtf) +
    geom_smooth(aes(colour = 'rlm', fill = 'rlm'), method = MASS::rlm) +
    # geom_smooth(aes(colour = 'loess', fill = 'loess'), method = 'loess') +
    geom_smooth(aes(colour = 'gam', fill = 'gam'), method = 'gam') +
    geom_smooth(aes(colour = 'lm', fill = 'lm'), method = 'lm') +
    scale_colour_discrete(name = 'Model type') +
    ylim(0, max(prep$dtf$y_var)) +
    guides(fill = 'none') +
    theme_fas()

  summary(bayesian_mods[[1]])$fixed['Intercept', 1]
  # source(file.path('~/antigenic_space', 'maarten-analyses', 'immune_editing',
  #                  'continuous_IE_detection_init.R'))
  bayesian_plots <-
    purrr::map(bayesian_mods, function(m) {
      plot_PS_vs_yr(dtf = prep$dtf, lm = m, reg_method = 'bayesian',
        shape_var = NULL, colour_vars = NULL) +
        ylim(c(0, summary(bayesian_mods[[1]])$fixed['Intercept', 1] * 2))
    })


  if (length(bayesian_mods) > 1) {
    loo_mods <- map(bayesian_mods, ~add_criterion(.x, 'loo')) %>%
      { . }

    # test_plot(print(bayesian_plots[[2]]))
    # res <- purrr::exec(brms::loo_compare, !!!unname(loo_mods))
    res <- brms::loo_compare(loo_mods[[1]], loo_mods[[2]], loo_mods[[3]])
    ## Ugly code. Extract index to retrieve model ordering
    dimnames(res)[[1]] <-
      gsub('.*\\[\\[(\\d+)\\]\\]$', '\\1', dimnames(res)[[1]]) %>%
      as.integer %>%
      { names(bayesian_mods)[.] }
  }

  fit <- bayesian_mods[['biased_fit']]
  editing_test <- tryCatch(suppressWarnings({
    brms::hypothesis(x = fit, hypothesis = '(ol / Intercept) < 0')
  }), error = function(e) { print(e); browser() })
  p_hypo <- plot_hypothesis(editing_test)

  editing_test <- tryCatch(suppressWarnings({
    brms::hypothesis(x = fit, hypothesis = '(ol) < 0')
  }), error = function(e) { print(e); browser() })
  p_hypo_alt <- plot_hypothesis(editing_test)

  ## Create filename from 'settings'
  plot_id <- args %>%
    modifyList(list(
        include_call = NULL,
        redo = NULL,
        z_normalize = NULL)) %>%
    imap(~glue('{.y}={.x}')) %>%
    paste(collapse = '-') %>%
    { gsub('\'|\\*', '', .) }

  plots <- bayesian_plots %>% append(list(p1, p_hypo, p_hypo_alt))
  img_fn <- file.path(img_loc,
    glue('bayesian_coef_panel-{plot_id}.pdf'))
  plot_panel_layout(plots = plots, filename = img_fn,
    nrow = 3, ncol =2)

  img_fn <- file.path(img_loc, glue('bayesian_pars-{plot_id}.png'))
  print_plot(plot(fit), fn = img_fn, w = 17.4, h = 20)
}


compare_methods_old <- function(
  idx = sample(1:nrow(sample_subs), 1),
  pick = NULL) {

  if (is.null(idx)) {
    stopifnot(!is.null(pick))
  } else {
    pick <- sample_subs[idx]
  }

  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    append(c('z_normalize' = F)) %>%
    # modifyList(list(redo = F, project_extended = pe))
    modifyList(list(redo = F))

  if (T) {
    ## This is a bit wasteful, we're computing regression for all tumor types,
    ## not just the one we're interested in right now. Oh well
    bayesian_res <- purrr::map_dfr(
      c('bayesian_biased_robust', 'bayesian_unbiased', 'bayesian_biased'),
      function(rm) {
        test <- rlang::exec(test_continuous_IE, !!!args, reg_method = rm)
        l <- test$stats[[pe]][[1]]
        return(l)
      })
  } else {
  }

  select_vars <- c('intercept', 'rc', 'p_val')
  rlm_test <- tryCatch({
    rlang::exec(test_continuous_IE, !!!args,
      reg_method = 'rlm_one_group')$stats[[pe]][[1]][select_vars]
    }, error = function(e) {
      print(e); return(map(auto_name(select_vars), ~NA))
    })

  if (!is.null(rlm_test$rc) && !is.na(rlm_test$rc)) {
    rlm_res <- tibble(
      rlm_intercept = rlm_test$intercept,
      rlm_rc = rlm_test$rc,
      intercept_change = (bayesian_res$intercept_estimate -
        rlm_test$intercept) / rlm_test$intercept,
      om_change = (bayesian_res$ol_estimate - rlm_test$intercept) /
        rlm_test$intercept,
      intercept_change_abs = bayesian_res$intercept_estimate -
        rlm_test$intercept,
      om_change_abs = bayesian_res$ol_estimate - rlm_test$intercept
    )
    res <- tibble(idx = idx, bayesian_res, rlm_res)
  } else {
    res <- tibble(idx = idx, bayesian_res)
  }
  return(res)
}


annotate_npc <- function(label, x, y, ...) {
  ggplot2::annotation_custom(grid::textGrob(
    x = unit(x, "npc"), y = unit(y, "npc"), label = label, ...))
}


lower_hex <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_abline(intercept = 0, slope = 1, color = 'grey10',
      linetype = 1) +
    geom_vline(xintercept = 0, color = 'grey80', linetype = 1) +
    geom_hline(yintercept = 0, color = 'grey80', linetype = 1) +
    geom_hex(bins = nrow(data) / 1000, ...) +
    scale_fill_viridis_c()
}


diag_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_vline(xintercept = 0, color = 'grey10') +
    geom_vline(xintercept = 0, color = 'grey80', linetype = 1) +
    geom_hline(yintercept = 0, color = 'grey80', linetype = 1) +
    geom_histogram(bins = 100)
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


#' Format the result of compile_all_coef_overview in preparation for factor
#' analysis
#'
#'
format_coef_overview <- function(
  dtf,
  delta_filter_error = NULL,
  filter_intercept = F,
  rc_filter_error = 2,
  pick_allele_method = 'none',
  intercept_filter_error = 2,
  plot_fishtails = NULL,
  pp_hl_threshold = .99) {

  dtf[, intercept := intercept_estimate]
  dtf[, rc := ol_estimate]
  dtf[, 'yr_fractional_change' := estimate]
  dtf[, log2_evid_ratio := log2(evid_ratio)]
  dtf[, log2_est_error := log2(est_error)]
  dtf[, log2_rc_error := log2(ol_est_error)]
  dtf[, log2_intercept_cov :=
    log2(intercept_est_error) - log2(abs(intercept_estimate))]
  dtf[, 'log2_rc_cov' := log2(ol_est_error / abs(ol_estimate))]
  dtf[post_prob <= pp_hl_threshold, star := 'G']
  dtf[post_prob >= (1-pp_hl_threshold), star := 'L']

  dtf <- filter_bayesian_coef_overview(
    dtf = dtf,
    filter_intercept = filter_intercept,
    delta_filter_error = delta_filter_error,
    rc_filter_error = rc_filter_error,
    intercept_filter_error = intercept_filter_error,
    plot_fishtails = plot_fishtails,
    pick_allele_method = pick_allele_method
  )

  stopifnot(all(dtf$sampling_convergence))
  dtf <- cbind(dtf, ds_param_grid[dtf[, analysis_idx], ])

  dtf$VE_threshold %<>% friendly_factor
  dtf$expression_threshold %<>% friendly_factor
  dtf$sts_filtering %<>% friendly_factor
  dtf$percentile_rank %<>% friendly_factor
  dtf$processing_threshold %<>% friendly_factor

  allowed_tumor_types <-
    dtf[log2_est_error <= -2, .N, by = project_extended][N > 50] %>%
    pull(project_extended)

  excluded_tumor_types <- dtf[, setdiff(project_extended, allowed_tumor_types)]
  dtf <- dtf[project_extended %in% allowed_tumor_types]

  po <- dtf[, wa_var(.SD, 'estimate'), by = project_extended] %>%
    .[order(V1), project_extended]
  dtf[, project_extended := factor(project_extended, levels = po)]
  po <- tumor_types[po]
  dtf[, tumor_type := factor(tumor_types[as.character(project_extended)],
    levels = po)]
  return(dtf)
}



sim_data <- function(
  prep, rc = 0, noise_bound = .1, lower_bound = F,
  scale_intercept = NULL, beta_x = F, ...) {

  dots <- list(...)

  b0 <- coef(lm(y_var ~ ol, prep$dtf))[1]
  if (!is.null(scale_intercept)) {
    b0 <- scale_intercept * b0
  }

  if (beta_x) {
    # ol <- rbeta(length(ol), 1, 1)
    prep$dtf$ol <- runif(length(prep$dtf$ol), 0, 1)
  }

  s_dtf <- with(prep$dtf, {
    noise <- runif(length(ol), min = -noise_bound, max = noise_bound)
    data.frame('ol' = ol, 'y_var' = (b0 + ol * rc) + noise)
  })

  if (lower_bound) {
    s_dtf$y_var <- pmax(s_dtf$y_var, 0)
  }
  return(s_dtf)
}


flatten_response <- function(dtf = NULL, prep = NULL) {
  if (is.null(dtf) && !is.null(prep$dtf)) {
    dtf <- prep$dtf
  }
  setDT(dtf)
  prelim_mod <- fit_rlm_model(dtf)
  correction <- predict(prelim_mod$lm) - prelim_mod$intercept
  dtf$y_var <- dtf$y_var - correction
  check_mod <- fit_rlm_model(dtf)
  stopifnot(maartenutils::eps(check_mod$yr_fractional_change, 1e-6))
  return(dtf)
}


#' Scale response variable of input data such that b0 will be 1
#'
#'
normalize_y_var <- function(dtf = NULL, prep = NULL) {
  if (is.null(dtf) && !is.null(prep$dtf)) {
    dtf <- prep$dtf
  }
  setDT(dtf)
  lm_mod <- fit_rlm_(dtf)
  dtf$y_var <- dtf$y_var / lm_mod$intercept
  check_mod <- fit_rlm_model(dtf)
  stopifnot(maartenutils::eps(check_mod$intercept, 1, 1e-6))
  return(dtf)
}


simulate_random_data <- function(b0 = .05, b1 = 0, eps = .05,
  zero_inflate = F, N_patients = 1000) {
  prep <- list(
    'dtf' = setDT(tibble(
      ol = runif(N_patients, 0, 1),
      y_var = b0 + b1 * ol + runif(N_patients, -eps, eps)
    ))
  )
  if (zero_inflate) {
    prep$dtf$y_var <- pmax(prep$dtf$y_var, 0)
  }
  return(prep)
}


rlm_plot_pick <- function(pick,
  label_pos_y = 'top',
  label_pos_x = 'left') {

  stopifnot(nrow(pick) == 1)
  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    modifyList(
      list(
        'z_normalize' = F,
        redo = F,
        project_extended = as.character(.$project_extended),
        tumor_type = as.character(.$project_extended)
      )
    ) %>%
    map(~if (class(.x) == 'factor') as.character(.x) else .x)

  ## Prepare and compile data
  prep <- call_func(prep_cont_IE_analyses, args)

  if (label_pos_y == 'top') {
    y_s <- c(.95, .85)
  } else if (label_pos_y == 'bottom') {
    y_s <- c(.15, .10)
  }

  if (label_pos_x == 'left') {
    x_s <- .05
    hjust <- 0
  } else if (label_pos_x == 'right') {
    x_s <- .95
    hjust <- 1
  }

  source(file.path('~/antigenic_space', 'maarten-analyses',
      'immune_editing', 'continuous_IE_detection_init.R'))
  mod <- fit_rlm_model(prep$dtf)$lm
  obs_yr <- coef(mod) %>% { .[2] / .[1] }
  if (F && !maartenutils::eps(pick$yr_fractional_change, obs_yr, 1e-3)) {
    args_wp <- args %>%
      modifyList(list(tumor_type = NULL, 
          # analysis_idx = NULL,
          project_extended = NULL))
    args_wp <- args %>%
      modifyList(list(#tumor_type = NULL, 
          #project_extended = NULL,
          analysis_idx = NULL))
    # test <- call_func(test_continuous_IE, args_wp)
    # test <- call_func(prep_continuous_param_grid, args_wp)

    # print(pick$yr_fractional_change)
    # print(obs_yr)
    # table(prep$dtf$project_extended)
    source(file.path('~/antigenic_space', 'maarten-analyses',
        'immune_editing', 'continuous_IE_detection_init.R'))
    prep <- call_func(prep_cont_IE_analyses, args)

    grid <- call_func(prep_continuous_param_grid,
      c(args, list(analysis_idxs = pick$analysis_idx)))
    grid <- call_func(prep_continuous_param_grid, args)
    grid$pipeline_param_titration_grid_dat %>%
      { .[project_extended == pick$project_extended,
          yr_fractional_change] }

    test <- call_func(test_continuous_IE, args)
    test$stats[[as.character(pick$project_extended)]][[1]] %>%
      { .$yr_fractional_change }
    print(nrow(prep$dtf) - pick$n_patients)

    table(prep$dtf$hla_allele_status_b)
    sum(is.na(prep$dtf$hla_allele_status_b))
    # mean(abs(residuals(mod))) / coef(mod)[1]
    # mod$scale
    # mad(residuals(mod))
  }

  p1 <- plot_PS_vs_yr(
      dtf = prep$dtf,
      plot_density = F,
      lm = mod,
      colour_vars = NULL,
      shape_var = NULL
    ) +
    annotate_npc(
      x = x_s, y = y_s[1],
      label = parse(text = sprintf('Delta==%.2f',
          pick$yr_fractional_change)),
      hjust = hjust, vjust = 1,
      gp = gpar(fontsize = 6)
    ) +
    annotate_npc(
      x = x_s, y = y_s[2],
      label = glue('MAD of norm. residuals',
        '= {signif(pick$norm_scale, 3)}'),
      hjust = hjust, vjust = 1,
      gp = gpar(fontsize = 6)
    )

  return(p1)
}


annotate_pipeline_settings <- function(dtf) {
  t_dat <- cbind(dtf, ds_param_grid[dtf[, analysis_idx], ])
  t_dat$VE_threshold %<>% friendly_factor
  t_dat$VE_threshold <- relevel(t_dat$VE_threshold, 'none')
  t_dat$expression_threshold %<>% friendly_factor
  t_dat$expression_threshold <-
    relevel(t_dat$expression_threshold, 'none')
  t_dat$sts_filtering %<>% friendly_factor
  t_dat$percentile_rank %<>% friendly_factor
  t_dat$percentile_rank <-
    factor(t_dat$percentile_rank, levels = c(4, 3, 1.9, 1))
  t_dat$processing_threshold %<>% friendly_factor
  t_dat$focus_allele %<>% friendly_factor
  t_dat$patient_inclusion_crit %<>% friendly_factor
  t_dat$patient_inclusion_crit <-
    factor(t_dat$patient_inclusion_crit,
      levels = c('', 'strict_TR', 'TR'),
      labels = c('none', 'FDR1', 'FDR10'))
  t_dat$LOH_HLA %<>% friendly_factor
  t_dat$LOH_HLA <- relevel(t_dat$LOH_HLA, 'no_LOHHLA')
  t_dat$overlap_var %<>% friendly_factor
  # browser(expr = any(is.na(t_dat)))
  # t_dat[which(is.na(t_dat), arr.ind = T)[1]]
  return(t_dat)
}


extract_levels <- function(v, start_i = 5) {
  if (length(v) == 1) {
    factor_name <- switch(v,
      'sts_filteringTRUE' = 'STS',
      '(Intercept)' = 'tumor_type',
      'overlap_varmean_score_AB' = 'overlap_HLAs'
    )
    factor_level <- switch(v,
      'sts_filteringTRUE' = 'STS',
      '(Intercept)' = '',
      'overlap_varmean_score_AB' = 'AB'
    )
    return(list(factor_level) %>% set_names(factor_name))
  }
  str_lengths <- stringr::str_length(v)
  min_length <- min(str_lengths)
  i <- start_i
  i <- 1
  while(i <= min_length && uniqueN(substr(v, 1, i)) == 1) {
    i <- i + 1
  }
  factor_name <- unique(substr(v, 1, i-1))
  factor_levels <- substr(v, i, str_lengths)
  print(factor_name)
  print(factor_levels)
  list(factor_levels) %>%
    set_names(factor_name) %>%
    return
}


fit_meta_lm <- function(dtf, ds_N = NULL, method = 'rem',
  grp_var = 'tumor_type',
  FE_vars = c('tumor_type', 'focus_allele', 'VE_threshold',
    'patient_inclusion_crit', 'expression_threshold', 'LOH_HLA',
    'analysis_name', 'sts_filtering', 'percentile_rank',
    'overlap_var'),
  RE_vars = c('patient_inclusion_crit', 'focus_allele')) {

  if (!is.null(ds_N)) {
    dtf <- dtf[sample(seq(1, .N), ds_N)]
  }

  allowed_vars <-
    map_lgl(FE_vars, function(v) dtf[, uniqueN(get(v))] > 1)
  FE_vars <- FE_vars[allowed_vars]
  if (length(FE_vars) == 0) return(NULL)

  FE_vars <- setdiff(FE_vars, grp_var)
  RE_vars <- setdiff(RE_vars, grp_var)
  # FE_vars <- setdiff(FE_vars, RE_vars)

  if (grepl('lm', method)) {
    if (method == 'interaction_lm') {
      c_formula <-
        sprintf('yr_fractional_change ~ (%s)^2',
          paste(FE_vars, collapse = ' + ')) %>%
        as.formula()
    } else if (method == 'lm') {
      c_formula <-
        sprintf('yr_fractional_change ~ (%s)',
          paste(FE_vars, collapse = ' + ')) %>%
        as.formula()
    }
    mod <- tryCatch(
      # lm(c_formula, data = dtf, weights = norm_scale),
      lm(c_formula, data = dtf),
      error = function(e) { print(e); NULL })
  } else if (method == 'rem') {
    library(lme4)
    library(lmerTest)

    c_formula <-
      sprintf('yr_fractional_change ~ %s + (1 + %s | %s)',
        paste(FE_vars, collapse = ' + '),
        paste(RE_vars, collapse = ' + '),
        grp_var) %>%
      as.formula()

    t_dat <- dtf
    time <- system.time(mod <- lme4::lmer(
      formula = c_formula,
      data = t_dat,
      REML = T,
      control = lme4::lmerControl(
        optimizer = 'nloptwrap'
        , optCtrl = list(
            maxeval = 10000
            # , xtol_abs = 5e-1
            # , xtol_rel = 5e-1
            # , xtol_abs = 1e-6
            # , xtol_rel = 5e-3
          )
      ),
      # weights = 1 / (t_dat$norm_scale + 1e-4)
      weights = t_dat$norm_scale + 1e-4
    ))
  } else if (grepl('brem', method)) {
    c_formula <-
      sprintf('yr_fractional_change|se(norm_scale) ~ %s + (%s | %s)',
        paste(FE_vars, collapse = ' + '),
        paste(RE_vars, collapse = ' + '),
        grp_var) %>%
      as.formula()

    time <- system.time(mod <- brms::brm(c_formula, data = dtf,
      prior = c(
        prior_string('normal(0, .2)', class = 'Intercept'),
        prior_string('normal(0, .5)', class = 'sd')
      ),
      cores = 16, chains = 16,
      sample_prior = 'yes', iter = 1000))
    attr(mod, 'time') <- time
    print(time)
  } else if (method == 'ridge') {
    allowed_vars <-
      map_lgl(FE_vars, function(v) dtf[, uniqueN(get(v))] > 1)
    FE_vars <- FE_vars[allowed_vars]
    if (length(FE_vars) == 0) return(NULL)

    c_formula <-
      sprintf('yr_fractional_change ~ %s',
        paste(FE_vars, collapse = ' + ')) %>%
      as.formula()
    c_formula <-
      sprintf('yr_fractional_change ~ 0 + (%s)^2',
        paste(FE_vars, collapse = ' + ')) %>%
      as.formula()
    MM <- model.matrix(c_formula, data = dtf)
    library(glmnet)
    mod <- cv.glmnet(
      y = dtf$yr_fractional_change,
      x = MM,
      alpha = 0,
      weights = 1 - dtf$norm_scale,
      family = "gaussian",
      lambda.min.ratio = 1e-6)
    # plot(model)
    attr(mod, 'MM') <- MM
  }
  alarm()
  return(mod)
}


#' Summarize and sum model coefficients
#'
#' The summing that's going on here is probably unsound; it doesn't
#' account for parameter correlation at all but rather assumes terms
#' to be independent of each other
#'
extract_tumor_type_estimates <- function(lm_mod) {
  if (is.null(lm_mod)) return(NULL)
  library(broom)
  library(broom.mixed)
  coef_df <- tidy(lm_mod, conf.int = T) %>%
    dplyr::rename(std_error = std.error, t = statistic)
  # , p = p.value
  intercept <- coef_df[coef_df$term == '(Intercept)', ]
  exp_VE <- coef_df[coef_df$term == 'VE_threshold5', ]
  LOHHLA <- coef_df[coef_df$term == 'LOH_HLAstrict_LOHHLA', ]
  sts <- coef_df[coef_df$term == 'sts_filteringTRUE', ]
  percentile <- coef_df[coef_df$term == 'percentile_rank1', ]
  patient_inclusion <- coef_df[coef_df$term == 'patient_inclusion_critFDR10', ]

  regex = '^tumor_type[^:]*$'

  formatted_coefs <- coef_df %>%
    dplyr::filter(grepl(regex, term)) %>%
    dplyr::mutate(term = gsub('tumor_type', '', term)) %>%
    dplyr::mutate(term = tumor_types_inv[term]) %>%
    dplyr::filter(!is.na(term)) %>%
    sum_estimates(intercept) %>%
    # debug_pipe %>%
    sum_estimates(exp_VE) %>%
    sum_estimates(LOHHLA) %>%
    sum_estimates(sts) %>%
    sum_estimates(percentile) %>%
    sum_estimates(patient_inclusion_crit) %>%
    { . }

  return(formatted_coefs)
}


change_factor_value <- function(fac, lev) {
  factor(lev, levels = levels(fac))
}


run_lmer_preds <- function(pick, me_mod, version = '',
  N_sims = 1e4) {
  predictions <- merTools::predictInterval(
    merMod = me_mod,
    which = 'full',
    # which = 'fixed',
    newdata = pick,
    returnSims = T,
    stat = 'median',
    .parallel = F,
    type = 'linear',
    n.sims = N_sims)

  pacman::p_load('qwraps2')
  means <- attr(predictions, 'sim.results') %>%
    apply(1, qwraps2::mean_ci)
  quartiles <- attr(predictions, 'sim.results') %>%
    apply(1, quantile, probs = c(.25, .75))
  perc_negative <- attr(predictions, 'sim.results') %>%
    apply(1, function(x) mean(x < 0))
  out <- cbind('version' = version, predictions, t(means),
    t(quartiles), perc_negative = perc_negative) %>%
    set_names(c('version', 'pred_median', 'pred_l', 'pred_h', 'mean',
        'mean_l', 'mean_h', 'q25', 'q75', 'perc_negative'))

  return(out)
}

top_left_legend <-
  theme(legend.position = c(0.05, .95), legend.justification = c(0, 1))

top_right_legend <-
  theme(legend.position = c(0.95, .95), legend.justification = c(1, 1))

right_legend <-
  theme(legend.position = 'right', legend.justification = c(.5, .5))

no_legend <-
  theme(legend.position = 'none')
