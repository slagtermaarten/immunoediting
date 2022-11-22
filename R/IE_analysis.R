theme_set(
  theme_fas(
    legend.box.background = element_blank(),
    legend.key.width = unit(5, 'mm'),
    panel.spacing = unit(.4, 'lines'),
    base_size = 6,
    legend.direction = 'vertical',
    legend.position = c(.95, .95),
    legend.justification = c(1, 1),
    panel.background = element_rect('grey93', colour = NULL,
      size = NULL, linetype = NULL, color = NULL),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
)

CV_IE_settings <- c('focus_allele',
  'overlap_var', 'patient_inclusion_crit', 'LOH_HLA',
  'processing_threshold', 'sts_filtering', 'percentile_rank',
  # 'expression_threshold', 'VE_threshold',
  'RNA_expression',
  'analysis_name')

IE_settings_pal <-
  gen_color_vector(arg = CV_IE_settings, name = 'FantasticFox1')

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
    breaks = analysis_names[c(1, 2, 3, 4, 5)],
    labels = c('All variants', 'Rooney O/E', 'Clonal',
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


var_names <- list(
  'Variant\nselection' = 'analysis_name',
  'Focus HLA-allele' = 'focus_allele',
  'Presentation\nscore' = 'overlap_var',
  'IT-resistant\npatients' = 'patient_inclusion_crit',
  'Adapting PS\nto LOH in HLA' = 'LOH_HLA',
  'Similarity\nto self-filter' = 'sts_filtering',
  'Affinity rank\npercentile threshold' = 'percentile_rank',
  'Bulk expression' = 'expression_threshold',
  'Variant level expression' = 'VE_threshold',
  'Expression\nfiltering' = 'exp_filter') %>%
  { as.list(set_names(names(.), .)) }
var_names[['b0']] <- var_names[['intercept']] <- 'hat(b[0])'
var_names[['b1']] <- var_names[['rc']] <- 'hat(b[1])'
var_names[['delta_mean']] <-
  var_names[['yr_fractional_change']] <- 'hat(Delta)'
var_names[['delta_SE']] <- 'SE(hat(Delta))'
var_names[['scale']] <- 'nMAE'

IE_labels <- list(
  '-log10(scale)' = '-log10(MAD of residuals)',
  '-scale' = '-MAD of residuals',
  'scale' = 'nMAE',
  '-log10(norm_scale)' = '-log10(MAD of residuals)',
  '-norm_scale' = '-MAD of residuals',
  '-log2(est_error)' = '-log2(estimation error)',
  'NMAD_prob_5' = expression(P(NMADR <= .5)),
  'NMAD_q75' = expression(P(NMADR <= .5)),
  'delta_mean' = expression(Delta*'  between  '*italic(h)==0~' and '*italic(h)==1),
  'intercept' = expression(hat(beta[0])),
  'log10_n_patients' = 'log10 patient number',
  'perm_delta_mean_c_mean' =
    expression('Patient permutation-corrected  '*Delta),
  'expression_threshold' = 'Gene-level expression threshold',
  'VE_threshold' = 'Variant-level expression threshold'
)


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
pan_IE_res_to_call_args <- function(pick) {
  ## Complement with all other required arguments to
  ## test_continuous_IE
  out <- as.list(pick) %>%
    { .[intersect(names(.), names(formals(test_continuous_IE)))] } %>%
    c(list('include_call' = T, 'redo' = F)) %>%
    map(~if (class(.x) == 'factor') as.character(.x) else .x)

  if (is.na(out$patient_inclusion_crit)) {
    out$patient_inclusion_crit <- 'FDR1'
  }

  return(out)
}


#' The fundamental and elementary visualisation of our 'continuous'
#' IE-detection method
#'
plot_PS_vs_yr <- function(
  dtf, lm,
  reg_method = 'rlm_one_group',
  # colour_vars = 'anchor_present_f',
  colour_vars = 'weight',
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
    dtf$weight <- lm$w
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
    } else if (colour_vars %in%
      c('anchor_present_f', 'hla_allele_status')) {
      col_title <- 'Patient carries\nfocus allele'
    } else if (colour_vars == 'weight') {
      col_title <- 'Patient\nregression weight'
    }
  }

  p <- ggplot(scatter_dat, aes_string(x = 'ol', y = 'y_var',
      colour = colour_vars, shape = shape_var))

  if (!plot_density) {
    p <- p + geom_point(alpha = point_alpha)
    if (!is.null(colour_vars) && colour_vars != 'weight') {
      l_cols <- pipeline_evaluation_cols %>%
        { .[1:(length(unique(dtf[[colour_vars]])))] } %>%
        setNames(unique(dtf[[colour_vars]]))
      l_shapes <- setNames(c(15:17),
        as.character(unique(dtf[[colour_vars]])))
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
        mapping = aes_string(
          ymin = 'y_var_CI_l', ymax = 'y_var_CI_h',
          size = NULL, colour = NULL, shape = NULL
        ), fill = 'firebrick3', show.legend = F, alpha = .2) +
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
    res_oi_l <- map(focus_hlas,
      ~dplyr::mutate(res_oi_l, focus_allele = .x))
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
plot_bayesian_error_vs_effect <- function(
  dtf,
  x_var = '-log2_rc_error',
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
  if (missing(args) || is.null(args) || !is.list(args))
    stop('No valid args')
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

  x_label <- if (x %in% names(IE_labels)) IE_labels[[x]] else x
  y_label <- if (y %in% names(IE_labels)) IE_labels[[y]] else y

  p <- dtf %>%
    ggplot(mapping = aes_string(x, y)) +
    geom_hline(yintercept = 0, color = 'grey20', size = .5) +
    # geom_vline(xintercept = 2, color = 'grey20', size = .5) +
    # geom_rect(xintercept = 2, color = 'grey20', size = .5) +
    geom_hex(bins = binsize) +
    scale_fill_viridis_c() +
    xlab(x_label) +
    ylab(y_label) +
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
    # dplyr::mutate(p_recon = 2 * pt(t_recon, dof,
    # lower.tail = t < 0)) %>%
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
  if (null_dat(dtf)) return(NULL)
  if ('term' %in% colnames(dtf) && 'coef' %nin% colnames(dtf)) {
    dtf[, 'coef' := term]
  }

  out <- dtf %>%
    .[grepl('Intercept', coef), 'coef_type' := 'Intercept'] %>%
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
    .[grepl('focus_allele', coef),
      'coef_type' := 'Focus HLA allele'] %>%
    .[, coef := gsub('focus_allele', '', coef)] %>%
    .[grepl('sts_filtering', coef),
      'coef_type' := 'STS filtering'] %>%
    .[, coef := gsub('sts_filtering_TRUE', 'STS', coef)] %>%
    # .[, coef := gsub('sts_filtering_FALSE', 'STS', coef)] %>%
    .[grepl('patient_inclusion', coef),
      'coef_type' := 'Exclusion of IT resistant patients'] %>%
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
    .[grepl('mean_score', coef),
      'coef_type' := 'Presentation score computation'] %>%
    .[, coef := gsub('overlap_varmean_score_', '', coef)] %>%
    .[grepl(':', coef), 'coef_type' := 'interaction term'] %>%
    { . }

  out$coef[out$coef_type == 'Focus HLA allele'] %<>%
    quickMHC::ppHLA()
  out$coef[out$coef_type == 'IE focus allele'] %<>%
    quickMHC::ppHLA()
  out$coef[out$coef == 'sts_filteringTRUE'] <- 'STS'
  out$coef[out$coef == 'rooney'] <-
    'Observed/expected (Rooney et al.)'
  out$coef[out$coef == 'marty'] <-
    'Recurrent oncogenic (Marty et al.)'
  out$coef[out$coef == 'driv_ess'] <-
    'No drivers/essential passengers'
  out$coef[out$coef == 'clon'] <- 'Clonal mutations only'
  out$coef[out$coef == 'twoD_sens_analysis'] <- 'All SNVs'
  out$coef[out$coef == 'strict_LOHHLA'] <-
    'Fully assessable patients only'
  out$coef[out$coef == 'LOHHLA'] <- 'LOH HLA for assessable alleles'
  out$coef[out$coef == 'none'] <- 'None'
  out$coef[out$coef == 'FDR10'] <- 'Stringent'
  out$coef[out$coef == 'FDR1'] <- 'Lenient'
  out$coef[out$coef == 'mean_score_AB'] <- 'A, B alleles'
  out$coef[out$coef == 'mean_score'] <- 'A, B, C alleles'
  out$coef <- gsub('lumb', 'LumB', out$coef)
  out$coef <- gsub('luma', 'LumA', out$coef)
  out$coef <- gsub('min', '-', out$coef)
  out$coef <- gsub('plus', '+', out$coef)
  out$coef <- gsub('an cancer', 'an-cancer', out$coef)
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
    ggrepel::geom_text_repel(data = label_dat,
      force = 10, show.legend = F) +
    scale_x_continuous(limits = c(-max_x, max_x)) +
    scale_colour_discrete(name = 'Coefficient type')
}


forest_plot <- function(x, ...) UseMethod('forest_plot')


forest_plot.data.frame <- function(t_dat,
  x_var = 'project_extended', col_var = NULL, bg_var = NULL) {

  if ('V1' %in% colnames(t_dat)) {
    y <- 'V2'; ymin <- 'V1'; ymax <- 'V3'
  } else if ('q50' %in% colnames(t_dat)) {
    y <- '`50%`'; ymin <- '`5%`'; ymax <- '`95%`'
    y <- 'q50'; ymin <- 'q05'; ymax <- 'q95'
  } else if ('estimate_corr' %in% colnames(t_dat)) {
    y <- 'estimate_corr'
    ymin <- 'estimate_corr - 1.96 * std_error_corr'
    ymax <- 'estimate_corr + 1.96 * std_error_corr'
  } else if (all(c('estimate', 'conf.low') %in% colnames(t_dat))) {
    y <- 'estimate'
    ymin <- 'conf.low'
    ymax <- 'conf.high'
  } else if (all(c('Q2.5', 'Est.Error') %in% colnames(t_dat))) {
    y <- 'Estimate'
    ymin <- 'Q2.5'
    ymax <- 'Q97.5'
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

  if ('id' %nin% colnames(t_dat)) {
    t_dat$id <- seq_along(t_dat$project_extended)
  }

  t_dat <- arrange(t_dat, tumor_type, id)
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

  pos <- position_dodge(width = .5)

  p <- p +
    geom_point(pch = 20, position = pos) +
    geom_hline(yintercept = 0, colour = 'grey50') +
    geom_linerange(position = pos) +
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

forest_plot.setting_titration <-
  function(t_dat, x_var = 'project_extended') {
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
      name = IE_labels[['delta']],
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

summarize_by_project <- function(dtf,
  method = c('delta_SE', 'brm', 'norm_scale')) {

  method <- match.arg(method)
  setDT(dtf)
  by_vars <- intersect(c('project_extended', 'id'), colnames(dtf))

  if (method == 'brm') {
    t_dat <- dtf %>%
      pipe_row_print %>%
      .[, .(
        wa_var(.SD, 'ci_lower', error_var = 'est_error'),
        wa_var(.SD, 'estimate', error_var = 'est_error'),
        wa_var(.SD, 'ci_upper', error_var = 'est_error')),
      by = by_vars] %>%
    .[order(V2)] %>%
      { . }
  } else if (method == 'norm_scale') {
    t_dat <- dtf[, {
      reldist::wtd.quantile(
        x = yr_fractional_change,
        q = c(.05, .5, .95), na.rm = FALSE,
        weight = 1 - norm_scale
      ) %>% set_names(c('q05', 'q50', 'q95')) %>% as.list()
    }, by = by_vars]
  } else if (method == 'delta_SE') {
    t_dat <- dtf %>%
      pipe_row_print %>%
      .[, .(
        wa_var(.SD, 'delta_CI_l', error_var = 'delta_SE'),
        wa_var(.SD, 'delta_mean', error_var = 'delta_SE'),
        wa_var(.SD, 'delta_CI_h', error_var = 'delta_SE')),
      by = by_vars] %>%
    .[order(V2)] %>%
      { . }
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
data_permutation_neg_control <- result_cacher(
  f = function(args, N_repeats = 1000) {
    prep <- call_func(prep_cont_IE_analyses, args)

    neg_c <- furrr::future_map_dfr(1:(N_repeats+1), function(i) {
      l_dtf <- prep$dtf
      if (i != 1) {
        l_dtf <- l_dtf %>% mutate(y_var = sample(y_var))
      }
      out <- fit_rlm_model(l_dtf)
      out$lm <- NULL
      out
    }, .id = 'i')

    args %>%
      append(neg_c[1]) %>%
      append(list(
        'perm_p' = valid_perms[, mean(.SD[1, estimate] >= .SD[-1, estimate])],
        'N_perms' = as.integer(nrow(valid_perms))
      ))
  },
  filename = function() {
    file.path(rds_dir, 'data_permutation_neg_control',
      paste0(paste(args, collapse = '-'),
        make_flag(N_repeats), '.rds'))
  },
  min_mod_time = '2021-03-04 11:00'
)


format_coef_overview_ <- function(dtf) {
  dtf$focus_allele %<>% friendly_factor
  dtf$patient_inclusion_crit %<>% friendly_factor
  # dtf$patient_inclusion_crit <-
  #   factor(dtf$patient_inclusion_crit,
  #     levels = c('', 'strict_TR', 'TR'),
  #     labels = c('none', 'FDR1', 'FDR10'))
  dtf$patient_inclusion_crit <-
    factor(dtf$patient_inclusion_crit,
      levels = rev(c('FDR10', 'FDR1', 'FDR0.01', 'none')),
      labels = rev(c('stringent', 'moderate', 'lenient', 'none'))
    )
  dtf$LOH_HLA %<>% friendly_factor
  dtf$LOH_HLA <- relevel(dtf$LOH_HLA, 'no_LOHHLA')
  dtf$sts_filtering %<>% friendly_factor
  dtf$sts_filtering <- relevel(dtf$sts_filtering, 'FALSE')
  dtf$overlap_var %<>% friendly_factor
  return(dtf)
}


#' Simulate response var with presentation score taken from observed
#' real data
#'
#' @param beta_X Also simulate presentation score when TRUE
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


rlm_plot_pick <- plot_rlm_pick <- function(pick,
  label_pos_y = 'top',
  label_pos_x = 'left') {

  stopifnot(nrow(pick) == 1)
  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    modifyList(
      list(
        # project_extended = as.character(.$project_extended),
        # tumor_type = as.character(.$project_extended)
				# tumor_type = NULL,
				# project_extended = NULL,
        'z_normalize' = F,
        redo = F
      )
    )

  ## Prepare and compile data
  prep <- call_func(prep_cont_IE_analyses, args)

  # source(file.path('~/antigenic_space', 'maarten-analyses',
  #     'immune_editing', 'continuous_IE_detection_init.R'))
  mod <- fit_rlm_(prep$dtf, simple_output = T)
  prep$dtf$weight <- mod$weights
  obs_yr <- coef(mod) %>% { .[2] / .[1] }
  if (T && !maartenutils::eps(pick$yr_fractional_change,
      obs_yr, 1e-3)) {
    stop('Unexpected delta in rlm_plot_pick')
  }

  p1 <- plot_PS_vs_yr(
      dtf = prep$dtf,
      plot_density = F,
      lm = mod,
      colour_vars = 'weight',
      shape_var = NULL
    ) + guides(color = 'none')

  if (F && grepl('rooney', pick$analysis_name)) {
    label_pos_y <- 'bottom'
  }

  if (label_pos_y == 'top') {
    y_s <- c(.95, .85, .75, .65)
  } else if (label_pos_y == 'bottom') {
    y_s <- c(.20, .15, .10, .05)
  }

  if (label_pos_x == 'left') {
    x_s <- .05
    hjust <- 0
  } else if (label_pos_x == 'right') {
    x_s <- .95
    hjust <- 1
  }

	p1 <- p1 +
    annotate_npc(
      x = x_s, y = y_s[1],
      label = parse(text = sprintf('Delta==%.2f',
          pick$yr_fractional_change)),
      hjust = hjust, vjust = 1, gp = gpar(fontsize = 6)
    ) +
    annotate_npc(
      x = x_s, y = y_s[2],
      label = glue('nMAE',
        '= {signif(pick$scale, 3)}'),
      hjust = hjust, vjust = 1, gp = gpar(fontsize = 6)
    ) +
    annotate_npc(
      x = x_s, y = y_s[3],
      label = glue('q75 of nMAE',
        '= {signif(pick$NMADR_q75, 3)}'),
      hjust = hjust, vjust = 1,
      gp = gpar(fontsize = 6)
    ) +
    annotate_npc(
      x = x_s, y = y_s[4],
      label = parse(text = sprintf('SE(hat(Delta))==%.2f',
          signif(pick$delta_SE, 3))),
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
  # browser(expr = any(is.na(t_dat)))
  # t_dat[which(is.na(t_dat), arr.ind = T)[1]]

  if ('elapsed_time' %in% colnames(t_dat) &&
      is.character(t_dat$elapsed_time)) {
    t_dat$elapsed_time %<>%
      str_replace(' secs', '') %>%
      as.numeric
  }
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


filter_level_count <- function(dtf, vars, min_level_count = 100,
  max_quant = .1) {
  level_count <- map_dfr(auto_name(vars), function(x) {
    enframe(table(dtf[[x]])) %>% mutate(var = x)
  })

  min_level_count <-
    max(quantile(level_count$value, max_quant), min_level_count)

  exclude_levels <- level_count %>%
    dplyr::filter(value < min_level_count)

  if (nrow(exclude_levels)) {
    for (i in 1:nrow(exclude_levels)) {
      dtf <- dtf[get(exclude_levels[i, ]$'var') !=
        exclude_levels[i, ]$'name']
    }
  }
  return(dtf)
}


fit_meta_lm <- function(dtf, ds_N = NULL, method = 'rem',
  grp_var = 'tumor_type',
  min_level_count = 100,
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

  dtf <- filter_level_count(dtf, vars = FE_vars,
    min_level_count = min_level_count)
  if (nrow(dtf) == 0) return(NULL)

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


to_factor <- function(fac, lev) {
  factor(lev, levels = levels(fac))
}


run_preds <- run_lmer_preds <-
  function(pick, me_mod, version = '', N_sims = 1e4) {
  if (class(me_mod) == 'lme4') {
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
  } else if (class(lm_mod) == '') {
  }

  return(out)
}


meta_lm_plot <- function(
  tumor_type = 'melanoma',
  method = c('lm', 'brem', 'interaction_lm')) {

  method <- match.arg(method)

  if (T) {
    t_dat <- f_setting_dtf
  } else {
    t_dat <- copy(f_setting_dtf[order(norm_scale), .SD[1],
      by = analysis_grp_vars])
  }

  t_dat <- t_dat %>%
    dplyr::filter(tumor_type == .env[['tumor_type']])

  if (method == 'brem') {
    fn <- file.path(rds_dir,
      glue('meta_analysis{make_flag(tumor_type)}\\
        {make_flag(method)}.rds'))
    if (IE_requires_computation(fn)) {
      mod <- fit_meta_lm(t_dat,
        FE_vars = c('VE_threshold', 'patient_inclusion_crit',
          'expression_threshold', 'LOH_HLA', 'analysis_name',
          'sts_filtering', 'percentile_rank', 'overlap_var'),
        RE_vars = c('VE_threshold', 'patient_inclusion_crit',
          'expression_threshold', 'LOH_HLA', 'percentile_rank'),
        grp_var = 'focus_allele',
        method = 'brem')
      saveRDS(mod, fn)
    } else {
      mod <- readRDS(fn)
      browser()
    }
  } else {
    mod <- fit_meta_lm(t_dat, method = method)
    if (is.null(mod)) return(NULL)
  }

  if (class(mod) == 'brmsfit') {
    coefs <-
      fixef(mod) %>%
      as.data.frame %>%
      rownames_to_column('term')
    allele_intercepts <-
      ranef(mod)$focus_allele[, , 'Intercept'] %>%
      as.data.frame %>%
      rownames_to_column('term') %>%
      dplyr::mutate(term = paste0('focus_allele', term))
    formatted_coefs <- rbind(coefs, allele_intercepts) %>%
      format_coefs()
    p_title <- parse(text = glue('{tumor_types_inv[tumor_type]}'))
  } else {
    formatted_coefs <-
      broom.mixed::tidy(mod, conf.int = T) %>%
      format_coefs() %>%
      { . }
    p_title <- parse(text = glue('{tumor_types_inv[tumor_type]}~',
        '"- adj. "*italic(R)^2=={\\
        signif(summary(mod)$adj.r.squared, 3)}'))
  }

  p1 <- formatted_coefs %>%
    forest_plot(x_var = 'label', bg_var = 'coef_type') +
    scale_x_discrete(name = '') +
    ylab(expression('Additive contribution to  '~Delta)) +
    no_legend +
    ggtitle(p_title) +
    scale_color_brewer(palette = 'Paired')

  return(p1)
}


iterative_delta <- function(model_name = 'RF', PIC = 'FDR10',
  LOHHLA = 'strict_LOHHLA', redo = T) {
  id <- glue('{make_flag(model_name)}\\
    {make_flag(LOHHLA)}{make_flag(PIC)}')

  base_pick <- f_setting_dtf[, grp_vars, with = F] %>%
    map(~to_factor(.x, levels(.x)[1])) %>%
    as.data.table %>%
    dplyr::mutate(tumor_type = NULL) %>%
    cbind(tumor_type =
      to_factor(f_setting_dtf$tumor_type,
        levels(f_setting_dtf$tumor_type)))

  fn <- file.path(rds_dir, glue('effect_size_storytelling{id}.rds'))

  if (redo || !file.exists(fn)) {
    before <- base_pick %>%
      run_preds(me_mod = me_mod, version = 'Lenient')

    base_pick$percentile_rank %<>% to_factor('1')
    base_pick$LOH_HLA %<>% to_factor(LOHHLA)
    stageI <- base_pick %>%
      run_preds(me_mod = me_mod,
        version = 'Strict affinity filter & HLA LOH')

    base_pick$VE_threshold %<>% to_factor('1')
    stageII <- base_pick %>%
      run_preds(me_mod = me_mod, version = 'Variant expression')

    base_pick$patient_inclusion_crit %<>% to_factor(PIC)
    stageIII <- base_pick %>%
      run_preds(me_mod = me_mod, version = 'IT treatment resistance')

    # base_pick$focus_allele %<>% to_factor('A1101')
    # stageIV <- base_pick %>%
    #   run_preds(me_mod = me_mod, version = 'Opportunistic focus allele')
    # base_pick$analysis_name %<>% to_factor('rooney_param_titration')

    t_dat <- rbind(before, stageI, stageII, stageIII) %>%
      dplyr::mutate(project_extended =
        rep(tumor_types_inv[as.character(base_pick$tumor_type)], 4)) %>%
      dplyr::mutate(project_extended = factor(project_extended,
          levels = unique(project_extended)[order(stageIII$mean)])) %>%
      dplyr::mutate(version = factor(version, levels = unique(version)))
    class(t_dat) <- c('setting_titration', class(t_dat))
    saveRDS(t_dat, fn)
  } else {
    t_dat <- readRDS(fn)
  }

  print_plot(
    forest_plot(t_dat),
    fn = file.path(img_loc,
      glue('tumor_type_lmer_version_titration{make_flag(id)}.png')),
    w = 17.4, h = 20)
}


fit_forest <- function(tumor_type) {
  if (!is.null(tumor_type)) {
    t_dat <- f_setting_dtf %>%
      dplyr::filter(tumor_type == .env[['tumor_type']])
  }
  t_dat_split <- initial_split(t_dat, prop = 3/4)
  t_dat_cv <- vfold_cv(training(t_dat_split))

  grp_vars <- c('focus_allele', 'VE_threshold',
    'patient_inclusion_crit', 'expression_threshold', 'LOH_HLA',
    'analysis_name', 'sts_filtering', 'percentile_rank',
    'overlap_var')

  c_formula <-
    sprintf('yr_fractional_change ~ %s',
      paste(grp_vars, collapse = ' + ')) %>%
    as.formula()

  rf_model <-
    rand_forest() %>%
    set_args(mtry = tune()) %>%
    set_engine('ranger', importance = 'impurity') %>%
    set_mode('regression')

  t_recipe <- recipe(c_formula, data = t_dat)

  rf_workflow <- workflow() %>%
    add_recipe(t_recipe) %>%
    add_model(rf_model)

  param_final <- rf_workflow %>%
    tune_grid(
      resamples = t_dat_cv,
      grid = expand_grid(mtry = 3:length(grp_vars))
    ) %>%
    select_best(metric = 'rmse')

  rf_workflow <- rf_workflow %>%
    finalize_workflow(param_final)

  rf_fit <- rf_workflow %>%
    last_fit(t_dat_split)

  RMSE <- rf_fit %>%
    collect_metrics() %>%
    filter(.metric == 'rmse') %>%
    pull(.estimate)

  R2 <- rf_fit %>%
    collect_metrics() %>%
    filter(.metric == 'rsq') %>%
    pull(.estimate)

  if (plot_reconstruction_error) {
    test_predictions <- rf_fit %>% collect_predictions()
    test_predictions <- test_predictions %>%
      mutate(frac_error = (.pred - yr_fractional_change) /
        yr_fractional_change)
    test_predictions <- test_predictions %>%
      mutate(abs_error = abs(.pred - yr_fractional_change))
    test_predictions$n_patients <-
      t_dat[test_predictions$.row, n_patients]
    test_predictions$norm_scale <-
      t_dat[test_predictions$.row, norm_scale]
    test_predictions %>%
      dplyr::filter(abs(abs_error) > .4) %>%
      print(n = 1000, max.lines = 1000, width = 100)

    RMSE_lab <- parse(text = glue('italic(RMSE)=={signif(RMSE, 3)}'))
    R2_lab <- parse(text = glue('italic(R)^2=={signif(R2, 3)}'))
    print_plot(
      qplot(y = .pred, x = yr_fractional_change, color = n_patients,
        data = test_predictions) +
        annotate_npc(RMSE_lab, .05, .95, hjust = 0, vjust = .5) +
        annotate_npc(R2_lab, .05, .85, hjust = 0, vjust = .5) +
        geom_abline(slope = 1, intercept = 0, linetype = 2) +
        xlab(expression(Delta)) +
        ylab(expression('Random forest '~hat(Delta))) +
        scale_colour_viridis_c(name = '# patients') +
        right_legend
      , fn = file.path(img_loc,
        glue('RF_test_set_error-{tumor_type}.png'))
      , w = 8.7, h = 8.7
    )
  }

  final_model <- fit(rf_workflow, t_dat)
  ranger_obj <- pull_workflow_fit(final_model)$fit
  var_imp <-
    ranger_obj$variable.importance %>%
    enframe('term', 'imp') %>%
    dplyr::mutate(tumor_type = tumor_type) %>%
    dplyr::mutate(RMSE = RMSE) %>%
    dplyr::mutate(R2 = R2) %>%
    dplyr::mutate(norm_imp = imp / sum(imp, na.rm = T)) %>%
    dplyr::mutate(project_extended = tumor_types_inv[tumor_type]) %>%
    dplyr::mutate(project_extended = factor(project_extended,
        levels = levels(f_setting_dtf$project_extended))) %>%
    { . }

  return(var_imp)
}


id_vars <- c('overlap_var',
  'patient_inclusion_crit', 'LOH_HLA', 'analysis_name',
  'focus_allele', 'reg_method', 'z_normalize',
  'expression_threshold', 'VE_threshold', 'sts_filtering',
  'percentile_rank', 'processing_threshold', 'tumor_type',
  'project_extended')


#' Count the number of sub-analysis per group, as defined by a unique
#' combination of all of the id_vars, and filter out incomplete groups
#'
#'
subset_subanalysis_groups <- function(t_dat, id_vars,
  filter_mode = 'max') {
  setDT(t_dat)
  worthwhile_settings <- t_dat[, .N, keyby = id_vars]

  if (F) {
    setkeyv(t_dat, id_vars)
    ass <- t_dat[worthwhile_settings[4], ]
    table(t_dat$sts_filtering)
    ass[, .(analysis_idx, sts_filtering, opt_string)]
    setting_dtf[, .(analysis_idx, percentile_rank, focus_allele, sts_filtering, opt_string)]
    f_setting_dtf[, .(analysis_idx, percentile_rank, focus_allele, sts_filtering, opt_string)]
    p_setting_dtf[, .(analysis_idx, percentile_rank, focus_allele, sts_filtering, opt_string)]
  }

  # worthwhile_settings <- worthwhile_settings[N > 1]
  if (filter_mode == 'max') {
    worthwhile_settings <-
      worthwhile_settings[N == max(worthwhile_settings$N)]
  }
  invisible(worthwhile_settings[, 'group_id' := 1:.N])
  setkeyv(t_dat, id_vars)
  t_dat <- t_dat[worthwhile_settings]
  # t_dat$N <- NULL

  # t_dat[, .N, by = id_vars][, table(N)]
  # t_dat

  return(t_dat)
}


#' Filter a data.frame of sub-analyses for matching groups of
#' sub-analyses that only differ in one of independent variables, such
#' as a neo-antigen filtering setting (e.g. the gene expression
#' filtering threshold). By analyzing such data subsets, the effect of
#' this var can be directly assessed from the data.
#'
#'
filter_NN_subs <- function(
  t_dat = f_setting_dtf,
  track_var,
  response_vars = 'delta_mean',
  max_groups = 1e4,
  identify_relative_outliers_by = NULL,
  identify_outliers_by = NULL,
  twoD_sens_analysis_only = T,
  restrict_highest_level = T,
  no_STS = T,
  verbose = FALSE) {

  stopifnot(length(track_var) == 1)

  if (all(c('delta_mean', 'perm_delta_mean_mean') %in%
      colnames(t_dat))) {
    ## perm_delta_mean_c_mean := delta_mean - perm_delta_mean_mean
    t_dat[, perm_delta_mean_mean :=
      delta_mean - perm_delta_mean_c_mean ]
    # t_dat[, summary(perm_delta_mean_mean)]
  }

  if (twoD_sens_analysis_only) {
    t_dat <- t_dat[analysis_name == 'twoD_sens_analysis']
  }

  if (no_STS) {
    t_dat <- t_dat[sts_filtering == F]
  }

  # table(t_dat$RNA_expression)
  if (track_var == 'VE_threshold') {
    t_dat <- t_dat[expression_threshold == 'none']
  } else if (track_var == 'expression_threshold') {
    t_dat <- t_dat[VE_threshold == 'none']
  }

  l_id_vars <-
    c(id_vars, 'project_extended', 'RNA_expression') %>%
    setdiff(c('VE_threshold', 'expression_threshold')) %>%
    intersect(colnames(t_dat)) %>%
    setdiff(track_var) %>%
    unique()
  # setdiff(id_vars, l_id_vars)
  setdiff(colnames(t_dat), l_id_vars)
  head(colnames(t_dat))
  t_dat_f <-
    subset_subanalysis_groups(t_dat = t_dat, id_vars = l_id_vars)
  # t_dat_f[, .N, .(group_id)][, table(N)]
  table(t_dat_f$group_id)
  t_dat_f[, .(group_id, analysis_idx)]
  head(t_dat_f, 2)

  t_dat <- t_dat_f

  if (F) {
    t_dat <-
      filter_coef_overview(dtf = t_dat
        , adaptive_delta_mean_c_SE_filter = T
        , print_messages = T
      )
  }

  if (!is.null(identify_relative_outliers_by)) {
    t_dat <- t_dat[, {
      m <- mean(get(identify_relative_outliers_by))
      s <- sd(get(identify_relative_outliers_by))
      Z <- (get(identify_relative_outliers_by) - m) / s
      .SD[abs(Z) < 5]
    }, by = project_extended]
  }

  if (!is.null(identify_outliers_by)) {
    t_dat[, summary(get(identify_outliers_by))]
    ## TODO adjust threshold to other kinds of variables if needed
    stopifnot(grepl('delta_mean', identify_relative_outliers_by))
    t_dat[, 'believable' := abs(get(identify_outliers_by)) <= 1]
    # browser()
    t_dat[, .N, believable]
    t_dat[, wilcox.test(wilcox_n_patients ~ believable)]
    t_dat[, wilcox.test(perm_delta_mean_c_SE ~ believable)]
    t_dat[, boxplot(perm_delta_mean_c_mean ~ believable)]
    t_dat[, boxplot(perm_delta_mean_c_SE ~ believable)]
    t_dat[believable == T, max(perm_delta_mean_c_SE)]
    t_dat <- t_dat[abs(get(identify_outliers_by)) <= 1]
    w_dat <- t_dat[abs(get(identify_outliers_by)) > 1]
  }

  ## Downsample groups to keep it manageable
  N_groups <- length(unique(t_dat$group_id))
  if (is.null(max_groups)) max_groups <- N_groups
  sel_groups <- sort(sample(1:N_groups, min(max_groups, N_groups)))
  t_dat <- t_dat[group_id %in% sel_groups]

  t_dat <- t_dat[, {
    l_dtf <- .SD[, c(track_var, response_vars, 'analysis_idx'), with = F]
    out <- suppressWarnings(melt(l_dtf,
        measure.vars = response_vars,
        variable.name = 'var'))
  }, by = c('group_id', l_id_vars)]

  pe_levs <- levels(t_dat$project_extended)
  t_dat <- t_dat %>%
    dplyr::mutate(
      project_extended =
        factor(tumor_types_inv[as.character(tumor_type)], levels =
          pe_levs)) %>%
    # dplyr::mutate(var = map_chr(as.character(var), ~var_names[[.x]]))
    # dplyr::mutate(var = map_chr(as.character(var), ~.x))
    { . }
    # t_dat$var
  # t_dat %>%
    # dplyr::mutate(var = map_chr(var, ~.x))
  # t_dat %>%
    # dplyr::mutate(var = map_chr(var, ~var_names[[.x]]))

  allowed_tumor_types <-
    t_dat[, uniqueN(group_id), by = project_extended] %>%
    .[V1 > 100, as.character(project_extended)]
  t_dat <- t_dat[project_extended %in% allowed_tumor_types]
  if (verbose) {
    print(t_dat[, uniqueN(group_id), project_extended])
  }
  if (null_dat(t_dat)) return(NULL)

  if (restrict_highest_level) {
    ## Extract all sub-analyses with the highest observed level of the
    ## titrated variable (VE_threshold) per tumor type First identify
    ## which sub-analyses and which donor summaries we need
    subset_df <-
      t_dat[, .('level' = ll(get(track_var))), by = project_extended]
    setnames(subset_df, 'level', track_var)
    setkeyv(t_dat, colnames(subset_df))
    t_dat <- t_dat[subset_df]
    stopifnot(
      t_dat[, uniqueN(get(track_var)),
        by = project_extended][, V1] == 1)
  }
  # t_dat$project_extended

  return(t_dat)
}


#' Track response variables such as delta and delta_SE as a function
#' of some other variable of interest (the track_var), ensuring a
#' valid comparison by comparing sub-analyses which are otherwise
#' identical
#'
#'
track_sum_stats_as_function_of_var <- function(
  t_dat = f_settings_dtf,
  track_var = 'VE_threshold',
  redo = F,
  identify_relative_outliers_by = NULL,
  twoD_sens_analysis_only = F,
  restrict_highest_level = F,
  no_STS = F,
  h = 10,
  identify_outliers_by = NULL,
  response_vars = c('perm_delta_mean_pc', 'delta_mean',
    'scale', 'delta_SE'),
  sort_var = response_vars[1],
  fn = NULL) {

  if (is.null(fn)) {
    response_var_string <- paste(response_vars, collapse = '_')
    fn <- file.path(IE_img_dir,
      glue('controlled_setting_titration\\
        {make_flag(track_var)}\\
        {make_flag(response_var_string)}.pdf'))
        if (!IE_requires_computation(fn) && !redo)
          return(NULL)
  }

  t_dat <- filter_NN_subs(
    t_dat = t_dat,
    track_var = track_var,
    restrict_highest_level = restrict_highest_level,
    no_STS = no_STS,
    identify_relative_outliers_by = identify_relative_outliers_by,
    twoD_sens_analysis_only = twoD_sens_analysis_only,
    identify_outliers_by = identify_outliers_by,
    response_vars = response_vars,
    max_groups = 1e4
  )

  if (null_dat(t_dat)) return(NULL)

  meds <-
    t_dat[, list('group_id' = 1, 'value' =  median(value)),
      by = list(project_extended, var, get(track_var))]
  setnames(meds, 'get', track_var)

  if (F) {
    p_dat <- meds[, .('value' = median(diff(value))),
      by = list(project_extended, var)]
    po <- p_dat[var == sort_var, .(project_extended, value)] %>%
      .[order(value)] %>%
      pull(project_extended)
    p_dat[, project_extended := factor(project_extended, levels = po)]

    p1 <- ggplot(p_dat,
      aes(x = project_extended, y = value, group = var)) +
      geom_hline(yintercept = 0) +
      geom_point() +
      geom_line() +
      right_legend +
      ggtitle(glue('Median difference between\nconsecutive levels of\n{track_var}')) +
      facet_grid(var ~ ., scale = 'free_y') +
      rotate_x_labels(45) +
      scale_x_discrete(name = '')
    lfn <- file.path(img_loc,
        glue('controlled_setting_titration-median_value_diff_between_levels\\
          {make_flag(track_var)}\\
          {make_flag(response_var_string)}.png'))
    print_plot(p1, fn = lfn, w = 8.4, h = 12)
  }

  # anns <-
  #   map(auto_name(unique(p_dat$var)), function(v) {
  #     ComplexHeatmap::HeatmapAnnotation(
  #       df = p_dat[var == v, mean_diff_sign_between_levels],
  #       which = 'row'
  #     )
  #   })
  # purrr::reduce(anns, `+`)
  # M <- dcast(p_dat, project_extended ~ var) %>%
  #   set_rownames(.$project_extended) %>%
  #   dplyr::select(-project_extended)
  # ann <- ComplexHeatmap::HeatmapAnnotation(
  #   df = M,
  #   which = 'row'
  # )
  # Heatmap(rep(0, nrow(M))) + ann
  # draw(ann, show_annotation_legend = T)

  if (F) {
    highest_lev <- meds[order(get(track_var)),
      last(unique(get(track_var)))]
    p_dat <- meds[get(track_var) == highest_lev]
    po <- p_dat[var == last(response_vars),
      .(project_extended, value)] %>%
      .[order(value)] %>%
      pull(project_extended)
    p_dat[, project_extended := factor(project_extended, levels = po)]
    p2 <- ggplot(p_dat,
      aes(x = project_extended, y = value, group = var)) +
      geom_point() +
      geom_line() +
      right_legend +
      ggtitle(glue('Value at max level of\n{track_var}')) +
      facet_grid(var ~ ., scale = 'free') +
      rotate_x_labels(45) +
      scale_x_discrete(name = '')
    lfn <- file.path(img_loc,
        glue('controlled_setting_titration-value_at_max_level\\
          {make_flag(track_var)}\\
          {make_flag(response_var_string)}.png'))
    print_plot(p2, fn = lfn, w = 8.4, h = 12)
  }

  my_grey <- 'grey50'
  t_dat[, .N, keyby = c('project_extended', track_var)]
  p <- ggplot(t_dat,
    aes_string(x = track_var, y = 'value', group = 'group_id')) +
    # geom_point(alpha = .2, size = .5, color = my_grey) +
    ggrastr::rasterise(geom_line(alpha = .1, color = my_grey),
      dpi = 300) +
    geom_hline(yintercept = 0, color = 'grey20', size = .5) +
    geom_point(data = meds, alpha = .8, size = .7, color = 'black') +
    geom_line(data = meds, alpha = .8, size = .5, color = 'black') +
    # scale_x_discrete(name = '', expand = c(0, 0)) +
    # facet_grid(var ~ project_extended, scale = 'free_y',
    #   labeller = label_parsed()) +
    right_legend +
    rotate_x_labels(45) +
    xlab(IE_labels[[track_var]])

  if (length(response_vars) == 1) {
    p <- p + ylab(IE_labels[[response_vars[1]]])
  }

  if (length(response_vars) > 1) {
    p <- p + facet_grid(var ~ project_extended, scale = 'free',
      labeller = label_parsed)
    w <- 25 / 8 * uniqueN(t_dat$project_extended) + 1
    if (w == 0) return(NULL)
    cat(track_var, 'width:', w, '\n')
  } else {
    p <- p + facet_wrap(~project_extended,
      labeller = label_parsed,
      scale = 'free', ncol = 4)
    w <- 17.4
  }

  print_plot_eval(print(p),
    width = w, height = h, filename = fn)
}



IE_source_data_permutation <- function(
  track_var = 'VE_threshold'
  , args = NULL
  , test_i = 1:nrow(t_dat)
  , check_donor_summaries = F
  , debug = F
  , redo_regressions = F
  , filter_N_levels = T
  # , debug = T
  ) {

  t_dat <- call_func(filter_NN_subs,
    modifyList(args, list(
        max_groups = NULL,
        track_var = track_var,
        twoD_sens_analysis_only = T, restrict_highest_level = T))
  )

  if (track_var == 'VE_threshold') {
    perturbs <- list(
        'vanilla' = c()
        , 'affinity' = c('percentile_rank')
        , 'VE' = c('rna_alt_read_count')
        , 'VE_affinity' = c('rna_alt_read_count', 'percentile_rank')
      )
  } else if (track_var == 'expression_threshold') {
    perturbs <- list(
        'vanilla' = c()
        , 'affinity' = c('percentile_rank')
        , 'GE' = c('rna_expresssion')
        , 'GE_affinity' = c('rna_expression', 'percentile_rank')
      )
  } else {
    stop('Not implemented')
  }

  if (check_donor_summaries) {
    ## Create the donor summaries we'll need
    t_dat[, .(hla_allele = focus_allele, analysis_idx)] %>%
      unique() %>%
      pipe_row_print %>%
      arrange(hla_allele, analysis_idx) %>%
      expand_grid(permute_input = perturbs) %>%
      purrr::pwalk(function(hla_allele, analysis_idx, permute_input) {
        create_donor_summary(
          hla_allele = hla_allele,
          analysis_idx = analysis_idx,
          include_non_NMD = F,
          # ncores = 1,
          # ncores = 4,
          permute_input = permute_input
        )
      })
    # grp_vars <- c('focus_allele', 'VE_threshold',
    #   'patient_inclusion_crit', 'expression_threshold', 'LOH_HLA',
    #   'analysis_name', 'sts_filtering', 'percentile_rank',
    #   'overlap_var', 'tumor_type')
  }

  res <-
    furrr::future_map_dfr(test_i, function(i) {
      full_res <-
        purrr::map(perturbs, function(l) {
          pan_IE_res_to_call_args(t_dat[i]) %>%
            modifyList(
              list(
                permute_input = l
                , redo = redo_regressions
                , project_extended = t_dat[i, project_extended]
              )
            ) %>%
            { call_func(test_continuous_IE, .) } %>%
            { .[['stats']][[t_dat[i, project_extended]]][[1]] }
        })
      # out <- as.data.frame(map(full_res, 'perm_delta_mean_c_mean'))
      # out <- full_res %>% unlist(recursive = F)
      out <- rbindlist(full_res, fill = T) %>%
        dplyr::mutate(permute_input =
          names(full_res)[which(!sapply(full_res, is.null))]) %>%
        dplyr::mutate(i = i)
      return(out)
    })

  res$permute_input %<>% factor(levels = names(perturbs))

  if (filter_N_levels) {
    N_levels <- res[, .N, by = i][, unique(N)]
    allowed_i <-
      res[perm_delta_mean_c_mean <= 1 &
      message == 'OK' &
      perm_delta_mean_c_mean >= -1] %>%
      .[, .N, by = i] %>%
      .[N >= (N_levels-1), i]
    res <- res[i %in% allowed_i]
  }

  return(res)
}


plot_source_data_neg_perm <- function(
  res, track_var = 'VE_threshold') {

  ## Compare all against all others
  comps <- unique(res$permute_input) %>%
    as.factor %>%
    { expand_grid(a = ., b = .) } %>%
    dplyr::filter(as.integer(a) < as.integer(b)) %>%
    pmap(function(a, b) c(as.character(a), as.character(b)))

  sum_stats <-
    res[,
      # as.list(maartenutils::mean_CI(perm_delta_mean_c_mean))
      as.list(mean_CI(perm_delta_mean_c_mean))
    , by = permute_input] %>%
    dplyr::mutate(grp = 'a') %>%
    { . }

  y_lims <-
    { ceiling(range(res$perm_delta_mean_c_mean, na.rm = T) * 10) / 10 } %>%
    abs %>%
    pmin(1) %>%
    { . * c(-1, 1) }

  # y_lims <- sum_stats %>%
  #   select(CI_l, mean, CI_h) %>% range %>% { . * 1.1 }

  p <- res %>%
    ggplot(aes(
        x = permute_input,
        alpha = perm_delta_mean_c_mean,
        y = perm_delta_mean_c_mean,
        group = i)
      ) +
    geom_line(alpha = .01, colour = 'grey50') +
    geom_point(
      data = sum_stats,
      mapping = aes(x = permute_input, y = mean),
      colour = 'grey20', size = 1,
      inherit.aes = F
    ) +
    geom_linerange(
      data = sum_stats,
      mapping = aes(x = permute_input, y = mean,
        ymin = CI_l, ymax = CI_h),
      colour = 'grey20', size = 1,
      inherit.aes = F
    ) +
    # geom_line(
    #   data = sum_stats,
    #   mapping = aes(x = permute_input, y = mean, group = grp),
    #   colour = 'grey20',
    #   inherit.aes = F
    # ) +
    ggpubr::stat_compare_means(
      paired = T
      , method = 'wilcox',
      , size = 2.5
      , comparisons = comps
      , label.y = map_dbl(seq_along(comps) - 1, ~.5 - .x * .1)
      # , label.y.npc = map_dbl(seq_along(comps), ~1 - .x * .1)
      # , label.y = c(.9, .7, .7)
      # , label.y.npc = .8
    ) +
    coord_cartesian(ylim = y_lims) +
    ylab(IE_labels[['perm_delta_mean_c_mean']]) +
    xlab('Source data permutation')

  if ('project_extended' %in% colnames(res)) {
    p <- p + facet_grid(. ~ project_extended, scale = 'free',
      labeller = label_parsed())
  }

  return(p)
}


cov_to_cor <- function(M) {
  cor_mat <- matrix(NA, nrow = nrow(M), ncol = nrow(M))
  for (i in 1:nrow(M)) {
    for (j in i:nrow(M)) {
      cor_mat[i, j] <- cor_mat[j, i] <-
        M[i, j] / prod(sqrt(diag(M)[c(i,j)]))
    }
  }
  dimnames(cor_mat) <- dimnames(M)
  return(cor_mat)
}


#' Extract highest observed level of a factor
#'
#'
ll <- function(v, min_size = 10) {
  t <- table(v)
  t <- t[t > min_size]
  names(t[length(t)])
}


top_left_legend <- theme(
  legend.position = c(0.05, .95), legend.justification = c(0, 1))

top_right_legend <- theme(
  legend.position = c(0.95, .95), legend.justification = c(1, 1))

right_legend <- theme(
  legend.position = 'right', legend.justification = c(.5, .5))

bottom_legend <- theme(
  legend.direction = 'horizontal',
  legend.position = 'bottom', legend.justification = c(.5, .5))

no_legend <- theme(legend.position = 'none')

ncat <- function(...) {
  cat(..., '\n')
}

calcGseaStat <- function (stats, selectedStats, gseaParam = 1,
  returnAllExtremes = FALSE,
    returnLeadingEdge = FALSE, scoreType = c("std", "pos", "neg"))
{
    scoreType <- match.arg(scoreType)
    S <- selectedStats
    r <- stats
    p <- gseaParam
    S <- sort(S)
    m <- length(S)
    N <- length(r)
    if (m == N) {
        stop("GSEA statistic is not defined when all genes are selected")
    }
    rAdj <- abs(r[S])^p
    NR <- sum(rAdj)
    if (NR == 0) {
        rCumSum <- seq_along(rAdj)/length(rAdj)
    }
    else {
        rCumSum <- cumsum(rAdj)/NR
    }
    tops <- rCumSum - (S - seq_along(S))/(N - m)
    if (NR == 0) {
        bottoms <- tops - 1/m
    }
    else {
        bottoms <- tops - rAdj/NR
    }
    maxP <- max(tops)
    minP <- min(bottoms)
    switch(scoreType,
      std = geneSetStatistic <- ifelse(maxP == -minP, 0,
        ifelse(maxP > -minP, maxP, minP)),
      pos = geneSetStatistic <- maxP,
      neg = geneSetStatistic <- minP)
    if (!returnAllExtremes && !returnLeadingEdge) {
        return(geneSetStatistic)
    }
    res <- list(res = geneSetStatistic)
    if (returnAllExtremes) {
        res <- c(res, list(tops = tops, bottoms = bottoms))
    }
    if (returnLeadingEdge) {
        leadingEdge <- if (maxP > -minP) {
            S[seq_along(S) <= which.max(bottoms)]
        }
        else if (maxP < -minP) {
            rev(S[seq_along(S) >= which.min(bottoms)])
        }
        else {
            NULL
        }
        res <- c(res, list(leadingEdge = leadingEdge))
    }
    res
}


#' Test factor enrichment
#'
#' Test enrichment of a discrete variable's level among the top or
#' bottom of a rank-ordered list
test_fe <- function(
  dtf,
  vn = 'focus_allele',
  reverse = FALSE) {
  # library(fgsea)
  if (!vn %in% colnames(dtf)) return(NULL)

  if (is.factor(dtf[[vn]])) {
    levs <- levels(dtf[[vn]]) %>%
      { setNames(., .) }
    levs <- levs[names(which(table(dtf[[vn]]) > 0))]
  } else {
    levs <- unique(dtf[[vn]]) %>%
      { setNames(., .) } %>%
      sort()
  }

  stopifnot(is.logical(reverse))

  if (F) {
    q1 <- function(x) quantile(x, .1)
    q9 <- function(x) quantile(x, .9)
    dtf[, list(q1(V1), median(V1), q9(V1)), by = vn]
  }

  out <-
    purrr::map_dbl(levs, function(lev) {
    suppressWarnings(calcGseaStat(
      # stats = setNames(1:nrow(dtf), dtf[[vn]]),
      # stats = as.integer(dtf[[vn]]),
      stats = if (!reverse) nrow(dtf):1 else 1:nrow(dtf),
      # stats = nrow(dtf):1,
      selectedStats = which(dtf[[vn]] == lev),
      # scoreType = 'pos'
      scoreType = 'std'
    ))}) %>%
    tibble::enframe(vn, 'ES') %>%
    dplyr::mutate(ES = ES - mean(ES)) %>%
    dplyr::arrange(desc(ES)) %>%
    { . }

  return(out)
}
# test_fe()


test_fe_dtf <- function(dtf, id_vars, sort_var = 'delta_n', reverse = F) {
  dtf_o <-
    dtf[, median(get(sort_var)),
      by = setdiff(id_vars, 'project_extended')][order(V1)]

  if (F) {
    sort(compute_median_rank(dtf_o, 'focus_allele'))
    sort(compute_median_rank(dtf_o, 'analysis_name'))
    sort(compute_median_rank(dtf_o, 'LOH_HLA'))
    sort(compute_median_rank(dtf_o, 'overlap_var'))
    sort(compute_median_rank(dtf_o, 'patient_inclusion_crit'))
  }

  dtf_o <-
    dtf_o %>%
    # dplyr::filter(abs(V1) <= 1) %>%
    .[, median(V1), by = setdiff(colnames(dtf_o), 'V1')] %>%
    .[order(V1)] %>%
    { . }

  if (any(!colnames(ds_param_grid) %in% colnames(dtf_o))) {
    dtf_o <- maartenutils::clean_columns(fh = dtf_o,
      col_names = colnames(ds_param_grid))
    dtf_o <- cbind(dtf_o, ds_param_grid[dtf_o[, analysis_idx], ])
    dtf_o <- pretty_overview_dtf(dtf_o)
  }

  ET <-
    setdiff(id_vars, 'analysis_idx') %>%
    auto_name() %>%
    purrr::map(function(vn) {
      output <- tryCatch(
        test_fe(vn = vn, dtf = dtf_o, reverse = reverse),
        error = function(e) { NULL })
      output
    }) %>%
    purrr::discard(is.null)

  return(ET)
}


CV_editing <- function(
    dtf = f_setting_dtf,
    tt = 'Thyroid',
    id_vars = c(CV_IE_settings, 'analysis_idx'),
    # sort_var = 'delta_n',
    sort_var = 'perm_delta_pq',
    verbose = T,
    reverse = F,
    min_N = 100
  ) {

  ## Prevent f_setting_dtf from getting messed up by setkeyv()
  dtf <- copy(dtf)

  all_tts <- setdiff(levels(dtf$project_extended),
    "Pan*'-'*cancer") %>%
    maartenutils::auto_name()

  test_opt_string_consistency(dtf)
  test_opt_string_consistency(f_setting_dtf)
  dtf <- pretty_overview_dtf(dtf)
  test_opt_string_consistency(dtf)
  test_opt_string_consistency(f_setting_dtf)

  other_tts <- setdiff(all_tts, c(tt, "Pan*'-'*cancer"))

  if (TRUE) {
    ## Limit to settings that are observed both within the tumor type
    ## in question and the combination of all the other tumor types
    l_id_vars <- setdiff(id_vars, 'project_extended')
    tallies <- purrr::map(list(tt, other_tts),
      function(.x) {
        dtf_l <- dtf[project_extended %in% .x]
        ## Also require each level to be observed at least min_N times
        ## al := abandoned levels
        al <- perform_settings_tally(dtf_l, l_id_vars)[N < min_N]

        out <- dtf_l[, .N, by = l_id_vars]
        if (!maartenutils::null_dat(al)) {
          for (i in seq(1, nrow(al))) {
            # out <- out[get(al$var[i]) != al$lev[i]]
            clev <- as.character(al$lev[i])
            if (is.na(clev)) browser()
            out <- out[get(al$var[i]) != clev]
            # out[get(al$var[i]) == al$lev[i]]
          }
        }
        out$N <- NULL
        return(out)
      })
    test_opt_string_consistency(f_setting_dtf)
    combined_observed_settings <-
      dplyr::inner_join(tallies[[1]], tallies[[2]], by = l_id_vars)
    test_opt_string_consistency(f_setting_dtf)
    # table(dtf$RNA_expression)
    # table(combined_observed_settings$RNA_expression)
    test_opt_string_consistency(f_setting_dtf)
    ## This is the culprit of f_setting_dtf getting messed up
    setkeyv(dtf, l_id_vars)
    test_opt_string_consistency(f_setting_dtf)
    dtf <- dtf[combined_observed_settings]
    test_opt_string_consistency(f_setting_dtf)
  }

  test_opt_string_consistency(dtf)
  test_opt_string_consistency(f_setting_dtf)

  dtf_i <- dtf[project_extended == tt]
  dtf_o <- dtf[project_extended %in% other_tts]
  test_opt_string_consistency(dtf_o)
  test_opt_string_consistency(dtf_i)
  test_opt_string_consistency(f_setting_dtf)
  # dtf_i$sts_filtering
  # dtf_o$sts_filtering

  ET_o <- test_fe_dtf(dtf_o,
    id_vars = id_vars,
    sort_var = sort_var,
    reverse = reverse)
  test_opt_string_consistency(dtf_o)
  test_opt_string_consistency(dtf_i)
  test_opt_string_consistency(f_setting_dtf)

  if (!length(ET_o)) {
    stop('No enrichment found at all')
  }

  ## Sort by max ES such that we will prioritize the strongest ES's
  ## when filtering
  ET_o <- ET_o[rev(order(map_dbl(ET_o, ~dplyr::slice_head(.x, n = 1)$ES)))]
  # ET_o['patient_inclusion_crit']
  # table(dtf_i[['patient_inclusion_crit']])

  ET_f <-
    purrr::map(ET_o, function(.x) {
      # dplyr::filter(.x, ES >= .05) %>%
      # dplyr::filter(.x, ES >= .025) %>%
      dplyr::filter(.x, ES >= .2) %>%
      dplyr::slice_head(n = 1)
    }) %>%
    purrr::keep(~nrow(.x) > 0) %>%
    map(~unlist(.x[, 1])) %>%
    { . }
  ## Manually exclude the focus allele from being a filtering
  ## criterion
  ET_f[['focus_allele']] <- NULL

  if (verbose) {
    print(tt)
    print(ET_f)
  }

  # n_subs_uf <- dtf[!project_extended %in% tt]
  # subs_uf <- copy(dtf[project_extended %in% tt])
  # table(subs_uf$VE_threshold)
  # table(n_subs_uf$VE_threshold)

  subs <- copy(dtf_i)
  if (length(ET_f) > 0) {
    for (i in 1:length(ET_f)) {
      idx <- which(subs[[names(ET_f)[i]]] == unlist(ET_f[[i]]))
      if (length(idx) > 0) {
        subs <- subs[idx, ]
      } else {
        if (verbose) {
          cat('Filtering issue. No entries for', tt,
            names(ET_f)[i], ET_f[[i]], '\n')
        }
      }
    }
    ## Inside tumor type (i) and selected (s)
    dtf_i_s <- subs
    dtf_i_s$type = 'Selected'
    ## ... unselected (us)
    ## Lazy fix for bug with Thyroid, I'm not interested in dtf_i_us
    ## ATM
    # dtf_i_us <- dtf_i[!unique(subs[, ..id_vars])]
    dtf_i_us <- dtf_i
  } else {
    dtf_i_s <- NULL
    dtf_i_us <- dtf_i
  }
  test_opt_string_consistency(dtf_i)
  dtf_i_us$type = 'Unselected'

  # table(subs$RNA_expression)
  if (verbose) {
    if (!null_dat(dtf_i_s)) {
      print(dtf_i_s[, summary(delta_n)])
    }
    print(dtf_i_us[, summary(delta_n)])
  }

  if (F) {
    ## This may change column types
    setkeyv(subs_uf, id_vars)
    stopifnot(subs_uf[unique(subs[, ..id_vars]), .N] == nrow(subs))
  }

  ET_i <- test_fe_dtf(dtf_i,
    id_vars = id_vars, sort_var = sort_var, reverse = reverse)
  test_opt_string_consistency(dtf_i)
  if (F) {
    stopifnot(is.null(subs$type))
    subs <- subs[, !duplicated(colnames(subs)), with = F]
    subs_uf <- subs_uf[, !duplicated(colnames(subs_uf)), with = F]
    which(map_chr(subs, class) != map_chr(subs_uf, class))
  }

  out <- list(
    'dtf' =
      dplyr::bind_rows(dtf_i_s, dtf_i_us) %>%
      dplyr::mutate(
        ranking = c('top', 'bottom')[as.integer(reverse)+1]),
    'ET_o' = ET_o,
    'ET_i' = ET_i
  )

  return(out)
}


compute_median_rank <- function(dtf, vn, f = median) {
  tapply(1:nrow(dtf), dtf[[vn]], f)
}


median_based_ranking <- function(dtf = analysis_idx_o, vn = 'LOH_HLA') {
  medians <- sort(compute_median_rank(dtf, vn))
  (medians[1] - medians[2]) / nrow(dtf)
}


#' Turn an enrichment table into a grob
ET2grob <- function(ET) {
  if (is.null(ET)) return(NULL)
  ET_print <-
    ET %>%
      purrr::discard(null_dat) %>%
      purrr::keep(~any(abs(.x$ES) > .05)) %>%
      imap(~bind_cols(
          type = c(.y, base::rep('', nrow(.x)-1)), .x)) %>%
      purrr::imap_dfr(function(.x, .y) {
        colnames(.x)[2] <- 'level'
        .x
      }) %>%
      # dplyr::mutate(col = 'dark') %>%
      # dplyr::mutate(seq_along(which(type != '')) %% 2)
      { . }

  p <-
    gridExtra::tableGrob(
      ET_print, rows = NULL,
      theme = gridExtra::ttheme_default(base_size = 6)
    )

  return(p)
}


#' Enrichment table to data.frame
#'
#' Convert list of enrichment factor tables to a single table
ET2dtf <- function(ET, ES_name = 'ET_i') {
  ET %>%
    imap(~bind_cols(type = rep(.y, nrow(.x)), .x)) %>%
    purrr::imap_dfr(function(.x, .y) {
      colnames(.x)[2] <- 'level'
      .x
    }) %>%
    dplyr::rename({{ES_name}} := ES) %>%
    { . }
}


compile_arg_string <- function(v) {
  stopifnot(is.list(v))
  stopifnot(!is.null(names(v)))
  purrr::imap_chr(v, ~glue('{.y}={.x}')) %>% paste(collapse = '-')
}

plot_delta_overview <- function(
  t_dat,
  version_string =
    compile_arg_string(attr(t_dat, 'filtering_args')),
  fn = file.path(IE_img_dir,
    glue::glue('{datef}_filtered_overview_{version_string}.pdf'))) {

  p_dat <-
    t_dat %>%
    dplyr::mutate(project_extended =
      forcats::fct_rev(project_extended)) %>%
    { . }

  # max_median_delta_n <-
  #   max(p_dat[, abs(median(delta_n)), by = project_extended]$V1)
  # max_median_delta_n <-
  #   max(p_dat[, max(abs(quantile(delta_n, c(.1, .9)))),
  #     by = project_extended]$V1)
  max_median_delta_n <-
    max(p_dat[, max(abs(quantile(delta_n, c(.25, .75)))),
      by = project_extended]$V1)
  # p_dat$project_extended

  p1 <-
    ggplot(p_dat, aes(y = project_extended, x = delta_n)) +
    geom_vline(xintercept = 0) +
    # geom_violin(draw_quantiles = c(.25, .5, .75),
    #   # color = 'indianred3',
    #   # fill = 'indianred3',
    #   # show.legend = FALSE,
    #   alpha = .2) +
    geom_boxplot(outlier.shape = NA) +
    # stat_summary(
    #   fun = 'median', colour = 'indianred3',
    #   geom = 'point',
    #   size = .2,
    #   show.legend = FALSE
    # ) +
    coord_cartesian(xlim = max_median_delta_n * c(-1, 1) * 1.05) +
    # scale_x_discrete(name = expression(Delta[n])) +
    xlab(expression(Delta[n])) +
    scale_y_discrete(
      name = '',
      # labels = parse(text = levels(p_dat$project_extended))) +
      labels = function(v) parse(text = v)) +
    theme()

  print_plot_eval(print(p1),
    width = 8.7, height = 10,
    filename = fn)
}


plot_p_overview <- function(
  t_dat,
  effect_var = 'delta',
  version_string =
    compile_arg_string(attr(t_dat, 'filtering_args')),
  fn = file.path(IE_img_dir,
    glue::glue('{datef}_p_overview_{version_string}.pdf'))) {

  max_median_effect_var <-
    max(t_dat[, max(abs(quantile(get(effect_var), c(.25, .75)))),
      by = project_extended]$V1)

  p1 <- t_dat[, .N,
    keyby = .(project_extended,
      # 'sign' = sign(get(effect_var)),
      'sign' = cut(
        get(effect_var),
        breaks = sort(unique(c(0, seq(-max_median_effect_var,
          max_median_effect_var, length.out = 7))))
      ),
      'p_val_bin' = cut(p_val, c(0, 0.05, 1)))] %>%
    .[!is.na(sign)] %>%
      # 'p_val_bin' = cut(p_val, seq(0, 1, by = .05)))] %>%
    # { . } %>%
    .[p_val_bin == '(0.05,1]', N := N/19] %>%
    # as_tibble() %>%
    # dplyr::mutate(N = case_when(
    #   p_val_bin == '(0.05,1]' ~ N / 19,
    #   TRUE ~ N,
    # )) %>%
    # { . }
    # .[, 'x_var' := paste0(p_val_bin)] %>%
    # .[, 'sign' := factor(sign, levels = c(-1, 1),
    #   labels = c('editing', 'enrichment'))] %>%
    ggplot(aes(x = sign, y = N, fill = p_val_bin)) +
      geom_col() +
      facet_wrap(~project_extended, labeller = label_parsed) +
      theme(legend.position = 'right') +
      rotate_x_labels(45) +
      scale_fill_discrete(name = 'Binned p-value') +
      xlab('Sign of delta') +
      ylab('# sub-analyses')

  print_plot_eval(print(p1),
    width = 17.4, height = 10,
    filename = fn)
}


param_sensitivity <- function(dtf,
  method = 'iqr',
  fn = file.path(IE_img_dir,
    glue::glue('cov_parameter_sensitivity.pdf'))) {

  p_dat <- id_vars %>%
    # setdiff('project_extended') %>%
    setdiff('tumor_type') %>%
    purrr::map_dfr(function(c_var) {
      # c_var = 'z_normalize'
      # c_var = 'project_extended'
      # c_var = 'analysis_name'

      l_id_vars <- id_vars %>%
        setdiff(c_var) %>%
        setdiff('tumor_type') %>%
        # c('project_extended') %>%
        unique()

      t_dat <- subset_subanalysis_groups(
        t_dat = dtf[!is.na(delta_n)],
        id_vars = l_id_vars)

      if (method == 'cov') {
        t_dat <- t_dat[,
          .(
            # 'cov' = sd(delta_n, na.rm = T) / mean(delta_n, na.rm = T),
            'cov_n' = sd(delta_n, na.rm = T) /
              abs(mean(delta_n, na.rm = T)),
            'N' = .N), keyby = l_id_vars]
      } else if (method == 'iqr') {
        t_dat <- t_dat[,
          .(
            # 'cov' = sd(delta_n, na.rm = T) / mean(delta_n, na.rm = T),
            # 'cov_n' = sd(delta_n, na.rm = T) /
            #   abs(mean(delta_n, na.rm = T)),
            'cov_n' = IQR(delta_n, na.rm = T) /
              abs(median(delta_n, na.rm = T)),
            'N' = .N), keyby = l_id_vars]
      }

      ## Debug
      if (F && c_var == 'analysis_name') {
        pick <- t_dat[N == 5][1]
        setkeyv(dtf, intersect(colnames(dtf), colnames(pick)))
        setkeyv(pick, intersect(colnames(dtf), colnames(pick)))
        setkeyv(setting_dtf, intersect(colnames(dtf), colnames(pick)))
        ## They're all the same, strange..!
        dtf[pick]
        setting_dtf[pick]

        source(file.path(IE_root, 'continuous_IE_detection_init.R'))
        call_func(test_continuous_IE, as.list(pick))

      } else {
        t_dat <- maartenutils::clean_columns('', t_dat, l_id_vars)
      }


      # table(dtf$project_extended)
      # any(t_dat$N > 1)
      if (all(is.na(t_dat$cov_n))) return(NULL)
      # summary(t_dat$cov)
      t_dat$var <- c_var
      return(t_dat)
    })

  setDT(p_dat)
  var_titles_l <- setNames(names(var_titles), var_titles)
  p_dat$var_print <- forcats::fct_recode(p_dat$var, !!!var_titles_l)
  # p_dat[, mean(cov, na.rm = T), by = var]
  # p_dat[, quantile(cov, .9, na.rm = T), by = var]
  p_dat$var_print <-
    factor(p_dat$var_print,
      levels = p_dat[, quantile(cov_n, .5, na.rm = T), by =
        var_print][order(V1)]$var_print)

  # unique(p_dat$var)

  if (F) {
    p_dat[, quantile(cov_n, .99, na.rm = T), keyby = var_print][order(V1)]
    p_dat[, quantile(cov_n, .75, na.rm = T), keyby = var_print][order(V1)]

    v = 'analysis_name'
    l_p_dat <- p_dat[var %in% v]
    table(l_p_dat$N)
    summary(l_p_dat$cov_n)
    l_p_dat[N == 5]
    print_plot_eval(with(l_p_dat, plot(cov_n, N)),
      width = 17.4, height = 10,
      filename = file.path(IE_img_dir, 'test.pdf'))
  }

  eps = .001
  eps = .01
  eps = .1
  eps = .05
  eps = .025

  max_x <- p_dat[, quantile(cov_n, .75, na.rm = T),
    by = var_print][, max(V1) * 1.1]
  # unique(p_dat$var)
  p <-
    # list(c('project_extended'),
    #   setdiff(unique(p_dat$var), 'project_extended')) %>%
    IE_var_groups %>%
    purrr::imap(function(v, i) {
      # v = 'overlap_var'
      l_p_dat <- p_dat[var %in% v]
      setDT(l_p_dat)
      # summary(l_p_dat$cov_n)
      # max_x <- l_p_dat[, quantile(cov_n, .75, na.rm = T),
      #   by = var_print][, max(V1) * 1.1]
      # print(max_x)
      p <- ggplot(l_p_dat,
        aes(y = var_print, x = cov_n)) +
        # geom_violin(
        #   # adjust = .1,
        #   bw = 1e-6,
        # #   draw_quantiles = c(.25, .5, .75)) +
        #   draw_quantiles = c(.10, .5, .9)) +
        # geom_vline(xintercept = 0) +
        coord_cartesian(
          # xlim = quantile(p_dat$cov_n, c(eps, 1-eps), na.rm = T)
          # xlim = c(0, quantile(l_p_dat$cov, 1-eps, na.rm = T))
          xlim = c(0, max_x)
        ) +
        geom_boxplot(outlier.shape = NA, fill = IE_pal[1], alpha = 1) +
        xlab('') +
        ylab('') +
        theme()
        # ylab('Variable')
      if (i == last(names(IE_var_groups))) {
        p <- p + xlab('Coefficient of variation')
      }
      return(p)
    })

  library(patchwork)
  print_plot_eval(print(patchwork::wrap_plots(p, ncol = 1,
        heights = map_dbl(IE_var_groups, length) %>% { ./sum(.) }
        # heights = c(.1, .9)
    )),
    width = 8.7,
    height = 15,
    filename = fn)

}


filtering <- function(dtf, args_list = NULL) {
  dtf <- purrr::exec(filter_coef_overview,
    dtf = dtf,
    # min_patients = 20,
    print_messages = T, !!!args_list
  )
  attr(dtf, 'filtering_args') <- args_list
  return(dtf)
}


row2list <- function(pick) {
  dput(map(pick, ~ifelse(is.factor(.x), as.character(.x), .x)))
}


IE_var_groups <- list(
 'tumor_type' = 'project_extended',
 'general' = c('focus_allele', 'overlap_var'),
 'patient_inclusion' = c('patient_inclusion_crit', 'LOH_HLA'),
 'variant_selection' = c('analysis_name'),
 'neo_antigen_prediction' = c('expression_threshold',
   'VE_threshold', 'sts_filtering', 'RNA_expression',
   'percentile_rank')
)

IE_var_groups <- list(
 'tumor_type' = 'project_extended',
 'general' = c('focus_allele', 'overlap_var',
   'patient_inclusion_crit', 'LOH_HLA',
   'analysis_name'),
 'neo_antigen_prediction' = c('expression_threshold',
   'VE_threshold', 'sts_filtering', 'RNA_expression',
   'percentile_rank')
)


annotate_IE_var_groups <- function(dtf) {
  setDT(dtf)
  dtf[, 'var_group' := '']
  for (vg in names(IE_var_groups)) {
    dtf[type %in% IE_var_groups[[vg]], 'var_group' := rep(vg, .N)]
  }
  return(dtf)
}


plot_subanalysis <- function(
  pick,
  dtf_mod_f = NULL,
  title_by_acf = FALSE,
  f = pryr::partial(fit_glm_log_,
    simple_output = TRUE,
    fit_CYT = FALSE, fit_offset = TRUE, fit_PS = TRUE),
  return_alt = FALSE) {

  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    { . }
  prep1 <- call_func(prep_cont_IE_analyses,
    modifyList(args, list('permute_y' = F)))
  dtf <- prep1$dtf
  if (!is.null(dtf_mod_f)) {
    dtf <- dtf_mod_f(dtf)
  }

  dtf <- dtf[order(i)]
  fit <- f(dtf)

  if (F && !is.null(pick$delta)) {
    stopifnot(fit$delta == pick$delta)
  }

  cols <- c('#C2DEE7', '#FF7C69')

  p1 <- ggplot(fit$data,
    aes(x = i, y = c, colour = ol)) +
    ggrastr::rasterise(geom_point(alpha = .3, size = .3), dpi = 300) +
    scale_colour_gradient(low = cols[1], high = cols[2],
      name = 'Presentation\nscore', breaks = c(0, 0.1, .5, 1)) +
    ylab('Predicted\nneo-antigens\n(log10 + 1)') +
    xlab('#SNVs\n(log10 + 1)')

  lt = 'dashed'
  if (all(names(coef(fit)) == c('(Intercept)',  'i', 'i:ol'))) {
    p1 <- p1 +
      geom_abline(
        intercept = coef(fit)[['(Intercept)']],
        slope = coef(fit)[['i']],
        colour = cols[1],
        linetype = lt
      ) +
      geom_abline(
        intercept = coef(fit)[['(Intercept)']],
        slope = coef(fit)[['i']] + coef(fit)[['i:ol']],
        colour = cols[2],
        linetype = lt
      ) +
      theme(legend.position = 'right') +
      theme()
  }

  return(p1)
}


plot_subanalysis_extended <- function(
  pick,
  dtf_mod_f = NULL,
  title_by_acf = FALSE,
  f = pryr::partial(fit_glm_log_,
    simple_output = TRUE,
    fit_CYT = FALSE, fit_offset = TRUE, fit_PS = TRUE),
  f_alt = pryr::partial(fit_glm_log_,
    simple_output = TRUE,
    fit_CYT = FALSE, fit_offset = TRUE, fit_PS = FALSE),
  return_alt = FALSE) {

  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    { . }
  prep1 <- call_func(prep_cont_IE_analyses,
    modifyList(args, list('permute_y' = F)))
  dtf <- prep1$dtf
  if (!is.null(dtf_mod_f)) {
    dtf <- dtf_mod_f(dtf)
  }

  dtf <- dtf[order(i)]
  # pacman::p_load('tseries')
  fit <- f(dtf)

  if (!is.null(f_alt)) {
    fit_alt <- f_alt(dtf)
  }

  # if (F) {
  #   test_res <- call_func(test_continuous_IE,
  #     c(args)[['stats']][[as.character(pick$project_extended)]][[1]])
  #   test_res <- call_func(test_continuous_IE,
  #     c(args, list(include_perms = TRUE)))
  # }

  if (F && !is.null(pick$delta)) {
    stopifnot(fit$delta == pick$delta)
  }

  fit$data$ol_b <- cut(fit$data$ol, seq(0, 1, by = .25))

  # acf_tit <- tryCatch(
  #   round(acf(residuals(fit), pl=FALSE, lag.max = 2)$acf[, , 1][2], 3),
  #   error = function(e) { print(e); NULL })

  p1 <- ggplot(fit$data,
    aes(x = i, y = c, colour = ol)) +
    ggrastr::rasterise(geom_point(alpha = .3), dpi = 300)

  # if (title_by_acf) {
  #   p1 <- p1 + ggtitle(acf(residuals(fit), pl=FALSE)['1'])
  # }

  if (all(names(coef(fit)) == c('(Intercept)',  'i', 'i:ol'))) {
    p1 <- p1 +
      geom_abline(
        intercept = coef(fit)[['(Intercept)']],
        slope = coef(fit)[['i']],
        colour = '#C2DEE7'
      ) +
      geom_abline(
        intercept = coef(fit)[['(Intercept)']],
        slope = coef(fit)[['i']] + coef(fit)[['i:ol']],
        # colour = 'purple'
        colour = '#FF7C69'
      ) +
      theme(legend.position = 'right') +
      theme()

    if (return_alt) {
      p1a <- ggplot(fit$data,
        aes(x = i, y = c, colour = ol)) +
        ggrastr::rasterise(geom_point(alpha = .3), dpi = 300) +
        geom_quantile(quantiles = .5, colour = 'grey10') +
        # geom_smooth(method = 'lm') +
        theme(legend.position = 'right') +
        facet_wrap(~ol_b, nrow = 1)

      if (title_by_acf) {
        p1a <- p1a + ggtitle(acf_tit)
      }

      if (F) {
        print_plot_eval(print(p1a),
          width = 17.4, height = 10,
          filename = file.path(IE_img_dir, 'test_scatter.pdf'))
      }

      return(p1a)
    }

    if (F) {
      ## Sanity checks
      yhat1 <- unname(predict(fit))
      yhat2 <- coef(fit)[['(Intercept)']] +
        coef(fit)[['i']] * fit$data$i +
        coef(fit)[['i:ol']] * fit$data$i * fit$data$ol
      stopifnot(maartenutils::eps(yhat1, yhat2))
      stopifnot(maartenutils::eps(residuals(fit),
          fit$data$c - yhat2))
    }

    if (F) {
      yhat_basic <- coef(fit)[['(Intercept)']] +
        coef(fit)[['i']] * fit$data$i
      yhat_ps <- yhat_basic +
        coef(fit)[['i:ol']] * fit$data$i * fit$data$ol
      fit$data$improv <-
        with(fit$data, (yhat_ps - yhat_basic) / yhat_basic)
    } else {
      fit$data$improv <-
        (residuals(fit) - residuals(fit_alt)) / residuals(fit_alt)
      # mean(residuals(fit_alt)) - mean(residuals(fit))
      # summary(fit$data$improv)
      # deviance(fit) - deviance(fit_alt)
    }

    p_eps = .03
    # fit$data$ol_b <- cut(fit$data$ol, seq(0, 1, by = .25))
    fit$data$ol_b <- cut(fit$data$ol, seq(0, 1, by = .1))
    p2 <- ggplot(fit$data, aes(x = ol_b, y = improv)) +
      geom_boxplot() +
      coord_cartesian(
        # ylim = c(quantile(fit$data$improv, p_eps), 0)
        ylim = quantile(fit$data$improv, c(p_eps, 1-p_eps))
      ) +
      theme(legend.position = 'right') +
      # gg_tabula_rasa +
      scale_y_continuous(
        name = 'Fractional change in residuals',
        expand = c(0, 0)
      ) +
      scale_x_discrete(name = 'Binned presentation score',
        expand = c(0, 0)) +
      theme()
    if (F) {
      print_plot_eval(print(p2),
        width = 17.4, height = 10,
        filename = file.path(IE_img_dir, 'test_residuals.pdf'))
    }

    test_out <- glm_fit_model(dtf, include_perms = TRUE)
    c_dat <- attr(test_out, 'perms')

    p_eps = .01

    marginal_boxplot <-
      ggplot(c_dat,
      aes_string(x = '1', y = 'y_var')) +
      geom_boxplot(alpha = .5, width = .5, outlier.size = .5,
        fill = 'firebrick3') +
      scale_x_continuous(name = 'Anchor carriers', breaks = c(1)) +
      scale_y_continuous(name = '', limits = c(0, y_max)) +
      theme_fas(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

    if (T) {
      p_scat <- ggplot(c_dat, aes(x = intercept, y = rc)) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        ggrastr::rasterise(
          geom_point(size = 1, alpha = .5, colour = 'grey10'),
          dpi = 300) +
        coord_cartesian(
          # ylim = c(quantile(fit$data$improv, p_eps), 0)
          xlim = quantile(c_dat$intercept, c(p_eps, 1-p_eps))
        ) +
        annotate('point', x = test_out$intercept, y = test_out$rc,
          color = 'red') +
        theme()
      print_plot_eval(print(p_scat),
        width = 17.4, height = 10,
        filename = file.path(IE_img_dir, 'perm_scatter.pdf'))
    }

    p3 <- ggplot(fit$data, aes(x = ol, y = improv)) +
      # geom_abline(intercept = 0, slope = 1) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_density_2d_filled(show.legend = FALSE) +
      ggrastr::rasterise(
        geom_point(alpha = .5, colour = 'grey90'), dpi = 300) +
      scale_colour_viridis_c() +
      coord_cartesian(
        # ylim = c(quantile(fit$data$improv, p_eps), 0)
        ylim = quantile(fit$data$improv, c(p_eps, 1-p_eps))
      ) +
      # geom_smooth(span = 100) +
      theme(legend.position = 'right') +
      scale_x_continuous(name = 'Presentation score',
        expand = c(0, 0)) +
      scale_y_continuous(
        name = 'Fractional change in residuals',
        expand = c(0, 0)
      ) +
      theme()
    if (F) {
      print_plot_eval(print(p3),
        width = 17.4, height = 10,
        filename = file.path(IE_img_dir, 'test_residuals.pdf'))
    }

  } else {
    rlang::abort('Unrecognized fit object')
  }

  if (F && !is.null(pick$perm_rc_median)) {
    p1 <- p1 +
      geom_abline(
        intercept = fit$offset - pick$perm_offset_median,
        slope = (fit$intercept - pick$perm_intercept_median) +
          (fit$rc - pick$perm_rc_median),
        colour = 'red'
      )
  }

  return(list(p1, p2, p3))
}


perform_settings_tally <- function(dtf, id_vars,
  by_vars = 'project_extended') {

  if (maartenutils::null_dat(dtf)) return(NULL)
  setDT(dtf)

  ## Assume all possible levels are observed at least once
  total_N <-
    prod(map_dbl(id_vars, function(v) dtf[, uniqueN(get(v))]))

  out <-
    purrr::map(id_vars, function(lv) {
      N_levs <- dtf[, uniqueN(get(lv))]
      if (N_levs == 1) return(NULL)
      dtf[, .(
        'var' = lv,
        'N' = .N,
        'frac_observed' = .N * N_levs / total_N),
        keyby = c('lev' = lv, by_vars)] %>%
        .[, 'lev' := get(lv)] %>%
        .[, (lv) := NULL]
    }) %>%
    purrr::discard(is.null) %>%
    rbindlist() %>%
    { . }

  return(out)
}


#' As the VE_threshold and expression_threshold variables are mutually
#' exclusive (one and/or the other is 'none' is always true), combine
#' them into a single variable
#'
#'
recode_RNA_expression <- function(dtf) {
  setDT(dtf)
  vns <- c('expression_threshold', 'VE_threshold')
  prefixes <- c('TPM=', 'VE=')

  for (i in seq_along(vns)) {
    orig_levs <- levels(dtf[, get(vns[i])]) %>%
      setdiff('none')
    if (!any(grepl(prefixes[i], orig_levs))) {
      new_levs <- orig_levs %>%
        { setNames(., paste0(prefixes[i], .)) } %>%
        { c('none' = 'none', .) }
      dtf[, (vns[i]) := forcats::fct_recode(get(vns[i]), !!!new_levs)]
      # levels(dtf[[vns[i]]])
      # dtf[, table(get(vns[i]))]
    }
  }

  dtf[, 'RNA_expression' := ifelse(
    expression_threshold != 'none',
    as.character(expression_threshold),
    ifelse(VE_threshold != 'none', as.character(VE_threshold), 'none'))]
  levs <- c(levels(dtf$expression_threshold),
    levels(dtf$VE_threshold)) %>% unique()
  dtf[, 'RNA_expression' := factor(as.character(RNA_expression),
    levels = levs)]

  if (F) {
    with(dtf, table(RNA_expression, VE_threshold))
    with(dtf, table(RNA_expression, expression_threshold))
  }

  return(dtf)
}


plot_overview_volcano <- function(
  dtf,
  effect_var = 'delta',
  p_var = 'p_val',
  # p_var_trans_f = function(x) -log10(x + 1e-5),
  # p_var_trans_f = function(x) -log10(x),
  p_var_trans_f = function(x) x,
  geom_mod = function(x) ggrastr::rasterise(x, dpi = 300),
  # geom_mod = identity,
  add_sig = TRUE) {

  if (maartenutils::null_dat(dtf)) return(NULL)

  # dtf[which(dtf[[p_var]] < .05), ]


  p_eps = .01
  p_eps = .025
  p_eps = .05
  if (add_sig) {
    if (grepl('pq$', p_var)) {
      positions <- data.table(
        delta_sign = c(-1L, -1L, -1L, 1L, 1L, 1L),
        'sig' = c('neg', 'none', 'pos', 'neg', 'none', 'pos'),
        xpos = c(p_eps, p_eps, p_eps, 1-p_eps, 1-p_eps, 1-p_eps),
        ypos = c(p_eps, .5, 1-p_eps, p_eps, .5, 1-p_eps)
      )
    } else {
      positions <- data.table(
        delta_sign = c(-1L, -1L, 1L, 1L),
        'sig' = c('neg', 'none', 'neg', 'none'),
        xpos = c(p_eps, p_eps, 1-p_eps, 1-p_eps),
        ypos = c(p_eps, 1-p_eps, p_eps, 1-p_eps)
      )
    }
    quartile_tallies <-
      dtf[, .N, keyby = .(
        project_extended,
        'delta_sign' = sign(get(effect_var)),
        'sig' = ifelse(
          p_var_trans_f(get(p_var)) < p_var_trans_f(0.05),
          'neg',
          ifelse(p_var_trans_f(get(p_var)) > p_var_trans_f(0.95),
            'pos', 'none'))
        )] %>%
    .[, 'frac' := N / sum(N), by = project_extended]

    quartile_tallies$frac_f <-
      unlist(maartenutils::fancy_scientific(quartile_tallies$frac))
    setkeyv(quartile_tallies, c('delta_sign', 'sig'))
    ## Left-join annotation coordinates
    quartile_tallies <- quartile_tallies[positions]
    quartile_tallies <- quartile_tallies[!is.na(project_extended)]

    lower_thresh_dat <- quartile_tallies[delta_sign == -1 & sig == 'neg'] %>%
      .[, yv := p_var_trans_f(.05)]
    upper_thresh_dat <- quartile_tallies[delta_sign == 1 & sig == 'pos'] %>%
      .[, yv := p_var_trans_f(.95)]
  }

  dtf[['p_var_trans']] <- p_var_trans_f(dtf[[p_var]])

  p <- ggplot(dtf, aes_string(
      x = effect_var,
      y = 'p_var_trans'
    ))

  if (add_sig) {
    p <- p +
      geom_hline(data = lower_thresh_dat, mapping = aes(yintercept = yv), color = 'indianred3') +
      geom_hline(data = upper_thresh_dat, mapping = aes(yintercept = yv), color = 'indianred3')
  }

  p_eps_x <- .01
  p <- p +
    geom_vline(xintercept = 0, color = 'grey10') +
    # geom_density_2d_filled() +
    geom_mod(geom_point(size = .01, alpha = .1, color = 'grey20')) +
    # geom_density_2d_filled() +
    coord_cartesian(
      ylim = quantile(dtf[['p_var_trans']], c(0, 1)),
      xlim = quantile(dtf[[effect_var]], c(p_eps_x, 1-p_eps_x))
    ) +
    facet_wrap(~project_extended, labeller = label_parsed) +
    theme(
      legend.position = 'right',
      plot.margin = margin(.1, .1, .1, .1, 'cm')
    ) +
    ggpp::geom_text_npc(data = quartile_tallies,
      inherit.aes = FALSE, color = 'indianred3', size = 2.5,
      alpha = .95,
      parse = TRUE,
      mapping = aes(npcx = xpos, npcy = ypos, label = frac_f)) +
    theme()
  if (F) {
    print_plot_eval(print(p),
      width = 17.4, height = 8,
      filename = file.path(IE_img_dir,
        glue::glue('volcano_test.pdf')))
  }

  return(p)
}


plot_delta_perm_relationship <- function(dtf, p_eps = .025,
  title = NULL) {
  p <- ggplot(dtf, aes(x = delta, y = perm_delta_median)) +
    # geom_density_2d_filled() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    ggrastr::rasterise(
      geom_point(alpha = .05, size = .1, colour = 'grey10'), dpi = 300) +
    geom_density_2d() +
    coord_cartesian(
      xlim = quantile(dtf$delta, c(p_eps, 1-p_eps)),
      ylim = quantile(dtf$perm_delta_median, c(p_eps, 1-p_eps))
    ) +
    theme(legend.position = 'right')
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


tidy <- function(...) UseMethod('tidy')


tidy.lmp <- function (x, conf.int = TRUE, conf.level = 0.95,
  ...) {
  ret <- as_tibble(summary(x)$coefficients, rownames = 'term')
  colnames(ret) <- c('term', 'estimate', 'iter', 'p.value')
  coefs <- tibble::enframe(stats::coef(x), name = 'term',
    value = 'estimate')
  ret <- left_join(coefs, ret, by = c('term', 'estimate'))

  if (conf.int) {
    ci <- broom:::broom_confint_terms(x, level = conf.level)
    ret <- dplyr::left_join(ret, ci, by = 'term')
  }

  return(ret)
}


plot_CV_results <- function(CV_settings, yvar = 'delta_n', ranking = 'top') {
  CV_settings %>%
  # map(~.x[[ranking]]) %>%
  purrr::map(function(x) {
    tt <- as.character(x$dtf$project_extended)[1]

    nrow(x$dtf) == p_setting_dtf[project_extended == tt, .N]
    x$dtf$type <- factor(x$dtf$type, levels = c('unselected', 'selected'))

    p_eps = .025
    p1 <-
      ggplot(x$dtf, aes_string(x = 'type', y = yvar, fill = 'type')) +
      geom_hline(yintercept = 0) +
      geom_boxplot(outlier.shape = NA) +
      # geom_violin(draw_quantiles = c(.25, .5, .75),
      #   # color = 'indianred3',
      #   # fill = 'indianred3',
      #   # show.legend = FALSE,
      #   alpha = .2) +
      # stat_summary(
      #   fun = 'median', colour = 'red',
      #   geom = 'point',
      #   size = 2,
      #   show.legend = FALSE
      # ) +
      scale_fill_manual(
        values = c(
          'Remaining' = 'indianred3',
          'Filtered out' = 'grey70'),
        guide = 'none') +
      coord_cartesian(
        ylim = quantile(x$dtf$delta_n, c(p_eps, 1-p_eps))
      ) +
      ggtitle(parse(text = tt)) +
      # facet_grid(~ranking) +
      xlab('') +
      ylab(expression(Delta)) +
      theme(legend.position = 'right')
    # ggsave(p1,
    #   width = 17.4, height = 10,
    #   unit = 'cm',
    #   filename = file.path(IE_img_dir,
    #     glue::glue('ranking_comp_{ranking}.pdf')))

    # p2 <- ET2grob(x$ET_o)
    # p3 <- ET2grob(x$ET_i)

    merged <-
      dplyr::inner_join(
        ET2dtf(x$ET_o, ES_name = 'ES_o'),
        ET2dtf(x$ET_i, ES_name = 'ES_i'),
        by = c('type', 'level')
      )

    # merged %>%
    #   dplyr::mutate(comp = factor(sign(ES_o) == sign(ES_i))) %>%
    #   dplyr::group_by(comp) %>%
    #   dplyr::summarize(
    #     min = min(abs(ES_o), abs(ES_i)),
    #     max = max(abs(ES_o), abs(ES_i))
    #   )

    if (FALSE) {
      threshold <-
        merged %>%
          dplyr::mutate(comp = factor(sign(ES_o) == sign(ES_i))) %>%
          dplyr::filter(comp == FALSE) %>%
          dplyr::summarize(max(pmax(abs(ES_o), abs(ES_i)))) %>%
          unlist() %>%
          unname() %>%
          { . }
    } else {
      threshold <- .25
    }

    hl_dtf <-
      merged %>%
      dplyr::filter(abs(ES_o) >= threshold | abs(ES_i) >= threshold)

    l_type_pal <- IE_settings_pal
    lame_types <-
      merged %>%
      dplyr::group_by(type) %>%
      dplyr::summarize(max_ES = max(pmax(abs(ES_o), abs(ES_i)))) %>%
      # dplyr::summarize(max_ES = max(abs(ES_o))) %>%
      # dplyr::summarize(max_ES = max(ES_o)) %>%
      dplyr::filter(max_ES < .2) %>%
      pull(type)
    l_type_pal[lame_types] <- 'grey70'

    xrange <- max(abs(range(merged$ES_i)), .25) %>%
      { c(-1 *., .) }
    yrange <- max(abs(range(merged$ES_o)), .25) %>%
      { c(-1 *., .) }
    p4 <-
      ggplot(merged, aes(x = ES_i, y = ES_o, color = type,
          label = level)) +
        geom_rect(xmin = -threshold, xmax = threshold,
          ymin = -threshold, ymax = threshold,
          color = NA,
          alpha = .01, fill = 'indianred3') +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        geom_point() +
        coord_cartesian(xlim = xrange, ylim = yrange) +
        scale_colour_manual(name = 'Subanalysis setting',
          values = l_type_pal) +
        # xlab(parse(text = glue::glue('Enrichment factor\nwithin {tt}')) +
        xlab(glue::glue('Enrichment factor\nwithin tumor type')) +
        ylab('Enrichment factor\noutside of tumor type') +
        ggrepel::geom_text_repel(data = hl_dtf, show.legend = F,
          size = 2) +
        # geom_vline(xintercept = threshold, linetype = 2) +
        # geom_hline(yintercept = threshold, linetype = 2) +
        # geom_vline(xintercept = -threshold, linetype = 2) +
        # geom_hline(yintercept = -threshold, linetype = 2) +
        theme(legend.position = 'left')
    if (F) {
      ggsave(p4,
        width = 17.4, height = 10,
        unit = 'cm',
        filename = file.path(IE_img_dir,
          glue::glue('ranking_comp_{ranking}.pdf')))
    }

    library(patchwork)
    # design <- 'AA
    # BC'
    # return(wrap_plots(p1, p2, p3, design = design))
    return(patchwork::wrap_plots(p4, p1))
  })
}


plot_CV_overview <- function(CV_settings) {
  ## tumor type
  ## --> XXX top, bottom
  ##   --> ET_i, ET_o
  ##     --> var


  # CV_settings[[3]][[2]]$sts_filtering

  var_titles_l <- setNames(names(var_titles), var_titles)
  p_dat <-
    purrr::map(CV_settings, ~.x[['ET_i']]) %>%
    purrr::discard(is.null) %>%
    purrr::map(~ET2dtf(.x, ES_name = 'ES_i')) %>%
    imap_dfr(~cbind(.x, 'project_extended' = .y)) %>%
    setDT() %>%
    .[, var_print := forcats::fct_recode(type, !!!var_titles_l)] %>%
    annotate_IE_var_groups() %>%
    { . }

  if (any(is.na(p_dat$level))) {
    browser()
  }

  # CV_settings[[1]]$ET_i$percentile_rank[, 1]
  # CV_settings[[1]]$ET_i$sts_filtering
  # p_setting_dtf$percentile_rank
  # levels(p_setting_dtf$sts_filtering)

  level_ordering <-
    names(var_titles) %>%
    # c('project_extended') %>%
    intersect(colnames(p_setting_dtf)) %>%
    map(~p_setting_dtf[, levels(get(.x))]) %>%
    unlist(recursive = F) %>%
    unique()
  # levels(p_dat$level)
  p_dat$level <- factor(p_dat$level, levels = level_ordering)
  if (F) {
    p_dat$n_level <- factor(p_dat$level, levels = level_ordering)
    p_dat[is.na(n_level)]
  }

  if (any(is.na(p_dat$level))) {
    browser()
  }
  p_dat$project_extended <- factor(p_dat$project_extended,
    levels = levels(p_setting_dtf$project_extended))
  # p_dat[var_print == 'Affinity rank percentile threshold']
  p_dat[var_print == 'Similarity to self filter']
  unique(p_dat$var_print)

  stopifnot(nrow(p_dat[var_group == '']) == 0)

  max_plots <- p_dat[, uniqueN(var_print), by = var_group][, max(V1)]

  N_var_groups <- length(unique(p_dat$var_group))

  pacman::p_load('colorspace')
  #' Make a row of plots for each group
  #' Each plot shows one variable (var_print), with multiple levels of
  #' that variable (horizontal axis) and tumor types on the vertical
  #' axis. Aim to adjust widths of plots to the number of levels in
  #' them, such that all levels get assigned identical widths.
  #'
  plot_lists <- purrr::imap(unique(p_dat$var_group), function(vg, j) {
    l_p_dat <- p_dat[var_group == vg]
    N_plots <- l_p_dat[, uniqueN(var_print)]
    N_levs <- map_int(unique(l_p_dat$var_print),
      ~length(unique(l_p_dat[var_print == .x]$level)))
    plots <- purrr::imap(unique(l_p_dat$var_print), function(varp, i) {
      p <- ggplot(l_p_dat[var_print == varp],
        aes(x = level, y = project_extended,
          color = ES_i, size = abs(ES_i))) +
        geom_point() +
        # facet_wrap(~var_print, scale = 'free_x', nrow = 3)  +
        # scale_colour_diverging() +
        scale_colour_continuous_diverging(
          name = 'Enrichment\nstatistic',
          limits = quantile(p_dat$ES_i, c(0, 1))
        ) +
        scale_size_continuous(
          guide = 'none',
          range = c(.1, 2),
          limits = quantile(p_dat$ES_i, c(0, 1))
        ) +
        scale_y_discrete(labels = function(x) parse(text = x)) +
        rotate_x_labels(45) +
        xlab('') +
        ggtitle(varp) +
        ylab('') +
        theme(legend.position = 'left') +
        theme(plot.margin = margin(.1, .1, .1, .1, 'cm'))
      if (i > 1) {
        p <- p + theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = 'none'
        )
      } else {
        # p <- p + ggtitle(varp)
      }
      if (j < N_var_groups) {
        p <- p + theme(legend.position = 'none')
      }
      # print_plot_eval(print(p),
      #   width = 17.4, height = 15,
      #   filename = file.path(IE_img_dir,
      #     glue::glue('settings_IE_GSEA_overview.pdf')))
      return(p)
    })
    if (N_plots < max_plots) {
      white_space <- map(seq(1, max_plots - N_plots), ~
        patchwork::plot_spacer())
    } else {
      white_space <- NULL
    }
    N_levs_norm <- N_levs %>% { . / sum(.) }
    patchwork::wrap_plots(c(plots, white_space),
      guides = 'collect',
      widths = c(1/max_plots * N_levs_norm,
        rep(1/max_plots, max_plots - N_plots))) +
      patchwork::plot_annotation(title = vg, caption = 'test')
  })

  return(plot_lists)
}


#' Find groups of sub-analyses that differ by only one setting and go
#' down in delta, or any other response var, with every level increase
#' in the setting
#'
#'
find_neg_groups <- function(
  track_var = 'VE_threshold') {

  redo = F
  identify_relative_outliers_by = NULL
  twoD_sens_analysis_only = FALSE
  restrict_highest_level = FALSE
  no_STS = FALSE
  identify_outliers_by = NULL
  response_vars = c('delta')
  # sort_var = response_vars[1]
  sort_var = 'delta'

  f_setting_dtf$RNA_expression <- NULL
  f_setting_dtf <- drop_high_VE(f_setting_dtf)
  # f_setting_dtf[, RNA_expression := droplevels(RNA_expression)]
  # f_setting_dtf$analysis_idx <- NULL

  all_levs <- f_setting_dtf[, levels(get(track_var))]
  offset_lev <- all_levs[1]
  top_lev <- all_levs[length(all_levs)]
  N_levs <- f_setting_dtf[, length(levels(get(track_var)))]

  set.seed(3)
  t_dat <- filter_NN_subs(
    t_dat = f_setting_dtf,
    track_var = track_var,
    restrict_highest_level = restrict_highest_level,
    no_STS = no_STS,
    identify_relative_outliers_by = identify_relative_outliers_by,
    twoD_sens_analysis_only = twoD_sens_analysis_only,
    identify_outliers_by = identify_outliers_by,
    response_vars = response_vars,
    max_groups = 1e4
  )
  # t_dat$var <- NULL
  # t_dat[, .N, by = group_id][, table(N)]
  # t_dat[group_id == 17]
  # f_setting_dtf$intercept > 0

  ## Identify groups that go down with increase of the track var
  all_neg_groups <-
    t_dat[var == 'delta',
  .(.N, 
    'all_neg' = .SD[order(get(track_var)), all(sign(diff(value)) == -1)]),
  by = group_id] %>%
    .[all_neg == TRUE]
  stopifnot(length(unique(all_neg_groups$N)) == 1)

  N_picks = 1000
  # message(glue::glue('{response_vars[1]} will be used to'))
  picks <-
    t_dat[group_id %in% all_neg_groups[, group_id]] %>%
    # .[var == response_vars[1]] %>%
    .[get(track_var) == offset_lev] %>%
    .[var == 'delta'] %>%
    .[sample(1:nrow(.), min(N_picks, nrow(.)))]
  # t_dat[group_id == 36][order(get(track_var))]

  l_id_cols <- setdiff(colnames(picks), 
    c('var', 'value', track_var, 'analysis_idx', 'RNA_expression')) %>%
    intersect(colnames(t_dat))

  setkeyv(t_dat, l_id_cols)
  setkeyv(picks, l_id_cols)
  picks <- 
    t_dat[group_id %in% all_neg_groups[, group_id]] %>%
    .[get(track_var) == top_lev] %>%
    .[picks] %>%
    # .[VE_threshold != i.VE_threshold]
    { . }
  picks$var <- picks$i.var <- NULL
  # picks$analysis_idx <- picks$i.analysis_idx <- NULL

  colnames(picks) <- stringr::str_replace(colnames(picks), 'i\\.', 'cons_')

  return(picks)
}


diagnose_regression <- function(
  pick = NULL, 
  l_dtf = NULL, 
  do_fit = TRUE,
  fit_func = fit_glm_log_,
  prefix = NULL) {

  if (!is.null(pick)) {
    prep <- call_func(prep_cont_IE_analyses,
      pan_IE_res_to_call_args(pick))
    l_dtf <- prep$dtf
  } else {
    stopifnot(!is.null(l_dtf))
  }

  if (do_fit) {
    out <- fit_func(l_dtf, simple_output = FALSE)
  } else {
    out <- list()
  }

  ## Determine the yr 
  out$mean_yr <- l_dtf[, mean(c / i)]

  if (TRUE) {
    fit <- fit_func(l_dtf, simple_output = TRUE)
    leverage <- hatvalues(fit)
    cooksd <- cooks.distance(fit)
    out$leverage_evenness <- compute_evenness(abs(leverage))
    out$cooksd_evenness <- compute_evenness(abs(cooksd))
    # out$baseline_yr <- coef(fit)[['(Intercept)']] + 
    #       median(fit$data$i) * coef(fit)[['i']]
    median_TMB <- median(fit$data$i)
    out$baseline_NAYR <- coef(fit)[['(Intercept)']] / median_TMB + 
      coef(fit)[['i']]
    out$slope_NAYR <- 1/median_TMB * coef(fit)[['i:ol']]
  }

  if (TRUE) {
    ps_fit <- lm(i ~ ol, data = l_dtf, weights = leverage)
    # if (coef(ps_fit)['(Intercept)'] < 0) browser()
    out$ol_vs_i_rc <- coef(ps_fit)['ol']
    out$ol_vs_i_rc_n <- coef(ps_fit)['ol'] / coef(ps_fit)['(Intercept)']
  }

  if (!is.null(prefix)) {
    out <- setNames(out, paste0(prefix, names(out)))
  }
  # print(out)

  return(out)
}


set_lower_yr <- function(pick,
  fit_func = fit_glm_log_,
  factors = sort(unique(c(10^seq(-4, 4, by = 1), 1)))) {

  ## The 'cons_' variables refer to the 'conservative' setting and non
  ## prefixed setting refer to the 'aggressive' setting in the
  ## comparison

  ## Determine name of the tracked variable
  tv_cons <- 
    stringr::str_subset(colnames(pick), 'cons_') %>%
    setdiff(c('cons_value', 'cons_analysis_idx'))
  tv <- stringr::str_replace(tv_cons, 'cons_', '')

  out_aggressive <- diagnose_regression(
    pick = pick, 
    fit_func = fit_func,
    prefix = 'aggressive_')

  conservative_pick <- copy(pick)
  conservative_pick[, (tv) := get(tv_cons)]
  conservative_pick[, analysis_idx := cons_analysis_idx]
  out_conservative <- diagnose_regression(
    pick = conservative_pick, 
    fit_func = fit_func,
    prefix = 'conservative_')

  ## Do we already see a negative delta at a constant neo-antigen
  ## yield rate with the current distribution of the data?
  conservative_prep <- call_func(prep_cont_IE_analyses,
    pan_IE_res_to_call_args(conservative_pick))
  l_dtf <- copy(conservative_prep$dtf)
  l_dtf$c <- round(l_dtf$i * out_aggressive$aggressive_mean_yr, 1)
  out_constant <- diagnose_regression(
    l_dtf = l_dtf, 
    fit_func = fit_func,
    prefix = 'constant_')

  ## Do we already see a negative delta with the conservative data
  ## scaled to have a similar yield rate as the aggressive data?
  fc <- out_aggressive$aggressive_mean_yr / out_conservative$conservative_mean_yr
  l_dtf <- copy(conservative_prep$dtf)
  l_dtf$c <- round(l_dtf$c * fc, 1)
  out_scaled <- diagnose_regression(
    l_dtf = l_dtf, 
    fit_func = fit_func,
    prefix = 'scaled_')

  out <- c(out_aggressive, out_conservative, out_constant, out_scaled)

  return(out)
}


plot_subanalysis <- function(
  pick,
  dtf_mod_f = NULL,
  title_by_acf = FALSE,
  fit_func = pryr::partial(fit_glm_log_,
    simple_output = TRUE,
    fit_CYT = FALSE, fit_offset = TRUE, fit_PS = TRUE),
  return_alt = FALSE) {

  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    { . }
  prep1 <- call_func(prep_cont_IE_analyses,
    modifyList(args, list('permute_y' = F)))
  dtf <- prep1$dtf
  if (!is.null(dtf_mod_f)) {
    dtf <- dtf_mod_f(dtf)
  }

  dtf <- dtf[order(i)]
  fit <- fit_func(dtf)

  if (F && !is.null(pick$delta)) {
    stopifnot(fit$delta == pick$delta)
  }

  cols <- c('#C2DEE7', '#FF7C69')

  p1 <- ggplot(fit$data,
    aes(x = i, y = c, colour = ol)) +
    ggrastr::rasterise(geom_point(alpha = .3, size = .3), dpi = 300) +
    scale_colour_gradient(low = cols[1], high = cols[2],
      name = 'Presentation\nscore', breaks = c(0, 0.1, .5, 1)) +
    ylab('Predicted\nneo-antigens\n(log10 + 1)') +
    xlab('#SNVs\n(log10 + 1)')

  lt = 'dashed'
  if (all(names(coef(fit)) == c('(Intercept)',  'i', 'i:ol'))) {
    p1 <- p1 +
      geom_abline(
        intercept = coef(fit)[['(Intercept)']],
        slope = coef(fit)[['i']],
        colour = cols[1],
        linetype = lt
      ) +
      geom_abline(
        intercept = coef(fit)[['(Intercept)']],
        slope = coef(fit)[['i']] + coef(fit)[['i:ol']],
        colour = cols[2],
        linetype = lt
      ) +
      theme(legend.position = 'right') +
      theme()
  }

  return(p1)
}


plot_subanalysis_yr_vs_ps <- function(
  pick,
  dtf_mod_f = NULL,
  title_by_acf = FALSE,
  add_margin = FALSE,
  fit_func = pryr::partial(fit_glm_log_,
    simple_output = TRUE,
    fit_CYT = FALSE, fit_offset = TRUE, fit_PS = TRUE),
  return_alt = FALSE) {

  args <- pick %>%
    pan_IE_res_to_call_args() %>%
    { . }
  prep1 <- call_func(prep_cont_IE_analyses,
    modifyList(args, list('permute_y' = F)))
  dtf <- prep1$dtf
  if (!is.null(dtf_mod_f)) {
    dtf <- dtf_mod_f(dtf)
  }

  dtf <- dtf[order(i)]
  fit <- fit_func(dtf)
  fit$data$yr <- fit$data$c / fit$data$i

  if (F && !is.null(pick$delta)) {
    stopifnot(fit$delta == pick$delta)
  }


  p1 <- ggplot(fit$data, aes(x = ol, y = yr)) +
    # ggrastr::rasterise(geom_point(alpha = .1, size = .3), dpi = 300) +
    geom_point(alpha = .1, size = .3, color = 'grey10') +
    ylab('Neo-antigen yield rate') +
    xlab('Presentation score') +
    theme()

  if (T) {
    if (all(names(coef(fit)) == c('(Intercept)',  'i', 'i:ol'))) {
      median_TMB <- mean(fit$data$i)
      median_TMB <- median(fit$data$i)
      delta <- with(as.list(coef(fit)), 
        `i:ol` / (`(Intercept)`/median_TMB + `i`)
      )
      baseline_NAYR <- with(as.list(coef(fit)), 
        (`(Intercept)`/median_TMB + `i`)
      )
      slope_NAYR <- with(as.list(coef(fit)), 
        1/median_TMB * `i:ol`
      )
      if (F) {
        print(fit$data[, summary(c/i)])
        print(mean(fit$data$c) / mean(fit$data$i))
        print(baseline_NAYR)
        print(summary(residuals(fit)))
      }
      # baseline_NAYR <- 
      #   predict(fit, tibble(ol = 0, i = median_TMB)) / median_TMB
      # edited_NAYR <- 
      #   predict(fit, tibble(ol = 1, i = median_TMB)) / median_TMB
      p1 <- p1 +
        geom_abline(
          intercept = baseline_NAYR,
          slope =  slope_NAYR,
          color = 'indianred3'
        ) +
        theme(legend.position = 'right')
    }

    if (!is.null(pick)) {
      p_tit <- glue::glue('Baseline NAYR: {round(baseline_NAYR, 4)} \
        delta NAYR: {round(delta, 4)}')
      p1 <- p1 + ggtitle(p_tit)
    }
  }

  if (F) {
    print_plot_eval(print(p1),
      width = 6, height = 6,
      filename = file.path(IE_img_dir, 'test.pdf'))
  }

  if (add_margin) {
    pacman::p_load('ggExtra')
    maartenutils::print_plot_eval(
      call = 
        { p1 <<- ggMarginal(p1, groupColour = FALSE, type = 'histogram',
          colour = NA, groupFill = FALSE, size = 9) },
      filename = paste0(tempfile(), '.pdf')
    )
  }

  return(p1)
}


plot_ps_vs_tmb <- function(pick, ann_p_val = FALSE) { 
  prep <- call_func(prep_cont_IE_analyses,
    pan_IE_res_to_call_args(pick))
  prep$dtf$i <- log10(prep$dtf$i + 1)

  fit <- broom::tidy(lm(data = prep$dtf, i ~ ol))
  fd <- fit$estimate[2] / fit$estimate[1]

  p <- ggplot(prep$dtf, aes(y = i, x = ol)) +
    geom_point(alpha = .1, size = .3, color = 'grey10') +
    geom_smooth(method = 'lm') +
    xlab('Presentation score') +
    ylab('TMB\n(log10(x+1) transformed)') +
    theme() +
    ggtitle(glue::glue('delta TMB: {round(fd, 4)}'))

  if (ann_p_val) {
    p_val <- dplyr::filter(fit, term == 'ol') %>%
      pull(p.value)
    p <- p + annotate_npc(label = parse(text = fancy_p(p_val)), 
      size = 2,
      x = .95, y = .95, hjust = 1, vjust = 1)
  }

  return(p)
}


capitalize <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}


pick2print <- function(pick) {
  pick <- recode_RNA_expression(pick)
  l_id_vars <- setdiff(c(id_vars, 'RNA_expression'), c('project_extended', 'reg_method',
      'z_normalize', 'processing_threshold', 'VE_threshold',
      'expression_threshold'))
  l_settings <- pick[, l_id_vars, with = F]
  colnames(l_settings) <- var_titles[colnames(l_settings)]
  l_settings <- unlist(map(l_settings, as.character))
  l_settings <- 
    l_settings[setdiff(match(var_titles, names(l_settings), nomatch = NULL), NA)]
  out <- tibble::enframe(l_settings, name = 'coef_type', value = 'coef')
  out <- format_coefs(out)
  out$label <- NULL
  colnames(out) <- c('setting', 'value')

  out$value <- gsub('_', ' ', out$value)

  if (F) {
    idx <- which(out$setting == 'Tumor type')
    # out$value[idx] <- stringr::str_to_title(out$value[idx])
    out$value[idx] <- capitalize(out$value[idx])
    idx <- which(out$setting == 'IT-resistant patients')
    # out$value[idx] <- stringr::str_to_title(out$value[idx])
    out$value[idx] <- capitalize(out$value[idx])
  }

  out$value <- purrr::map_chr(out$value, capitalize)
  out$value[out$value == 'TRUE'] <- 'Yes'
  out$value[out$value == 'FALSE'] <- 'No'

  return(out)
}
# pick2print(f_setting_dtf[1, ])


compute_delta_w_offset <- function(...) UseMethod('compute_delta_w_offset')
 

compute_delta_w_offset.lm <- function(fit) {
  compute_delta_w_offset(
    l = as.list(coef(fit)), 
    median_TMB = median(fit$data$i)
  ) 
}

compute_delta_w_offset.list <- function(l, median_TMB = 1) {
  with(l, 
    `i:ol` / (`(Intercept)`/median_TMB + `i`)
  )
}
