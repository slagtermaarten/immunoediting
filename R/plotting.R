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


