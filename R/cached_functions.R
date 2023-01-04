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

  source_m_times <- gen_file_overview(rds_dir, pat = 'param_titration')$mtime
  if (!IE_requires_computation(all_coef_fn) && !redo &&
      !any(source_m_times > file.mtime(all_coef_fn)) &&
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

  # if (!include_non_ds_tallies) {
  #   all_settings %<>% 
  #     dplyr::filter(analysis_name == 'twoD_sens_analysis')
  # }

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

  ## Make sure twoD_sens_anlalysis donor summaries are available for
  ## all hla_alleles and analysis_idxs
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
      patient_CYT = r[['patient_CYT']],
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
  patient_CYT = 'all',
  LOH_HLA = 'no_LOHHLA',
  analysis_name = 'twoD_sens_analysis',
  focus_allele = 'A0201',
  analysis_idx = 1,
  reg_method = 'glm_log',
  project_extended = NULL,
  hla_sim_range = NULL,
  z_normalize = FALSE,
  permute_id = NULL,
  return_grob = FALSE,
  ncores = 1,
  check_res = F,
  partition_vars = c(),
  stats_idx = 1,
  skip_missing = FALSE,
  verbose = F,
  permute_input = c(),
  include_call = FALSE,
  ds_frac = NULL,
  return_res = TRUE,
  redo = F,
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
    patient_CYT = patient_CYT,
    overlap_var = overlap_var
  )

  ds_fn <- gen_cont_IE_ds_fn(
    hla_allele = focus_allele,
    analysis_name = analysis_name,
    analysis_idx = analysis_idx,
    permute_input = permute_input,
    overlap_var = overlap_var
  )

  if (!file.exists(ds_fn)) {
    mywarning(msg = glue::glue('{ds_fn} does not exist'))
    return(NULL)
  }

  return_cached_results <- 
    (is.null(permute_id) && is.null(ds_frac) &&
    !IE_requires_computation(o_fn, verbose = F) &&
    file.mtime(o_fn) > file.mtime(ds_fn) &&
    !redo && !check_res)

  if (is.na(skip_missing) || is.na(return_cached_results)) 
    browser()
     
  if (skip_missing || return_cached_results) {
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
    patient_CYT = patient_CYT,
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

  if (F) {
    browser()
    file.mtime(obj_fn)
  }

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
  dtf <- determine_continuous_IE_yval(
    dtf,
    replace_zeroes = replace_zeroes
  )
  # print(summary(dtf$i))

  ## Extract the  tumor projects present in this object
  dtf <- subset_project(dtf, project_extended)

  if (patient_CYT != 'all') {
    dtf <- merge_CYT(dtf)
    CYT_threshold <- quantile(dtf$CYT, .75)
    if (patient_CYT == '<=75') {
      dtf <- dtf[CYT <= CYT_threshold]
    } else if (patient_CYT == '>75') {
      dtf <- dtf[CYT > CYT_threshold]
    }
  }

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
  patient_CYT = 'all',
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

  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  }

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
      patient_CYT = patient_CYT,
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
  }, .parallel = (T && ncores > 1)) %>% rbindlist(fill = T)

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


