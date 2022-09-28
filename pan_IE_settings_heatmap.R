do_test_grobs <- T
## {{{ Init
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
# pacman::p_load(ComplexHeatmap)
if (!'fasanalysis' %in% loadedNamespaces())
  devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
library(ComplexHeatmap)
# devtools::load_all('~/libs/ComplexHeatmap')
library(dplyr)
pacman::p_load(circlize)
# devtools::install('~/libs/maartenutils')
pacman::p_load(maartenutils)
pacman::p_load(data.table)
pacman::p_load(RColorBrewer)
# ma_dir <- file.path('~/antigenic_space', 'maarten-analyses')
# rds_dir <- file.path(ma_dir, 'rds')
## The relative cell size assigned to the analyses with the lowest power
# source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_helpers.R'))
# options(error = recover)
options(error = traceback)

## }}} Init

## {{{ Helper functions
trans_size_size <- function(x, min_cell_size = .1) {
  # min_cell_size + x * (1 - min_cell_size)
  min_x <- min(x, na.rm = T)
  max_x <- max(x, na.rm = T)
  min_cell_size + ((x - min_x) / (max_x - min_x)) * (1 - min_cell_size)
}


test_grob <- function(grob, doit = F) {
  if (!doit) return(invisible())
  reset_plotting_device()
  dev.new()
  pushViewport(viewport(width = 0.9, height = 0.9))
  grid.rect()
  if (any(class(grob) %in% c('grob', 'gtable', 'gTree', 'gList'))) {
    grid.draw(grob)
  } else {
    draw(grob)
  }
  popViewport()
}


legend_header <- function(legend_title = 'Test label') {
  Legend(title = legend_title, at = c(0), labels = '',
         gap = unit(0, 'mm'),
         title_gp = gpar(fontsize = base_size, fontface = 'bold'),
         grid_height = unit(0, 'mm'))
}


#' Edit a grid object's viewport such that it will be justified left, bottom
#' with respect to its parent viewport
#'
#'
just_left <- function(grob_object) {
  stopifnot(any(class(grob_object) == 'grob'))
  grob_object$vp <- viewport(x = unit(0, 'npc'), y = unit(0, 'npc'),
                             just = c(0, 0))
  grob_object
}


extract_legend_grob <- function(grob_object) {
  if (any(class(grob_object) == 'Legends')) {
    return(grob_object@grob)
  } else {
    return(grob_object)
  }
}


#' Hclust cannot handle matrices in which for some pairs of rows and columns,
#' only 1 or fewer shared values are non-NA. This function recurrently
#' identifies the most aggravating column/row, excludes that column/row and
#' checks whether more columns/rows need to be excluded. As a side-effect, it
#' also identifies columns/rows that are all-NA.
#'
#' @param mat Matrix to investigate
#' @param min_shared_fields Minimum number of positions that are not NA in both
#' vectors in order not to flag the vector pair as problematic
#'
identify_problematic_combs <- function(mat, min_shared_fields = 1) {
  exclude_rows <- NULL
  exclude_cols <- NULL
  stopifnot(is.matrix(mat))

  ## Loop over candidate removals
  for (k in 1:nrow(mat)) {
    candidate_rows <- setdiff(1:nrow(mat), exclude_rows)
    problem_row_combs <- NULL
    for (i in candidate_rows) {
      i_idx <- which(candidate_rows == i)
      for (j in candidate_rows[i_idx:length(candidate_rows)]) {
        if (sum(!is.na(mat[i, ]) & !is.na(mat[j, ])) <= min_shared_fields) {
          problem_row_combs <- rbind(problem_row_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_row_combs)) break
    exclude_rows <- c(exclude_rows,
      as.integer(names(which.max(table(problem_row_combs)))))
  }

  for (k in 1:ncol(mat)) {
    candidate_cols <- setdiff(1:ncol(mat), exclude_cols)
    problem_col_combs <- NULL
    for (i in candidate_cols) {
      i_idx <- which(candidate_cols == i)
      for (j in candidate_cols[i_idx:length(candidate_cols)]) {
        if (sum(!is.na(mat[, i]) & !is.na(mat[, j])) <= min_shared_fields) {
          problem_col_combs <- rbind(problem_col_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_col_combs)) break
    exclude_cols <- c(exclude_cols,
      as.integer(names(which.max(table(problem_col_combs)))))
  }

  return(list('row' = exclude_rows, 'column' = exclude_cols))
}


remove_problematic_combs <- function() {
  problematic_combs <- identify_problematic_combs(
    mat = mat, min_shared_fields = min_shared_fields)
  if (!is.null(problematic_combs$row)) {
    mat <- mat[-problematic_combs$row, ]
  }
  if (!is.null(problematic_combs$column)) {
    mat <- mat[, -problematic_combs$column]
  }
  return(mat)
}
formals(remove_problematic_combs) <- formals(identify_problematic_combs)

gen_cont_colors <- function(names = c(0, 50, 5000),
                            colors = c('white', 'darkseagreen4')) {
  fun <- colorRamp2(c(0, 1), colors = colors)
  fun(seq(0, 1, length.out = length(names) + 1)) %>%
    .[-1] %>%
    setNames(as.character(names))
}
## Helper functions }}}

## {{{ General color settings
my_grey <- 'grey90'
border_line_style <- gpar(lty = 1, lwd = 1, col = 'grey20')
# plot(gen_color_vector(name = 'Zissou1', n = 5))
# row_cols <- brewer.pal(9, 'Set1')
# source('~/libs/maartenutils/R/plotting.R')
row_cols <- gen_color_vector(name = 'Royal1', arg = 9L)
analysis_name_cols <- c(gen_color_vector(name = 'Royal1', arg = 4L),
  'darkseagreen4')
column_annotation_colors = list(
  et = gen_cont_colors(c(0, 50, 5000)),
  STS = setNames('darkorange3', 'STS'),
  PR = gen_cont_colors(c(perc_ranks, NA),
    c('darkslategray4', 'white'))[seq_along(perc_ranks)],
  # VE = gen_cont_colors(c(0, 10, 100), c('white', 'goldenrod3'))
  VE = gen_cont_colors(c(0, 1, 5), c('white', 'goldenrod3'))
)
## }}}

## {{{ Constants/settings
ID_cols <- c('overlap_var', 'patient_inclusion_crit', 'LOH_HLA',
             'hla_sim_range', 'analysis_name')
# unique(dtf, by = ID_cols)
model_formula <- analysis_name + overlap_var + patient_inclusion_crit +
  LOH_HLA + hla_sim_range ~ analysis_idx
model_formula <- analysis_name + overlap_var + patient_inclusion_crit +
  LOH_HLA ~ analysis_idx

point_size <- unit(.5, 'mm')
legend_spacer <- unit(2, 'mm')
base_size <- 8
annotation_par <- gpar(fontsize = 8)
axis_par <- gpar(fontsize = 10)
base_size <- 6
annotation_par <- gpar(fontsize = 6)
axis_par <- gpar(fontsize = 8)
default_legend_settings <- list(title_gp = gpar(fontsize = base_size,
                                                fontface = 'italic'),
                                labels_gp = gpar(fontsize = base_size),
                                grid_height = unit(3, 'mm'),
                                grid_width = unit(3, 'mm'))

ca_titles = c(
  et = 'Gene expression threshold (TPM)',
  STS = 'Similarity to self filter',
  PR = 'Affinity rank percentile threshold',
  VE = 'Variant coverage threshold')
ra_titles = c(
  regression_var = 'Response variable',
  analysis_name = 'Variant selection',
  overlap_var = 'Computation of\npresentation score',
  patient_inclusion_crit = 'Removing IT\nresistant patients',
  hla_sim_range = 'Presentation score\npatient filtering',
  LOH_HLA = 'Adapting presentation score\nto LOH in HLA')
ra_titles_short = c(
  # regression_var = 'Regression variable',
  analysis_name = 'Variant selection',
  overlap_var = 'Presentation score',
  patient_inclusion_crit = 'IT-resistant patients',
  hla_sim_range = 'PS filtering',
  LOH_HLA = 'LOH in HLA')
replace_name <- function(v, old, new) {
  names(v)[names(v) == old] <- new
  return(v)
}
var_titles <- c(ca_titles, ra_titles_short) %>%
  c('project_extended' = 'Tumor type') %>%
  c('focus_allele' = 'IE focus allele') %>%
  replace_name('et', 'expression_threshold') %>%
  replace_name('VE', 'VE_threshold') %>%
  replace_name('STS', 'sts_filtering') %>%
  replace_name('PR', 'percentile_rank') %>%
  { . }
perc_ranks <- c(1, 1.9, 3, 4)
output_dir <- '~/antigenic_space/maarten-analyses/img/20-09-24/IE_heatmaps'
## }}}

## {{{ Main function
plot_pan_IE_heatmap <- function(
  tumor_type = 'pan_cancer',
  restrict_to_powered_analyses = F,
  hla_allele = 'integrated',
  fill_var = 'rc',
  return_res = F,
  output_format = 'png',
  # cell_size = unit(1.5, 'mm'),
  cell_size = unit(2, 'mm'),
  legend_size = unit(3, 'mm'),
  ann_size = unit(2, 'mm'),
  adaptive_p_val_bound = T,
  PPV = .45,
  z_normalize = F,
  ncores = 1,
  # power_var = '0.8',
  cell_size_determinant = 'p_val.adj',
  redo = F,
  subset_perc = 1,
  p_val_bound = NULL,
  debug_func = F,
  pick_allele_method = 'DT_careful',
  include_non_ds_tallies = T,
  reg_method = 'rlm_one_group',
  filter_intercept = 'none',
  filter_error = NULL,
  parse_labels = F,
  pp_hl_threshold = .01,
  compress_rc = NULL,
  cluster_vals = F,
  filter_problematic_rows = F,
  cap_fill = T,
  fn_suffix = '',
  make_filtering_diagnostic_plots = F,
  min_cell_size = .1,
  do_test_grobs = F) {

  stopifnot(is.numeric(subset_perc) && subset_perc <= 1 && subset_perc > 0)

  include_power <- cell_size_determinant == 'power' &&
    !is.null(p_val_bound) &&
    is.finite(p_val_bound)

  bayesian_analysis <- grepl('bayesian', reg_method)
  if (!bayesian_analysis && fill_var != 'yr_fractional_change') {
    message('Changing fill_var to yr_fractional_change')
    fill_var <- 'yr_fractional_change'
  }

  ## prep_pan_IE_heatmap -> compile_all_coef_overview --> summarize HLA alleles
  dtf <- prep_pan_IE_heatmap(
    p_val_bound = (p_val_bound %||% NA),
    tumor_type = tumor_type,
    hla_allele = hla_allele,
    redo = redo,
    include_rooney = T,
    fn_suffix = '',
    PPV = PPV,
    pick_allele_method = pick_allele_method,
    reg_method = reg_method,
    include_non_ds_tallies = include_non_ds_tallies,
    z_normalize = z_normalize,
    filter_error = filter_error,
    repertoire_overlap_dat = NULL,
    filter_intercept = filter_intercept,
    fill_var = fill_var,
    ncores = ncores,
    pp_hl_threshold = pp_hl_threshold,
    make_filtering_diagnostic_plots = make_filtering_diagnostic_plots,
    adaptive_p_val_bound = adaptive_p_val_bound
  )
  if (null_dat(dtf)) return(NULL)

  if (debug_func) {
    browser()
  }

  if (fn_suffix == '-im_union' && 'hla_sim_range' %in% colnames(dtf)) {
    dtf <- dtf[hla_sim_range == 'all']
  }

  prob_subs <- copy(dtf[, .N, by = analysis_grp_vars][N > 1])
  sel <- prob_subs[, N := NULL]
  if (nrow(prob_subs) > 0) {
    warning('More than one row per analysis group detected')
    setkeyv(dtf, analysis_grp_vars)
    dtf[sel[1]]
  }

  fill_mat_comp <- tryCatch(reshape2::dcast(dplyr::ungroup(dtf),
      model_formula, value.var = fill_var),
    error = function(e) { print(e); browser() }) %>%
    dplyr::mutate(regression_var = ifelse(
        analysis_name == 'rooney_param_titration',
        'zrooney', 'yield_rate')) %>%
    dplyr::mutate(analysis_name = ifelse(
        analysis_name == 'rooney_param_titration',
        'twoD_sens_analysis', as.character(analysis_name))) %>%
    dplyr::mutate(analysis_name = factor(
        analysis_name, levels = levels(dtf$analysis_name))) %>%
    dplyr::mutate(regression_var = factor(regression_var,
        levels = c('yield_rate', 'zrooney'))) %>%
    setDT

  ## Encode the various variables to be plotted in the data into separate
  ## equally sized matrices for ComplexHeatmap
  fill_mat_id <- fill_mat_comp[, sapply(fill_mat_comp, class) != 'numeric',
    with = F]
  fill_mat_X <- fill_mat_comp[, sapply(fill_mat_comp, class) == 'numeric',
    with = F] %>%
    as.matrix()

  # table(cut(as.vector(fill_mat_X), c(-Inf, 0, 1, Inf))) %>%
  #   { . / sum(.) } %>%
  #   print
  N_ids <- ncol(fill_mat_id) - 1

  p_var <- determine_p_var(fill_var)
  p_var.adj <- sprintf('%s.adj', p_var)
  signif_effect.adj <- sprintf('%s_signif_effect.adj', p_var)
  if (grepl('bayesian', reg_method)) {
    signif_effect.adj <- 'star'
  } else {
    bools <- c(p_var, p_var.adj, signif_effect.adj) %in% colnames(dtf)
    if (!all(bools)) { browser() }
  }

  icon_mat <- reshape2::dcast(dtf, model_formula,
    value.var = signif_effect.adj) %>%
    as.matrix %>%
    .[, (N_ids+1):ncol(.)]
  icon_mat[icon_mat == '*'] <- ''
  icon_freqs <- table(icon_mat[!is.na(icon_mat)])
  icon_freqs <- setNames(c(icon_freqs['L'], icon_freqs['G']), c('L', 'G'))
  icon_freqs <- repl.na(icon_freqs)

  cell_size_mat <- reshape2::dcast(dtf, model_formula,
    value.var = cell_size_determinant) %>%
    as.matrix %>%
    .[, (N_ids+1):ncol(.)] %>%
    apply(2, as.numeric)
  # hist(as.vector(cell_size_mat), breaks = 100)

  ## Transform cell sizes such that bigger is interpretable
  ## as a stronger effect
  if (grepl('power', cell_size_determinant) &&
      power_var %in% colnames(dtf)) {
    cell_size_mat[cell_size_mat > 1] <- NA
    messagef('%s of cells annotated with power analysis data',
             scales::percent(1 - mean(is.na(as.vector(cell_size_mat)))))
    cell_size_mat <- 1 / cell_size_mat
  } else if (grepl('p_val', cell_size_determinant)) {
    cell_size_mat <- -log10(cell_size_mat)
  } else if (grepl('log2_est_error', cell_size_determinant)) {
    cell_size_mat <- -cell_size_mat
  } else {
    message('No cells annotated with power analysis data')
    cell_size_mat <- matrix(rep(NA, prod(dim(icon_mat))),
                            ncol = ncol(icon_mat))
  }
  n_cell_size_mat <- trans_size_size(cell_size_mat)

  ## NAs are okay, Infs, -Inf, and NaN not
  if (any(!is.finite(cell_size_mat) & !is.na(cell_size_mat))) {
    idxs <- which(!is.finite(cell_size_mat) &
                  !is.na(cell_size_mat), arr.ind = T)
    cell_size_mat[idxs[1, ]]
    perm_browser()
  }

  if (restrict_to_powered_analyses) {
    fill_mat_X[which(is.na(cell_size_mat), arr.ind = T)] <- NA
    icon_mat[which(is.na(cell_size_mat), arr.ind = T)] <- NA
    n_cell_size_mat[which(is.na(cell_size_mat), arr.ind = T)] <- NA
  }

  informative_cell_idx <-
    which(!is.na(fill_mat_X) & !eps(fill_mat_X, 0, 1e-12), arr.ind = T)
  if (fill_var == 'rc' && length(informative_cell_idx) == 0) {
    return(NULL)
  }

  ## Subselect informative (i.e. not all NA) rows and columns
  if (filter_problematic_rows) {
    problematic_combs <- identify_problematic_combs(fill_mat_X)
    i_rows <- setdiff(1:nrow(fill_mat_X), problematic_combs$row)
    i_rows_alt <- tryCatch(apply(fill_mat_X, 1, function(x) !all(is.na(x))),
      error = function(e) { print('Could not subset rows'); 1:nrow(fill_mat_X) })
    ## Subselect a percentage of rows
    i_rows <- i_rows[1:ceiling(length(i_rows) * subset_perc)]
    i_columns <- setdiff(1:ncol(fill_mat_X), problematic_combs$column)
    # i_columns <- tryCatch(apply(fill_mat_X, 2, function(x) !all(is.na(x))),
    #   error = function(e) { print('Could not subset columns');
    #   1:col(fill_mat_X) })
    ## Subselect a percentage of columns
    i_columns <- i_columns[1:ceiling(length(i_columns) * subset_perc)]
  } else if (T) {
    i_columns <- which(apply(fill_mat_X, 2, function(x) any(!is.na(x))))
    i_rows <- which(apply(fill_mat_X, 1, function(x) any(!is.na(x))))
  } else {
    i_columns <- 1:ncol(fill_mat_X)
    i_rows <- 1:nrow(fill_mat_X)
  }
  nc <- length(i_columns)
  nr <- length(i_rows)
  fill_mat_X <- fill_mat_X[i_rows, i_columns]
  cell_size_mat <- cell_size_mat[i_rows, i_columns]
  n_cell_size_mat <- n_cell_size_mat[i_rows, i_columns]
  icon_mat <- icon_mat[i_rows, i_columns]
  fill_mat_id <- fill_mat_id[i_rows, ]
  yr_idxs <- which(fill_mat_id$regression_var == 'yield_rate')
  rooney_idxs <- which(fill_mat_id$regression_var == 'zrooney')

  summary(as.vector(fill_mat_X), na.rm = T)

  if (cap_fill) {
    ## Retrieve the percentage of analyses with more than 100% neo-antigen
    ## enrichment (i.e. growth of 100%). These analyses are not to be trusted
    pos_perc <- as.vector(fill_mat_X) %>%
      { .[!is.na(.)] } %>%
      { ecdf(.)(1) }
    ## Cap extremely low and high values
    probs <- sort(c(pos_perc, 1 - pos_perc))
    cap_values <- quantile(as.vector(fill_mat_X), probs = probs, na.rm = T)
    fill_mat_X[which(fill_mat_X <= cap_values[1], arr.ind = T)] <-
      cap_values[1]
    fill_mat_X[which(fill_mat_X >= cap_values[2], arr.ind = T)] <-
      cap_values[2]
    stopifnot(cap_values == range(as.vector(fill_mat_X), na.rm = T))
  }

  X_max <- max(abs(fill_mat_X), na.rm = T)
  X_max_rooney <- max(abs(fill_mat_X[-yr_idxs, ]), na.rm = T)

  if (!is.null(compress_rc) && !is.na(compress_rc) && compress_rc != 0) {
    ## Transform to compressed space with sigmoidal function
    ## x_mid is the x-coordinate for which y == .5, the smaller it is, the more
    ## the input values get compressed
    x_mid <- 0
    # compress_rc = 1
    trans_to <- function(x, a = 2^compress_rc, x_mid = 0,
                         C = x_max, D = 2 * x_max) {
      (1 / (1 + exp(-a * x + x_mid)) * D) - C
    }
    if (F) {
      # x_min = 6
      # x_mid = 1
      # x_mid = 0
      with(list(x = seq(-x_min, x_min, by = 1e-5)),
        plot(x, trans_to(x, a = 1, x_mid), type = 'l',
          axes = F, xlab = 'Input', ylab = 'Output'))
      segments(c(-15, x_mid), c(.5, -5), c(x_mid, x_mid),
        c(.5, .5), col = 'red')
      p_eps <- 0.09
      points(x = x_mid, y = .5, pch = 19, col = 'red')
      text(x = x_mid + p_eps, y = .5 + p_eps, labels = 'x_mid', col = 'red',
        adj = c(0, 1))
      # axis(side = 1,
      #   at = c(seq(-x_min, x_min, by = 1), x_mid, x_max),
      #   labels = c(seq(-x_min, x_min, by = 1), 'x_mid', 'x_max'))
      axis(side = 1)
      axis(side = 2)
    }
  } else {
    trans_to <- trans_back <- function(x) x
  }

  X_max_norm <- round_up(trans_to(X_max), dec = 1)

  if (F) {
    fill_breaks <- X_max_norm %>% { seq(-1 * ., ., length.out = 5) }
  } else {
    obs_range <- range(fill_mat_X, na.rm = T)
    fill_breaks <- c(seq(obs_range[1], obs_range[2], length.out = 5), 0) %>%
      sort
  }
  if (do_test_grobs) hist(as.vector(fill_mat_X), breaks = 100)
  rc_labels_yr <- trans_back(fill_breaks)

  if (parse_labels) {
    rc_labels_yr %<>% fancy_scientific(rc_labels_yr)
    rc_labels_yr %<>% sapply(function(x) parse(text = x))
  } else {
    rc_labels_yr <- round(fill_breaks, 2)
  }

  ftr <- 1.8
  ftr <- 2.5
  ## Heat map specific color settings
  rc_col_fun_yr <- unname(X_max_norm) %>%
    { c(-1 * ., 0, .) } %>%
    { colorRamp2(colors = c(
        darken('#A30D1D', ftr), my_grey, darken('#3F9CB5', ftr)),
                 breaks = .) }

  rc_res_rooney <- ceiling(abs(log10(X_max_rooney)))
  if (F) {
    X_max_norm <- round_up(X_max_rooney, dec = rc_res_rooney - 1)
    rc_breaks_rooney <- trans_to(X_max_norm) %>%
      seq(-1 * ., ., length.out = 5)
    rc_labels_rooney <- trans_back(rc_breaks_rooney)
  } else {
    X_max_norm <- round_up(X_max_rooney, dec = rc_res_rooney)
    ## We want between 5 and 9 markers, more markers if compression is high
    rc_labels_rooney <- seq(-X_max_norm, X_max_norm,
                            length.out = min(max(5, compress_rc), 9)) %>%
      round(rc_res_rooney + 1)
    rc_breaks_rooney <- trans_to(rc_labels_rooney)
  }
  rc_col_fun_rooney <- trans_to(X_max_rooney) %>%
    setNames(NULL) %>%
    { c(-1 * ., 0, .) } %>%
    { colorRamp2(colors = c('#A30D1D', my_grey, '#3F9CB5'), breaks = .) }

  LOHHLA_vec <- fill_mat_id$LOH_HLA %>% { ifelse(grepl('no_', .), NA, .) }
  LOHHLA_levels <- setdiff(unique(LOHHLA_vec), NA)

  LOHHLA_cols <- darken(rep(row_cols[7], 2), c(1, 1.5)) %>%
    setNames(c('LOHHLA', 'strict_LOHHLA'))
  LOHHLA_cols <- LOHHLA_cols[LOHHLA_levels]
  # analysis_name_cols %>% as.color_vector %>% plot
  row_annotation_colors <- list(
    regression_var = setNames(row_cols[1:2], c('yield_rate', 'zrooney')),
    analysis_name = setNames(analysis_name_cols, analysis_names),
    LOH_HLA = LOHHLA_cols,
    patient_inclusion_crit = setNames(darken(row_cols[8], 1), 'TR'),
    overlap_var = setNames(darken(row_cols[9], 1), c('mean_score_AB')),
    hla_sim_range = setNames(darken(row_cols[1], 1), c('<= 1'))
  )

  row_annotation_settings <- list(
    regression_var = list(
      title = ra_titles['regression_var'],
      breaks = c('yield_rate', 'zrooney'),
      labels = c('Neo-antigen yield rate', 'Obs. / exp.\nneo-antigens'),
      values = fill_mat_id$regression_var),
    analysis_name = list(
      title = ra_titles['analysis_name'],
      breaks = analysis_names[c(1, 3, 4, 5)],
      labels = c('All variants', 'Clonal',
        'No drivers/essentials', 'Marty oncogenic'),
      values = fill_mat_id$analysis_name),
    LOH_HLA = list(
      title = ra_titles['LOH_HLA'],
      legend_direction = 'horizontal',
      nrow = 2,
      breaks = LOHHLA_levels,
      labels = setNames(c('LOHHLA', 'High confidence LOHHLA'),
        c('LOHHLA', 'strict_LOHHLA'))[LOHHLA_levels],
      values = LOHHLA_vec),
    patient_inclusion_crit = list(
      title = ra_titles['patient_inclusion_crit'],
      values = fill_mat_id$patient_inclusion_crit %>%
      { ifelse(. == 'TR', ., NA) }),
    overlap_var = list(
      title = ra_titles['overlap_var'],
      breaks = 'mean_score_AB',
      labels = 'Restricted to A & B alleles',
      values = fill_mat_id$overlap_var %>%
      { ifelse(. == 'mean_score_AB', ., NA) }),
    hla_sim_range = list(
      title = ra_titles['hla_sim_range'],
      breaks = '<= 1',
      labels = c('No pts \\w PS > 1'),
      values = fill_mat_id$hla_sim_range %>%
      { ifelse(. == '<= 1', ., NA) })
  )

  lgds_to_make <- setdiff(names(row_annotation_settings),
                          c('regression_var'))
  if (grepl('im_union', fn_suffix)) {
    lgds_to_make <- setdiff(lgds_to_make, 'hla_sim_range')
  }
  lgds_to_make %<>%
    { .[map_lgl(., ~ !all(is.na(row_annotation_settings[[.x]]$values)))] }
  row_annotation_settings %<>%
    { .[c('regression_var', lgds_to_make)] }

  row_annotation_legends <-
    map(lgds_to_make, function(an) {
      at_vals <- row_annotation_settings[[an]][['breaks']] %||%
        unique(row_annotation_settings[[an]][['values']])
      if (is.logical(at_vals)) at_vals <- setdiff(at_vals, FALSE)
      at_vals <- setdiff(at_vals, NA)
      label_vals <- at_vals
      if (!is.null(row_annotation_settings[[an]][['labels']]))
        label_vals <- row_annotation_settings[[an]][['labels']]
      cols <- row_annotation_colors[[an]][at_vals]
      if (length(at_vals) != length(cols)) {
        print('Not enough colors')
      }
      lgd <- tryCatch(Legend(
        at = at_vals, labels = label_vals,
        legend_gp = gpar(fill = cols),
        title = row_annotation_settings[[an]][['title']],
        nrow = row_annotation_settings[[an]][['nrow']],
        direction = row_annotation_settings[[an]][['legend_direction']],
        grid_width = legend_size,
        grid_height = legend_size,
        title_gp = default_legend_settings$title_gp,
        labels_gp = default_legend_settings$labels_gp),
      error = function(e) { print(e); NULL })
      return(lgd)
    })
  # test_grob(row_annotation_legends[[1]], do_test_grobs)
  # test_grob(row_annotation_legends[[2]], do_test_grobs)
  # test_grob(row_annotation_legends[[3]], do_test_grobs)
  # test_grob(row_annotation_legends[[4]], do_test_grobs)
  # test_grob(row_annotation_legends[[5]], do_test_grobs)

  row_annotation_args <- list(
    # regression_var = row_annotation_settings$regression_var$values,
    # regression_var = anno_block(gp = gpar(fill = row_cols[1:2]),
    #     labels = rev(row_annotation_settings$regression_var$labels),
    #     labels_gp = gpar(col = 'white', fontsize = 8),
    #     which = 'row'),
    analysis_name = row_annotation_settings$analysis_name$values,
    LOH_HLA = row_annotation_settings$LOH_HLA$values,
    patient_inclusion_crit = row_annotation_settings$patient_inclusion_crit$values,
    overlap_var = row_annotation_settings$overlap_var$values,
    which = 'row',
    na_col = 'white',
    gp = gpar(width = ann_size, lwd = 0, col = 'white'),
    col = row_annotation_colors,
    # simple_anno_size = ann_size,
    # annotation_width = c(40, 1, 1, 1, 1, 1) * ann_size,
    # annotation_width = c(1, 1, 1, 1, 1) * ann_size,
    width = 4 * ann_size,
    # height = ann_size,
    simple_anno_size = ann_size,
    gap = unit(c(0, 0, 0, 0, 0), 'mm'),
    # row_split = fill_mat_id$regression_var,
    # left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = row_cols[1:2]),
    #     labels = rev(row_annotation_settings$regression_var$labels),
    #     labels_gp = gpar(col = 'white', fontsize = 8))),
    show_annotation_name = F)
  if (!grepl('-im_union', fn_suffix)) {
    row_annotation_args[['hla_sim_range']] <-
      row_annotation_settings$hla_sim_range$values
  }
  # lt = lapply(1:3, function(x) cumprod(1 + runif(1000, -x/100, x/100)) - 1)
  # ha = rowAnnotation(foo = anno_horizon(lt))
  # ha <- HeatmapAnnotation(
  #   regression_var =
  #     anno_block(gp = gpar(fill = row_cols[1:2]),
  #                labels = row_annotation_settings$regression_var$labels,
  #                labels_gp = gpar(col = 'white', fontsize = 8),
  #                ))
  row_annotation <- do.call(HeatmapAnnotation, row_annotation_args)
  # test_grob(row_annotation, do_test_grobs)

  column_annotation_legend_settings <-
    with(as.list(ds_param_grid[as.integer(colnames(fill_mat_X))[i_columns], ]), {
      list(
        PR = list(
          at = unique(percentile_rank),
          # at = rev(unique(percentile_rank)),
          title = 'Affinity rank percentile threshold',
          labels = NULL,
          nrow = 1,
          legend_direction = 'horizontal'),
        et = list(
          at = setdiff(unique(expression_threshold), NA),
          title = 'Gene expression threshold',
          nrow = 1,
          labels = NULL,
          legend_direction = 'horizontal'),
        VE = list(
          at = setdiff(unique(VE_threshold), NA),
          title = 'Variant coverage threshold',
          labels = NULL,
          nrow = 1,
          legend_direction = 'horizontal'),
        STS = list(
          at = unique(sts_filtering),
          title = 'Similarity to self filter',
          nrow = 1,
          labels = 'STS',
          # breaks = 'STS',
          labels = 'Neo-antigens filtered for STS'))
  })

  column_annotation_legends <-
    # lapply(names(column_annotation_legend_settings)[1],
    lapply(names(column_annotation_legend_settings),
      function(an) {
        at_vals <- column_annotation_legend_settings[[an]][['at']]
        if (is.logical(at_vals))
          at_vals <- setdiff(at_vals, FALSE)
        label_vals <- at_vals
        if (!is.null(column_annotation_legend_settings[[an]][['labels']]))
          label_vals <- column_annotation_legend_settings[[an]][['labels']]
        lgd <- Legend(
          at = at_vals,
          labels = label_vals,
          legend_gp = gpar(fill = column_annotation_colors[[an]]),
          title = column_annotation_legend_settings[[an]][['title']],
          nrow = column_annotation_legend_settings[[an]][['nrow']],
          direction = column_annotation_legend_settings[[an]][['legend_direction']],
          grid_width = legend_size,
          grid_height = legend_size,
          title_gp = default_legend_settings$title_gp,
          labels_gp = default_legend_settings$labels_gp)
        test_grob(lgd, do_test_grobs)
        return(lgd)
      })
  test_grob(column_annotation_legends[[1]], do_test_grobs)
  test_grob(column_annotation_legends[[2]], do_test_grobs)
  test_grob(column_annotation_legends[[3]], do_test_grobs)
  test_grob(column_annotation_legends[[4]], do_test_grobs)

  column_annotation <-
    with(as.list(ds_param_grid[as.integer(colnames(fill_mat_X))[i_columns], ]),
      HeatmapAnnotation(
        PR = percentile_rank,
        et = expression_threshold,
        VE = VE_threshold,
        STS = ifelse(sts_filtering, 'STS', NA),
        annotation_name_gp = annotation_par,
        show_annotation_name = F,
        na_col = 'white',
        which = 'column',
        height = ann_size * 4,
        # width = cell_size * length(percentile_rank),
        # gap = unit(rep(0, 4), 'mm'),
        col = column_annotation_colors,
        simple_anno_size = ann_size,
        gp = gpar(lwd = 0, col = 'white'),
        annotation_name_side = 'right'
      )
    )
  test_grob(column_annotation, doit = do_test_grobs)

  main_heatmap <- ComplexHeatmap::Heatmap(
    fill_mat_X,
    name = 'rc',
    show_row_names = F,
    column_title = 'Neo-antigen prediction pipeline settings',
    column_title_side = 'bottom',
    show_column_names = F,
    cluster_rows = cluster_vals,
    cluster_columns = cluster_vals,
    clustering_method_columns = 'ward.D2',
    clustering_method_rows = 'ward.D2',
    column_title_gp = axis_par,
    cell_fun = function(j, i, x, y, w, h, col) {
      if (cell_size_determinant == 'none') {
        cs <- 1
      } else {
        cs <- n_cell_size_mat[i, j]
      }
      if (!is.na(cs) && !is.na(fill_mat_X[i, j])) {
        ## Sqrt transform the cell size indicator in order to get the required
        ## side length (of a square). I.e. surface area will be linearly
        ## proportional to cs
        w <- h <- cell_size * cs
        fill_col <- rc_col_fun_yr(trans_to(fill_mat_X[i, j]))
        if (F && fill_mat_id[i, 'regression_var'] == 'zrooney') {
          fill_col <- rc_col_fun_rooney(trans_to(fill_mat_X[i, j]))
        }
        grid.rect(x = x, y = y,
                  width = w, height = h,
                  gp = gpar(lwd = 0, fill = fill_col, col = 'white'))
      }
      if (!is.na(icon_mat[i, j]) && icon_mat[i, j] != '') {
        if (icon_mat[i, j] == 'L') {
          l_pch = 6
        } else if (icon_mat[i, j] == 'G') {
          l_pch = 2
        }
        grid.points(x, y, pch = l_pch, size = point_size)
      }
    },
    rect_gp = gpar(type = 'none'),
    bottom_annotation = column_annotation,
    show_column_dend = F,
    # row_title = row_annotation_settings$regression_var$labels,
    # row_title = '',
    row_title_gp = annotation_par,
    show_row_dend = F)
  if (do_test_grobs)
    print(main_heatmap)

  ## Legend
  if (cell_size_determinant != 'none') {
    if (grepl('p_val', cell_size_determinant)) {
      ## From weak to strong evidence against H0
      ## Compute locations in log10-space, label them with coords in linear
      ## space
      # table(10^(-1 * sort(as.vector(cell_size_mat))))
      size_legend_breaks <- range(as.vector(cell_size_mat), na.rm = T) %>%
        { seq(.[1], .[2], length.out = 5) } %>%
        { 10^(-1 * .) } %>%
        { . }
      if (parse_labels) {
        size_legend_labels <- size_legend_breaks %>%
          fancy_scientific(digits = 2) %>%
          sapply(function(x) parse(text = x))
      } else {
        size_legend_labels <- format(size_legend_breaks, digits = 2)
      }
      if (grepl('adj', cell_size_determinant)) {
        size_title <- glue('Adjusted p-value slope')
      } else {
        size_title <- glue('Unadjusted p-value slope')
      }
    } else if (grepl('est_error', cell_size_determinant)) {
      size_legend_breaks <- range(as.vector(cell_size_mat), na.rm = T) %>%
        { seq(.[1], .[2], length.out = 5) } %>%
        { . }
      if (parse_labels) {
        size_legend_labels <- size_legend_breaks %>%
          fancy_scientific(digits = 2) %>%
          sapply(function(x) parse(text = x))
      } else {
        size_legend_labels <- format(size_legend_breaks, digits = 2)
      }
      size_title <- as.character(glue('-log2(estimation error)'))
    } else {
      size_legend_breaks <- range(cell_size_mat, na.rm = T) %>%
        { . + 1 } %>%
        log10 %>%
        { seq(.[1], .[2], length.out = 5) } %>%
        { 10^. - 1 } %>%
        { -1 * . } %>%
        round(3)
      # c(0, fancy_scientific(size_legend_breaks[2:length(power_legend_breaks)],
      #   digits = 2)),
      # c(0, round(size_legend_breaks[2:length(power_legend_breaks)], 1)),
      size_title <- glue('IE detection power\n(required fractional\n',
                         'loss for sufficient\nstatistical power)')
    }
    sizes <- cell_size * seq(0, 1, length.out = 5) %>%
      trans_size_size
    cell_size_lgd <- tryCatch(Legend(
      at = size_legend_breaks,
      legend_gp = list(col = 'grey80'),
      size = sizes,
      pch = 15,
      type = 'points',
      grid_width = cell_size,
      grid_height = cell_size,
      title_gp = default_legend_settings$title_gp,
      labels_gp = default_legend_settings$labels_gp,
      # row_gap = .7 * cell_size,
      # cell_scaling_factors = sizes,
      labels = size_legend_labels,
      title = size_title
    ), error = function(e) { print(e); browser() })
    test_grob(cell_size_lgd, do_test_grobs)
    title_grob <- cell_size_lgd@grob$children[[1]]$label %>%
      textGrob(x = unit(0, 'npc'),
               y = unit(0, 'npc') + unit(2, 'mm'),
               hjust = 0, vjust = 0,
               gp = gpar(fontsize = 6, fontface = 'italic')
      )
    bottom_grob <- cell_size_lgd@grob$children[[2]]
    cell_size_lgd_grob <- gtable::gtable(
      widths = unit(c(convertUnit(widthDetails(title_grob), 'cm')),
                    c('cm')),
      heights = unit(c(1, convertUnit(heightDetails(bottom_grob), 'cm')),
                     c('lines', 'cm'))
      ) %>%
      gtable::gtable_add_grob(grobs = title_grob, t = 1,
        b = 1, l = 1, r = 1, clip = 'off') %>%
      gtable::gtable_add_grob(grobs = bottom_grob,
        t = 2, b = 2, l = 1, r = 1, clip = 'off')
    test_grob(cell_size_lgd_grob, doit = do_test_grobs)
  } else {
    cell_size_lgd <- NULL
  }

  if (parse_labels) {
    icon_labels <- mapply(sprintf, rep('"%s"~(italic(n)==%d)', 2),
      c('Negative slope', 'Positive slope'), icon_freqs[c('L', 'G')]) %>%
      setNames(NULL) %>%
      { . }
    icon_labels %<>% sapply(function(x) parse(text = x))
  } else {
    icon_labels <- mapply(sprintf, rep('%s (n=%d)', 2),
      c('Negative slope', 'Positive slope'), icon_freqs[c('L', 'G')]) %>%
      unname %>% { . }
  }

  if (!bayesian_analysis) {
    size_title <- 'Significant\nregression slope\n(FDR-corrected)'
  } else {
    size_title <- glue::glue('Posterior probability > \\
      {format(1 - pp_hl_threshold, digits = 2)}')
  }
  icon_lgd <- tryCatch(Legend(
    at = c('L', 'G'),
    labels = icon_labels,
    title = size_title,
    grid_height = legend_size,
    grid_width = legend_size,
    size = point_size,
    # legend_gp = c(#default_legend_settings,
    #   list(fill ='white', fontsize = .5,
    #        width = ann_size, lwd = 1, col = 'black')),
    type = 'points',
    ncol = 1,
    title_gp = default_legend_settings$title_gp,
    labels_gp = default_legend_settings$labels_gp,
    pch = rev(c(2, 6))
  ), error = function(e) { print(e); NULL })
  test_grob(icon_lgd, doit=do_test_grobs)

  # rc_lgd <- tryCatch(Legend(
  #   at = rc_breaks,
  #   legend_gp = c(list(fill = rc_col_fun(rc_breaks))),
  #   grid_height = legend_size,
  #   grid_width = legend_size,
  #   title_gp = default_legend_settings$title_gp,
  #   labels_gp = default_legend_settings$labels_gp,
  #   labels = rc_labels,
  #   title = 'Regression coefficient'
  # ), error = function(e) { print(e); browser() })
  # test_grob(rc_lgd, do_test_grobs)

  fill_lgd_yr <- tryCatch(Legend(
    at = fill_breaks,
    col_fun = rc_col_fun_yr,
    # legend_gp = c(list(fill = rc_col_fun(rc_breaks))),
    grid_height = legend_size,
    grid_width = legend_size,
    title_gp = default_legend_settings$title_gp,
    labels_gp = default_legend_settings$labels_gp,
    labels = rc_labels_yr,
    title = switch(fill_var,
      'rc' = 'Neo-antigen\nyield rate\nregression slope',
      'effective_editing_max' = 'FC neo-antigen yield rate\nbetween PS=0 and PS=1',
      'effective_editing_mean' = 'Neo-antigen loss\nwith PS in range\n[0, avg(PS)]',
      'yr_fractional_change' = 'relative change\nNAYR between\nPS=0 and PS=1',
      'depletion_full' = 'log2 neo-antigen\ndepletion between\nPS=0 and PS=1',
      'estimate' = 'log2 neo-antigen\ndepletion between\nPS=0 and PS=1',
      'rc_CYT' = 'Regression coefficient (CYT)')
  ), error = function(e) { print(e); browser() })
  test_grob(fill_lgd_yr, do_test_grobs)

  title_grob <- fill_lgd_yr@grob$children[[1]]$label %>%
    textGrob(x = unit(0, 'npc'),
             y = unit(0, 'npc') + unit(2, 'mm'),
             hjust = 0, vjust = 0,
             gp = gpar(fontsize = 6, fontface = 'italic')
  )
  bottom_grob <- fill_lgd_yr@grob$children[[2]]

  hist_grob <- local({
    p_dat <-
      data.frame(x = suppressWarnings(trans_to(as.vector(fill_mat_X)))) %>%
      mutate(sign_var = c('L', 'G')[as.integer(x > 0) + 1])
    # sum(maartenutils::eps(p_dat$x, 0), na.rm = T)
    # summary(as.vector(fill_mat_X))
    # summary(trans_to(as.vector(fill_mat_X)))
    # summary(p_dat$x)
    if (do_test_grobs) hist(p_dat$x, breaks = 100)
    # min(p_dat$x, na.rm = T)
    # format(fill_breaks, digits = 5)
    p <- ggplot(p_dat, aes(x = x, fill = NULL)) +
      geom_histogram(binwidth = 2 * diff(range(fill_breaks)) / sum(is.na(p_dat$x))) +
      geom_vline(xintercept = 0, size = .5, linetype = 3, color = 'grey98') + 
      scale_x_continuous(name = '', expand = c(0.0, 0.0),
                         # breaks = fill_breaks,
                         # limits = range(fill_breaks) + 0 * 1e-1 * c(-1, 1)) +
                         limits = range(fill_breaks)) +
      scale_y_continuous(name = '', expand = c(0, 0),
        # trans = 'identity') +
        trans = 'log') +
      # scale_fill_manual(values = rev(c('#A30D1D', '#3F9CB5'))) +
      theme_bw() +
      theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
            panel.spacing = unit(0, 'cm'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank())
    library(gtable)
    p_gtab <- maartenutils::to_g(p)
    ## Extract the panel, stripping away the axes and margins
    ret <- p_gtab$grob[[which(p_gtab$layout$name == 'panel')]]
    ret$vp <- viewport(angle = 90, height = unit(1, 'cm'),
                       width = heightDetails(bottom_grob))
    ret
  })
  # test_grob(hist_grob, do_test_grobs)

  rc_lgd_yr_grob <- gtable(
      widths = unit(
        c(1, convertUnit(widthDetails(title_grob), 'cm')),
        c('cm', 'cm')),
      heights = unit(
        c(convertUnit(heightDetails(title_grob), 'cm'), 
          convertUnit(heightDetails(bottom_grob), 'cm')), 
        c('cm', 'cm'))
    ) %>%
    gtable::gtable_add_grob(grobs = title_grob, t = 1,
      l = 1, r = 2, clip = 'off') %>%
    gtable::gtable_add_grob(grobs = bottom_grob, t = 2,
      l = 2, clip = 'off') %>%
    gtable::gtable_add_grob(grobs = hist_grob, t = 2,
      l = 1, clip = 'on')
  test_grob(rc_lgd_yr_grob, do_test_grobs)

  if (F) {
    rc_lgd_rooney <- tryCatch(Legend(
      at = rc_breaks_rooney,
      col_fun = rc_col_fun_rooney,
      # legend_gp = c(list(fill = rc_col_fun(rc_breaks))),
      grid_height = legend_size,
      grid_width = legend_size,
      title_gp = default_legend_settings$title_gp,
      labels_gp = default_legend_settings$labels_gp,
      labels = rc_labels_rooney,
      title = switch(fill_var,
        'rc' = 'Observed / expected\nneo-antigens\nregression slope',
        'effective_editing_max' = 'Neo-antigen loss\nwith PS in range\n[0, 1]',
        'effective_editing_max' = 'Observed / expected\nneo-antigens\nloss in PS-range [0, 1]',
        'effective_editing_mean' = 'Neo-antigen loss\nwith PS in range\n[0, avg(PS)]',
        'rc_CYT' = 'Regression coefficient (CYT)')
    ), error = function(e) { print(e); browser() })
    test_grob(rc_lgd_rooney, do_test_grobs)
  }

  lgds <- c(list(rc_lgd_yr_grob, cell_size_lgd_grob), row_annotation_legends)
  if (sum(icon_freqs) > 0) {
    lgds <- append(lgds, icon_lgd, 2)
  }
  ra_legends_grob <-
    gtable_matrix(name = 'vlegends',
      grobs = matrix(map(lgds, ~just_left(extract_legend_grob(.x))),
                     ncol = 1),
      heights = purrr::map_dbl(lgds,
                               ~convertUnit(heightDetails(extract_legend_grob(.x)),
                                            'cm')) %>%
        unit('cm'),
      widths = purrr::map_dbl(lgds,
                                ~convertUnit(widthDetails(extract_legend_grob(.x)),
                                             'cm')) %>%
        max %>% unit('cm')) %>%
    gtable::gtable_add_row_space(unit(2, 'mm'))
  test_grob(ra_legends_grob, do_test_grobs)

  ca_legends_grob <- do.call(packLegend, c(
      # list(legend_header('Neo-antigen prediction\npipeline settings')),
      column_annotation_legends,
      list(direction = 'horizontal')))
  test_grob(ca_legends_grob, do_test_grobs)

  comb <- row_annotation + main_heatmap

  annotation_par <- c(annotation_par, list(fontface = 'italic'))
  annotation_width <- ca_titles %>%
    map(~textGrob(.x, gp = do.call(gpar, annotation_par))) %>%
    map(~grobWidth(.x)) %>%
    map_dbl(~convertUnit(.x, 'cm')) %>%
    max
  ra_label_ext <- unit(15, 'mm')

  top_left_height <- cell_size * (nr / 1) + ra_label_ext
  top_height <- max(top_left_height, heightDetails(ra_legends_grob))
  plot_ca_names <-
    as.numeric(convertUnit(top_left_height, 'cm')) >
    (as.numeric(convertUnit(heightDetails(ra_legends_grob), 'cm')) + 2)
  top_left_width <- cell_size * (nc / 1) +
    ComplexHeatmap:::width(row_annotation)
  top_right_width <- do.call(max,
      list(#convertUnit(ComplexHeatmap:::width(rc_lgd_rooney), 'cm'),
           convertUnit(ComplexHeatmap:::width(fill_lgd_yr), 'cm'),
           unit(annotation_width, 'cm'))) + unit(10, 'mm')
  bottom_left_height <-
    # ComplexHeatmap:::height(ra_legends_grob) +
    ComplexHeatmap:::height(ca_legends_grob) + unit(0, 'cm')
  # top_left_width <- unit(1, 'npc') - top_right_width
  # top_left_height <- unit(1, 'npc') - bottom_left_height

  plot_height <- {
      top_height + bottom_left_height +
      ComplexHeatmap:::height(column_annotation) + unit(0, 'cm')
    } %>%
    min(unit(29, 'cm')) %>%
    max(unit(10, 'cm'))

  ## Plot width cannot be more than 17.4
  plot_width <- { cell_size * nc + unit(7.2, 'cm') } %>%
    min(unit(17.4, 'cm'))

  messagef('Plot dimensions %.1f [cm] x %.1f [cm]',
    as.numeric(convertUnit(plot_height, 'cm')),
    as.numeric(convertUnit(plot_width, 'cm')))

  plot_expr <- rlang::expr({
    vp_grid <-
      grid.layout(nrow = 2, ncol = 2,
                  widths = unit.c(top_left_width, top_right_width),
                  heights = unit.c(top_height, bottom_left_height))
    # grid.show.layout(vp_grid)
    pushViewport(viewport(layout = vp_grid))

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    tryCatch(
      draw(comb,
           row_title = 'Immunoediting analysis settings',
           row_split = fill_mat_id$regression_var,
           row_title_gp = axis_par,
           padding = unit(c(0, 0, 0, 0), 'mm'), ## bottom, left, top, right
           newpage = F,
           show_annotation_legend = F,
           show_heatmap_legend = F,
           merge_legend = T),
      error = function(e) { print(e) })
    if (plot_ca_names && nrow(fill_mat_id) >= 60) {
      for (an in names(ca_titles)) {
        decorate_annotation(an, {
          grid.text(ca_titles[an],
            x = unit(1, 'npc') + legend_spacer, just = 'left',
            gp = do.call(gpar, c(annotation_par, list(fontface = 'italic'))))
        })
      }
    }
    for (an in lgds_to_make) {
      decorate_annotation(an, {
        grid.text(ra_titles_short[an],
          y = unit(-2, 'mm'), just = 'right',
          gp = do.call(gpar,
            c(annotation_par, list(fontface = 'italic'))),
          rot = 90)
      }, slice = 2)
    }
    for (i in 1:2) {
      decorate_row_title('rc', {
        grid.rect(gp = gpar(lwd = 0, fill = row_cols[i]))
        grid.text(row_annotation_settings$regression_var$labels[i],
          x = unit(.5, 'npc'),
          y = unit(.5, 'npc'),
          just = 'center',
          rot = 90,
          gp = do.call(gpar,
            c(annotation_par, list(col = 'white'))))
      }, slice = i)
    }
    decorate_row_title('rc', {
      grid.text('Response variable',
                y = unit(-2, 'mm'), just = 'right',
                gp = do.call(gpar,
                             c(annotation_par, list(fontface = 'italic'))),
                rot = 90)
    }, slice = 2)
    decorate_heatmap_body('rc', {
      grid.lines(c(0.0, 0.0, 1, 1, 0), c(0, 1, 1, 0, 0),
                 gp = border_line_style)
    }, slice = 1)
    decorate_heatmap_body('rc', {
      grid.lines(c(0.0, 0.0, 1, 1, 0), c(0, 1, 1, 0, 0),
                 gp = border_line_style)
    }, slice = 2)
    popViewport()

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    viewport(just = c('left', 'top'),
             width = widthDetails(ra_legends_grob),
             height = heightDetails(ra_legends_grob),
             x = unit(0, 'npc') + legend_spacer,
             gp = gpar(fontsize = 6),
             y = unit(1, 'npc')) %>% pushViewport()
    grid.draw(ra_legends_grob)
    popViewport()
    popViewport()

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = c(1, 2)))
    draw(ca_legends_grob,
         just = c('center', 'top'),
         x = unit(0.5, 'npc'),
         y = unit(1, 'npc') - unit(0, 'mm'))
    popViewport()

    popViewport()
  })

  if (return_res) {
    op <- grid.grabExpr(eval(plot_expr), wrap = T, wrap.grobs = T,
                        width = convertUnit(plot_width, 'inch'),
                        height = convertUnit(plot_height, 'inch'))
    return(op)
  } else {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = T)
    tumor_type_simple <- tumor_types[tumor_type]
    of_base_name <- file.path(output_dir, glue::glue('{tumor_type_simple}\\
      {make_flag(hla_allele)}\\
      {make_flag(compress_rc)}\\
      {make_flag(fill_var)}\\
      {make_flag(z_normalize)}\\
      {make_flag(pick_allele_method)}\\
      {make_flag(cap_fill)}\\
      {make_flag(filter_intercept)}'))
    if (output_format == 'png') {
      o_fn <- paste0(of_base_name, '.png')
      png(o_fn,
          width = as.numeric(convertUnit(plot_width, 'inch')),
          height = as.numeric(convertUnit(plot_height, 'inch')),
          units = 'in', res = 400)
    } else if (output_format == 'pdf') {
      o_fn <- paste0(of_base_name, '.pdf')
      pdf(o_fn,
          width = as.numeric(convertUnit(plot_width, 'inch')),
          height = as.numeric(convertUnit(plot_height, 'inch')))
    } else {
      stop('Unrecognized output format')
    }
    tryCatch(eval(plot_expr), error = function(e) {
      print(e); success <<- F })
    dev.off()
    message(glue::glue('Wrote results to {o_fn}'))
    return(o_fn)
  }
}
## }}}
