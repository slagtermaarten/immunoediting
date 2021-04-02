#!/usr/bin/env Rscript

#' SETTINGS
# ro_fn <- file.path(rds_dir, glue::glue('focus_allele_avg_HLA_presentation-hla_sim\\
#                                  _method_-cutoff-integration_method_\\
#                                  -union.rds'))
# RO_file <- file.path(rds_dir,
#                      glue::glue('focus_allele_avg_HLA_presentation-hla_sim_method\\
#                           _-cutoff-integration_method_-union.rds'))
fn_suffix = ''
fn_suffix = '-im_union'
do_test_grobs <- T
do_test_grobs <- F
min_cell_size <- .2
options(error = traceback)
## Used to be T
filter_problematic_rows <- F
## Used to be T
cluster_vals <- F
## Used to be T
restrict_to_powered_analyses <- F

## {{{ Init
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
# pacman::p_load(ComplexHeatmap)
if (!'fasanalysis' %in% loadedNamespaces())
  devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
devtools::load_all('~/libs/ComplexHeatmap')
pacman::p_load(circlize)
# devtools::install('~/libs/maartenutils')
pacman::p_load(maartenutils)
pacman::p_load(data.table)
pacman::p_load(RColorBrewer)
# ma_dir <- file.path('~/antigenic_space', 'maarten-analyses')
# rds_dir <- file.path(ma_dir, 'rds')
## The relative cell size assigned to the analyses with the lowest power
source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_init.R'))
# source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_helpers.R'))
# options(error = recover)
## }}} Init

## {{{ Data
analysis_names <- c('twoD_sens_analysis', 'rooney_param_titration',
                    'clon_param_titration', 'driv_ess_param_titration',
                    'marty_param_titration')


ds_param_grid <- structure(
  list(
    expression_threshold = c(NA, 0, 50, 500, 5000, NA, NA, NA, NA, NA, NA, 0,
      50, 500, 5000, NA, NA, NA, NA, NA, NA, 0, 50, 500, 5000, NA, NA, NA, NA,
      NA, NA, 0, 50, 500, 5000, NA, NA, NA, NA, NA, NA, 0, 50, 500, 5000, NA,
      NA, NA, NA, NA, NA, 0, 50, 500, 5000, NA, NA, NA, NA, NA, NA, 0, 50, 500,
      5000, NA, NA, NA, NA, NA, NA, 0, 50, 500, 5000, NA, NA, NA, NA, NA),
    VE_threshold = c(NA, NA, NA, NA, NA, 0, 3, 10, 50, 100, NA, NA, NA, NA, NA,
      0, 3, 10, 50, 100, NA, NA, NA, NA, NA, 0, 3, 10, 50, 100, NA, NA, NA, NA,
      NA, 0, 3, 10, 50, 100, NA, NA, NA, NA, NA, 0, 3, 10, 50, 100, NA, NA, NA,
      NA, NA, 0, 3, 10, 50, 100, NA, NA, NA, NA, NA, 0, 3, 10, 50, 100, NA, NA,
      NA, NA, NA, 0, 3, 10, 50, 100),
    sts_filtering = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE),
    percentile_rank = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9,
      1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4)),
  class = c("data.table", "data.frame"),
  row.names = c(NA, -80L)
)
# }}}

## {{{ Helper functions
trans_power_size <- function(x) {
  min_cell_size + x * (1 - min_cell_size)
}

ceiling_decimal <- function(x, dec) {
  ceiling(x * 10^dec) / 10^dec
}
round_up <- ceiling_decimal


floor_decimal <- function(x, dec) {
  floor(x * 10^dec) / 10^dec
}
round_down <- floor_decimal


matrix_to_vec <- function(X, ei) {
  X <- n_power_mat
  r_i <- ei %% dim(X)[1]
  c_i <- floor(ei / dim(X)[2])
  # X_vec <- as.vector(X)
}


reset_plotting_device <- function() {
  while (length(dev.list()) > 0) {
    dev.off()
  }
}
# reset_plotting_device()


test_grob <- function(grob, doit = F) {
  if (!doit) return(invisible())
  reset_plotting_device()
  dev.new()
  pushViewport(viewport(width = 0.9, height = 0.9))
  grid.rect()
  draw(grob)
  popViewport()
}


test_legend <- function(grob, doit = F){
  if (!doit) return(invisible())
  reset_plotting_device()
  dev.new(width = 3, height = 3)
  pushViewport(viewport(width = 0.9, height = 0.9))
  grid.rect()  # border
  draw(grob,
    x = unit(.5, "npc"),
    y = unit(.5, "npc"),
    just = c("center", "center")
  )
  popViewport()
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
# identify_problematic_combs(rc_mat)
# remove_problematic_combs(rc_mat)

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
row_cols <- gen_color_vector(name = 'Royal1', arg = 9)
analysis_name_cols <- c(gen_color_vector(name = 'Royal1', arg = 4),
  'darkseagreen4')
column_annotation_colors = list(
  et = gen_cont_colors(c(0, 50, 5000)),
  STS = setNames('darkorange3', 'STS'),
  PR = gen_cont_colors(c(perc_ranks, NA),
    c('darkslategray4', 'white'))[seq_along(perc_ranks)],
  VE = gen_cont_colors(c(0, 10, 100), c('white', 'goldenrod3'))
)
## }}}

## {{{ Constants/settings
ID_cols <- c('overlap_var', 'patient_inclusion_crit', 'LOH_HLA',
             'hla_sim_range', 'analysis_name')
# unique(dtf, by = ID_cols)
model_formula <- analysis_name + overlap_var + patient_inclusion_crit +
  LOH_HLA + hla_sim_range ~ analysis_idx

point_size <- unit(.5, 'mm')
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
  analysis_name = 'Variant selection',
  overlap_var = 'Computation of\npresentation score',
  patient_inclusion_crit = 'Removing IT\nresistant patients',
  hla_sim_range = 'Presentation score\npatient filtering',
  LOH_HLA = 'Adapting presentation score\nto LOH in HLA')
ra_titles_short = c(
  analysis_name = 'Variant selection',
  overlap_var = 'Presentation score',
  patient_inclusion_crit = 'IT-resistant patients',
  hla_sim_range = 'PS filtering',
  LOH_HLA = 'LOH in HLA')
if (fn_suffix == '-im_union') {
  ra_titles_short <- ra_titles_short[names(ra_titles_short) != 'hla_sim_range']
}
perc_ranks <- c(1, 1.9, 3, 4)
output_dir <- '~/antigenic_space/maarten-analyses/img/cont_IE_heatmaps_V2'
# output_dir <- '~/Dropbox/cont_IE_heatmaps_V2'
## }}}

## {{{ Main function
plot_pan_IE_heatmap <- function(tumor_type = 'pan_cancer',
                                restrict_to_powered_analyses = T,
                                fill_var = 'rc',
                                return_res = F,
                                cell_size = unit(1.5, 'mm'),
                                adaptive_p_val_bound = T,
                                PPV = .45,
                                ann_size = 1.5 * cell_size,
                                ncores = 1,
                                redo = F,
                                subset_perc = 1,
                                p_val_bound = NULL,
                                debug_func = F,
                                compress_rc = NULL,
                                do_test_grobs = F) {
  crc_str <- format_flag(val = compress_rc, name = 'compress_rc')
  fvf <- ifelse(fill_var == 'rc', '', '-CYT_rc')
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = T)
  tumor_type_simple <- tumor_types[tumor_type]

  if (debug_func) {
     browser()
  }
  # messagef('Starting %s', tumor_type)
  ## prep_pan_IE_heatmap --> compile_all_coef_overview --> summarize HLA alleles
  dtf <- prep_pan_IE_heatmap(p_val_bound = p_val_bound,
                             tumor_type = tumor_type,
                             redo = redo,
                             fn_suffix = fn_suffix,
                             PPV = PPV,
                             fill_var = fill_var,
                             ncores = ncores,
                             adaptive_p_val_bound = adaptive_p_val_bound)

  if (fn_suffix == '-im_union' && 'hla_sim_range' %in% colnames(dtf)) {
    dtf <- dtf[hla_sim_range == 'all']
  }

  # if (F) {
  #   hist(dtf$es)
  #   dput(as.list(dtf[eps(es, 0.05)][1]))
  #   dtf[eps(es, 4)]
  #   prob_args <- as.list(dtf[eps(es, 4)][1])
  #   library(magrittr)
  #   prob_args %<>%
  #     { .[intersect(names(.), names(formals(optim_compute_power_continuous_IE)))] }
  #   prob_args %>% { do.call(optim_compute_power_continuous_IE, .) }
  # }

  if (null_dat(dtf)) return(NULL)
  stopifnot(is.numeric(subset_perc) && subset_perc <= 1 && subset_perc > 0)
  # browser()
  # plot(dtf[!is.na(signif_effect.adj), .(p_val, p_val.adj)])
  # plot(dtf[signif_effect.adj == 'L', .(p_val, p_val.adj)])
  # plot(dtf[signif_effect.adj == 'G', .(p_val, p_val.adj)])
  if (F) {
    if (!is.null(p_val_bound)) {
      dtf[p_val.adj > p_val_bound, signif_effect.adj := NA]
    }
  }

  RC_mat <- tryCatch(dcast(dtf, model_formula, value.var = fill_var),
                     error = function(e) { print(e); browser() })
  RC_mat <- rbind(
    cbind('regression_var' = 'yield_rate',
          copy(RC_mat[analysis_name != 'rooney_param_titration'])),
    cbind('regression_var' = 'rooney',
          copy(RC_mat[analysis_name == 'rooney_param_titration'])) %>%
      .[, analysis_name := 'twoD_sens_analysis'])
        
  RC_mat_IDs <- RC_mat[, 1:6]
  # rc_mat <- log10(rc_mat + 1)
  rc_mat <- as.matrix(RC_mat[, 7:ncol(RC_mat)]) %>%
    # { sign(.) * (abs(.)^1/10) }
    identity

  p_var <- switch(fill_var, 'rc' = 'p_val', 'rc_CYT' = 'p_val_CYT')
  p_var.adj <- sprintf('%s.adj', p_var)
  signif_effect.adj <- sprintf('%s_signif_effect.adj', p_var)
  # es_var <- switch(fill_var, 'rc' = 'p_val', 'rc_CYT' = 'p_val_CYT')

  icon_mat <- dcast(dtf, model_formula, value.var = signif_effect.adj) %>%
    as.matrix %>%
    .[, 6:ncol(.)]
  icon_freqs <- table(icon_mat[!is.na(icon_mat)])
  icon_freqs <- setNames(c(icon_freqs['L'], icon_freqs['G']), c('L', 'G'))
  icon_freqs <- repl.na(icon_freqs)

  # icon_mat[which(!is.na(icon_mat), arr.ind = T)]
  # mean(is.na(as.vector(icon_mat)))

  if ('es' %in% colnames(dtf)) {
    power_mat <- dcast(dtf, model_formula, value.var = 'es') %>%
      as.matrix %>%
      .[, 6:ncol(.)] %>%
      apply(2, as.numeric)
    messagef('%s of cells annotated with power analysis data',
      scales::percent(1 - mean(is.na(as.vector(power_mat)))))
    ## power_mat contains 'es' values, lower is better
  } else {
    message('No cells annotated with power analysis data')
    power_mat <- matrix(rep(NA, prod(dim(icon_mat))), ncol = ncol(icon_mat))
  }

  n_power_mat <- { power_mat + 1 } %>%
    log10 %>%
    {
      min_p <- min(., na.rm = T)
      max_p <- max(., na.rm = T)
      1 - ((. - min_p) / (max_p - min_p))
    } %>%
    trans_power_size

  ## NAs are okay, Infs, -Inf, and NaN not
  if (any(!is.finite(power_mat) & !is.na(power_mat))) {
    idxs <- which(!is.finite(power_mat) & !is.na(power_mat), arr.ind = T)
    power_mat[idxs[1, ]]
    perm_browser()
  }

  if (restrict_to_powered_analyses) {
    rc_mat[is.na(power_mat)] <- NA
    icon_mat[is.na(power_mat)] <- NA
    n_power_mat[is.na(power_mat)] <- NA
  }

  informative_cell_idx <-
    which(!is.na(rc_mat) & !eps(rc_mat, 0, 1e-12), arr.ind = T)
  if (fill_var == 'rc' && length(informative_cell_idx) == 0) {
    return(NULL)
  }

  ## Subselect informative (i.e. not all NA) rows and columns
  if (filter_problematic_rows) {
    problematic_combs <- identify_problematic_combs(rc_mat)
    i_rows <- setdiff(1:nrow(rc_mat), problematic_combs$row)
    i_rows_alt <- tryCatch(apply(rc_mat, 1, function(x) !all(is.na(x))),
      error = function(e) { print('Could not subset rows'); 1:nrow(rc_mat) })
    ## Subselect a percentage of rows
    i_rows <- i_rows[1:ceiling(length(i_rows) * subset_perc)]
    i_columns <- setdiff(1:ncol(rc_mat), problematic_combs$column)
    # i_columns <- tryCatch(apply(rc_mat, 2, function(x) !all(is.na(x))),
    #   error = function(e) { print('Could not subset columns'); 1:col(rc_mat) })
    ## Subselect a percentage of columns
    i_columns <- i_columns[1:ceiling(length(i_columns) * subset_perc)]
  } else {
    i_columns <- 1:ncol(rc_mat)
    i_rows <- 1:nrow(rc_mat)
  }
  nc <- length(i_columns)
  nr <- length(i_rows)
  rc_mat <- rc_mat[i_rows, i_columns]
  power_mat <- power_mat[i_rows, i_columns]
  icon_mat <- icon_mat[i_rows, i_columns]
  RC_mat_IDs <- RC_mat_IDs[i_rows, ]
  # col_pal <- maartenutils::darken(fas_discrete_colors(3), 1.1) %>% { .[c(1, 3)] }

  x_max <- max(abs(rc_mat), na.rm = T)
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
      segments(c(-15, x_mid), c(.5, -5), c(x_mid, x_mid), c(.5, .5), col = 'red')
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
    # seq(-3, 3, by = .25) %>% { eps(., trans_from(trans_to(.))) }
    trans_to <- function(x) x
    x_mid <- 0
  }

  # with(list(x = seq(-1, 1, by = .1)), plot(x, trans_from(x), type = 'p'))
  rc_res <- ceiling(abs(log10(x_max)))
  x_max_norm <- round_up(x_max, dec = rc_res)
  # rc_breaks <- seq(trans_to(-x_max_norm), trans_to(x_max_norm), length.out = 5)
  # rc_labels <- rc_breaks %>% trans_from
  ## We want between 5 and 9 markers, more markers if compression is high
  rc_labels <- seq(-x_max_norm, x_max_norm,
                   length.out = min(max(5, compress_rc), 9)) %>%
    round(rc_res + 1)
  rc_breaks <- trans_to(rc_labels)
  LOHHLA_vec <- RC_mat_IDs$LOH_HLA %>% { ifelse(grepl('no_', .), NA, .) }
  LOHHLA_levels <- setdiff(unique(LOHHLA_vec), NA)

  # rescale <- function(x) (x - min(abs(x), na.rm = T)) /
  #   (max(abs(x), na.rm = T) - min(abs(x), na.rm = T))
  # seq(-3, 3, by = .25) %>% rescale
  # seq(-3, 3, by = .25) %>% { trans_to(.) } %>% rescale

  ## Heat map specific color settings
  rc_col_fun <- trans_to(x_max) %>%
    setNames(NULL) %>%
    { c(-1 * ., 0, .) } %>%
    { colorRamp2(colors = c('#A30D1D', my_grey, '#3F9CB5'), breaks = .) }
  LOHHLA_cols <- darken(rep(row_cols[7], 2), c(1, 1.5)) %>%
    setNames(c('LOHHLA', 'strict_LOHHLA'))
  LOHHLA_cols <- LOHHLA_cols[LOHHLA_levels]
  row_annotation_colors <- list(
    # regression_var = setNames(analysis_name_cols, analysis_names),
    analysis_name = setNames(analysis_name_cols, analysis_names),
    LOH_HLA = LOHHLA_cols,
    patient_inclusion_crit = setNames(darken(row_cols[8], 1), 'TR'),
    overlap_var = setNames(darken(row_cols[9], 1), c('mean_score_AB')),
    hla_sim_range = setNames(darken(row_cols[1], 1), c('<= 1'))
  )

  row_annotation_settings <- list(
    regression_var = list(
      title = ra_titles['regression_var'],
      breaks = c('yield_rate', 'rooney'),
      labels = c('Neo-antigen yield rate', 'Observed / expected neo-antigens'),
      values = RC_mat_IDs$regression_var),
    analysis_name = list(
      title = ra_titles['analysis_name'],
      breaks = analysis_names[c(1, 3, 4, 5)],
      labels = c('All variants', 'Clonal',
        'No drivers/essentials', 'Marty oncogenic'),
      values = RC_mat_IDs$analysis_name),
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
      values = RC_mat_IDs$patient_inclusion_crit %>%
      { ifelse(. == 'TR', ., NA) }),
    overlap_var = list(
      title = ra_titles['overlap_var'],
      breaks = 'mean_score_AB',
      labels = 'Restricted to A & B alleles',
      values = RC_mat_IDs$overlap_var %>%
      { ifelse(. == 'mean_score_AB', ., NA) }),
    hla_sim_range = list(
      title = ra_titles['hla_sim_range'],
      breaks = '<= 1',
      labels = c('No pts \\w PS > 1'),
      values = RC_mat_IDs$hla_sim_range %>%
      { ifelse(. == '<= 1', ., NA) })
  )

  lgds_to_make <- setdiff(names(row_annotation_settings), 'regression_var')
  if (fn_suffix == '-im_union') {
    lgds_to_make <- setdiff(lgds_to_make, 'hla_sim_range')
  }

  row_annotation_legends <-
    lapply(lgds_to_make,
      function(an) {
        at_vals <- unique(row_annotation_settings[[an]][['values']])
        if (is.logical(at_vals)) at_vals <- setdiff(at_vals, FALSE)
        at_vals <- setdiff(at_vals, NA)
        label_vals <- at_vals
        if (!is.null(row_annotation_settings[[an]][['labels']]))
          label_vals <- row_annotation_settings[[an]][['labels']]
        cols <- row_annotation_colors[[an]][at_vals]
        if (length(at_vals) != length(cols)) {
          print('Not enough colors')
          browser()
        }
        lgd <- Legend(
          at = at_vals, labels = label_vals, legend_gp = gpar(fill = cols),
          title = row_annotation_settings[[an]][['title']],
          nrow = row_annotation_settings[[an]][['nrow']],
          direction = row_annotation_settings[[an]][['legend_direction']],
          grid_width = 2 * cell_size,
          grid_height = 2 * cell_size,
          title_gp = default_legend_settings$title_gp,
          labels_gp = default_legend_settings$labels_gp)
        # browser()
        # class(lgd)
        # test_legend(lgd, do_test_grobs)
        return(lgd)
      })
  # test_grob(row_annotation_legends[[1]], do_test_grobs)
  # test_grob(row_annotation_legends[[2]], do_test_grobs)
  # test_grob(row_annotation_legends[[3]], do_test_grobs)
  # test_grob(row_annotation_legends[[4]], do_test_grobs)
  # test_grob(row_annotation_legends[[5]], do_test_grobs)


  row_annotation_args <- list(
    # regression_var = row_annotation_settings$regression_var$values,
    analysis_name = row_annotation_settings$analysis_name$values,
    LOH_HLA = row_annotation_settings$LOH_HLA$values,
    patient_inclusion_crit = row_annotation_settings$patient_inclusion_crit$values,
    overlap_var = row_annotation_settings$overlap_var$values,
    which = 'row',
    na_col = 'white',
    gp = gpar(width = ann_size, lwd = 0, col = 'white'),
    col = row_annotation_colors,
    # simple_anno_size = ann_size,
    annotation_width = c(40, 1, 1, 1, 1, 1) * ann_size,
    width = ann_size,
    gap = unit(c(0, 0, 0, 0, 0), 'mm'),
    show_annotation_name = F)
  if (fn_suffix != '-im_union') {
    row_annotation_args[['hla_sim_range']] <-
      row_annotation_settings$hla_sim_range$values
  }
  row_annotation <- do.call(HeatmapAnnotation, row_annotation_args)
  test_grob(row_annotation, do_test_grobs)

  column_annotation_legend_settings <-
    with(as.list(ds_param_grid[as.integer(colnames(rc_mat))[i_columns], ]), {
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
          title = 'Variant coverage threshold\n(# reads)',
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
          grid_width = 2 * cell_size,
          grid_height = 2 * cell_size,
          title_gp = default_legend_settings$title_gp,
          labels_gp = default_legend_settings$labels_gp)
        test_legend(lgd, do_test_grobs)
        return(lgd)
      })
  test_grob(column_annotation_legends[[1]], do_test_grobs)
  test_grob(column_annotation_legends[[2]], do_test_grobs)
  test_grob(column_annotation_legends[[3]], do_test_grobs)
  test_grob(column_annotation_legends[[4]], do_test_grobs)

  column_annotation <-
    with(as.list(ds_param_grid[as.integer(colnames(rc_mat))[i_columns], ]),
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

  yield_rate_ri <- RC_mat_IDs[, which(regression_var != 'rooney')]
  rooney_ri <- RC_mat_IDs[, which(regression_var == 'rooney')]

  yield_rate_HM <- ComplexHeatmap::Heatmap(
    rc_mat[yield_rate_ri, ],
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
      r_power <- n_power_mat[yield_rate_ri, ][i, j]
      if (!is.na(r_power) && !is.na(rc_mat[i, j])) {
        ## Sqrt transform the normalized power in order to get the required side
        ## length (of a square)
        w <- h <- cell_size * sqrt(r_power)
        grid.rect(x = x, y = y,
                  width = w, height = h,
                  gp = gpar(lwd = 0,
                            fill = rc_col_fun(trans_to(rc_mat[yield_rate_ri, ][i, j])),
                            col = 'white'))
      }
      if (!is.na(icon_mat[yield_rate_ri, ][i, j]) && icon_mat[yield_rate_ri, ][i, j] != '') {
        if (icon_mat[yield_rate_ri, ][i, j] == 'L') {
          l_pch = 6
        } else if (icon_mat[yield_rate_ri, ][i, j] == 'G') {
          l_pch = 2
        }
        grid.points(x, y, pch = l_pch, size = point_size)
      }
    },
    rect_gp = gpar(type = 'none'),
    bottom_annotation = NULL,
    show_column_dend = F,
    show_row_dend = F)

  rooney_HM <- ComplexHeatmap::Heatmap(
    rc_mat[rooney_ri, ],
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
      r_power <- n_power_mat[rooney_ri, ][i, j]
      if (!is.na(r_power) && !is.na(rc_mat[i, j])) {
        ## Sqrt transform the normalized power in order to get the required side
        ## length (of a square)
        w <- h <- cell_size * sqrt(r_power)
        grid.rect(x = x, y = y,
                  width = w, height = h,
                  gp = gpar(lwd = 0,
                            fill = rc_col_fun(trans_to(rc_mat[rooney_ri, ][i, j])),
                            col = 'white'))
      }
      if (!is.na(icon_mat[rooney_ri, ][i, j]) && icon_mat[rooney_ri, ][i, j] != '') {
        if (icon_mat[rooney_ri, ][i, j] == 'L') {
          l_pch = 6
        } else if (icon_mat[rooney_ri, ][i, j] == 'G') {
          l_pch = 2
        }
        grid.points(x, y, pch = l_pch, size = point_size)
      }
    },
    rect_gp = gpar(type = 'none'),
    bottom_annotation = column_annotation,
    show_column_dend = F,
    show_row_dend = F)
  # test_grob(rooney_HM, T)
  # test_grob(rooney_HM, do_test_grobs)
  if (do_test_grobs) print(rooney_HM)

  ## Legend
  # devtools::load_all('~/libs/ComplexHeatmap')
  power_legend_breaks <- range(power_mat, na.rm = T) %>%
    { . + 1 } %>%
    log10 %>%
    { seq(.[1], .[2], length.out = 5) } %>%
    { 10^. - 1 } %>%
    { -1 * . }
  sizes <- seq(1, 0, length.out = 5) %>% trans_power_size
  # browser()
  power_lgd <- tryCatch(Legend(
    at = power_legend_breaks,
    legend_gp = c(list(fill = 'grey80')),
    grid_width = cell_size,
    grid_height = cell_size,
    row_gap = .7 * cell_size,
    title_gp = default_legend_settings$title_gp,
    labels_gp = default_legend_settings$labels_gp,
    cell_scaling_factors = sizes,
    labels = c(0,
      # fancy_scientific(power_legend_breaks[2:length(power_legend_breaks)],
      #   digits = 2)),
      round(power_legend_breaks[2:length(power_legend_breaks)], 1)),
    title = 'IE detection power\n(required slope for significance)'
  ), error = function(e) { print(e); browser() })
  test_legend(power_lgd, do_test_grobs)

  # do_test_grobs = T
  icon_labels <- mapply(sprintf, rep('"%s"~(italic(n)==%d)', 2),
    c('Negative', 'Positive'), icon_freqs[c('L', 'G')]) %>%
    setNames(NULL) %>%
    sapply(function(x) parse(text = x))
  icon_lgd <- tryCatch(Legend(
    at = c('L', 'G'),
    title = 'Significant\nregression slope\n(FDR-corrected)',
    grid_height = 2 * cell_size,
    grid_width = 2 * cell_size,
    labels = icon_labels,
    # gap = unit(40, 'mm'),
    size = point_size,
    legend_gp = c(default_legend_settings,
      list(fill = 'white', fontsize = .5,
           width = ann_size, lwd = 1, col = 'black')),
    type = 'points',
    title_gp = default_legend_settings$title_gp,
    labels_gp = default_legend_settings$labels_gp,
    pch = rev(c(2, 6))
  ), error = function(e) { print(e); browser() })
  test_legend(icon_lgd, do_test_grobs)
  

  # rc_lgd <- tryCatch(Legend(
  #   at = rc_breaks,
  #   legend_gp = c(list(fill = rc_col_fun(rc_breaks))),
  #   grid_height = 2 * cell_size,
  #   grid_width = 2 * cell_size,
  #   title_gp = default_legend_settings$title_gp,
  #   labels_gp = default_legend_settings$labels_gp,
  #   labels = rc_labels,
  #   title = 'Regression coefficient'
  # ), error = function(e) { print(e); browser() })
  # test_legend(rc_lgd, do_test_grobs)

  rc_lgd <- tryCatch(Legend(
      at = rc_breaks,
      col_fun = rc_col_fun,
      # legend_gp = c(list(fill = rc_col_fun(rc_breaks))),
      grid_height = 2 * cell_size,
      grid_width = 2 * cell_size,
      title_gp = default_legend_settings$title_gp,
      labels_gp = default_legend_settings$labels_gp,
      labels = rc_labels,
      title = switch(fill_var,
        'rc' = 'Regression coefficient',
        'rc_CYT' = 'Regression coefficient (CYT)')
  ), error = function(e) { print(e); browser() })
  test_legend(rc_lgd, do_test_grobs)

  legend_header <- function(legend_title = 'test') {
    Legend(title = legend_title, at = c(0), labels = '',
      gap = unit(0, 'mm'),
      title_gp = gpar(fontsize = base_size, fontface = 'bold'),
      grid_height = unit(0, 'mm'))
  }
  test_legend(legend_header(), do_test_grobs)

  legends_grob <- do.call(packLegend, c(
      list(legend_header(switch(fill_var,
                                'rc' = 'Immunoediting regression',
                                'rc_CYT' = 'Immunoediting regression (CYT)')),
           rc_lgd, icon_lgd, power_lgd),
      list(legend_header('Immunoediting\nanalysis settings')),
      row_annotation_legends,
      list(legend_header('Neo-antigen prediction\npipeline settings')),
      column_annotation_legends,
      list(direction = 'vertical')))
  test_grob(legends_grob, do_test_grobs)

  comb_yield_rate <- row_annotation[yield_rate_ri, ] + yield_rate_HM
  comb_rooney <- row_annotation[rooney_ri, ] + rooney_HM
  # comb <- main_heatmap
  # heatmap_width <- ComplexHeatmap:::width(draw(comb))
  # heatmap_width <- ComplexHeatmap:::width(column_annotation)
  legend_width <- ComplexHeatmap:::width(legends_grob) + unit(2, 'mm')
  legend_height <- ComplexHeatmap:::height(legends_grob)

  ## US letter size: 216 x 279 mm
  plot_height <- { cell_size * (nr / 2) + unit(4.2, 'cm') } %>%
    max(legend_height + unit(5, 'cm')) %>%
    min(unit(25, 'cm'))
  ## Plot width cannot be more than 17.4
  plot_width <- { cell_size * nc + unit(10, 'cm') } %>%
    min(unit(17.4, 'cm'))
  messagef('Plot dimensions %.1f [cm] x %.1f [cm]',
    as.numeric(convertUnit(plot_width, 'cm')),
    as.numeric(convertUnit(plot_height, 'cm')))

  plot_expr <- rlang::expr({
    left_width <- unit(1, 'npc') - legend_width
    left_vp <- viewport(
        just = c('left', 'bottom'),
        width = left_width,
        height = 1,
        x = 0, y = unit(0, 'npc'))
    right_vp <- viewport(
        just = c('left', 'bottom'),
        width = legend_width,
        height = 1,
        x = left_width, y = unit(0, 'npc'))
    column_anno_h <- unit(1, 'cm') + ComplexHeatmap::height(column_annotation)
    yield_rate_HM_h <- unit(length(yield_rate_ri) / nrow(rc_mat), 'npc') - 
      column_anno_h
    rooney_HM_h <- unit(length(rooney_ri) / nrow(rc_mat), 'npc') + 
      column_anno_h
    yield_rate_vp <- viewport(
        just = c('left', 'bottom'),
        width = 1,
        height = yield_rate_HM_h, 
        # just = rev(c('top', 'right')),
        # just = rev(c('bottom', 'left')),
        # just = rev(c('bottom', 'left')),
        x = 0,
        y = rooney_HM_h)
    rooney_vp <- viewport(
        just = c('left', 'bottom'),
        width = 1,
        height = rooney_HM_h,
        # just = rev(c('top', 'right')),
        # just = rev(c('bottom', 'left')),
        # just = rev(c('bottom', 'left')),
        x = 0,
        y = unit(0, 'npc'))
    pushViewport(left_vp)
    pushViewport(yield_rate_vp)
    tryCatch(
      draw(comb_yield_rate,
           # padding = unit(c(10, 5, 5, 5), 'mm'), ## bottom, left, top, right
           padding = unit(c(10, 10, 10, 10), 'mm'), ## bottom, left, top, right
           newpage = F,
           show_annotation_legend = F,
           show_heatmap_legend = F,
           merge_legend = T),
      error = function(e) { print(e) })
    popViewport()
    pushViewport(rooney_vp)
    tryCatch(
      draw(comb_yield_rate,
           # padding = unit(c(10, 5, 5, 5), 'mm'), ## bottom, left, top, right
           padding = unit(c(10, 10, 10, 10), 'mm'), ## bottom, left, top, right
           newpage = F,
           show_annotation_legend = F,
           show_heatmap_legend = F,
           merge_legend = T),
      error = function(e) { print(e) })
    for (an in names(ca_titles)) {
      decorate_annotation(an, {
        grid.text(ca_titles[an],
          x = unit(1, 'npc') + unit(2, 'mm'),
          just = 'left',
          gp = do.call(gpar,
            c(annotation_par, list(fontface = 'italic'))))
      })
    }
    for (an in names(ra_titles_short)) {
      decorate_annotation(an, {
        grid.text(ra_titles_short[an],
          y = unit(-2, 'mm'), just = 'right',
          gp = do.call(gpar,
            c(annotation_par, list(fontface = 'italic'))),
          rot = 90)
      })
    }
    decorate_heatmap_body('rc', {
      grid.lines(c(0.0, 0.0, 1, 1, 0), c(0, 1, 1, 0, 0),
                 gp = border_line_style)
    })
    decorate_annotation('analysis_name', {
      grid.text(label = 'Immunoediting analysis settings',
        x = unit(0, 'npc') - unit(4, 'mm'),
        gp = do.call(gpar, c(axis_par, list())),
        y = unit(0.5, 'npc'),
        rot = 90,
        just = 'center')
    })
    popViewport()
    pushViewport(right_vp)
    # grid.text('Legends')
    draw(legends_grob,
      just = c('right', 'top'),
      # x = unit(1, 'npc') - unit(5, 'mm'),
      # y = unit(1, 'npc') - unit(5, 'mm'))
      x = unit(1, 'npc') - unit(10, 'mm'),
      y = unit(1, 'npc') - unit(10, 'mm'))
      # x = unit(5, 'mm'),
      # y = unit(5, 'mm'),
      # just = 'center')
    # grid.rect()
    popViewport()
    # decorate_annotation('PR', {
    #   grid.lines(unit(c(0, 0), 'mm'), unit(c(1, 1), 'npc'))
    # })
    # legend_x <- unit(.78, 'npc')
    # draw(rc_lgd,
    #   x = legend_x,
    #   y = unit(1, 'npc'),
    #   just = c('left', 'top')
    # )
    # draw(power_lgd,
    #   x = legend_x,
    #   y = unit(.93, 'npc'),
    #   just = c('left', 'top')
    # )
    # draw(icon_lgd,
    #   x = legend_x,
    #   y = unit(.85, 'npc'),
    #   just = c('left', 'top')
    # )
  })
  # eval(plot_expr)
  # panB <- grid.grabExpr(eval(plot_expr))
  # str(panB)[[1]]
  # panB <- grid.grabExpr(plot_expr)
  # success <- T
  # sys_file_open(o_fn)
  # if (!success) { file.remove(o_fn) } else { sys_file_open(o_fn) }

  if (return_res) {
    op <- grid.grabExpr(eval(plot_expr))
    # attr(op, 'dim') <- c(nr, nc)
    return(op)
  } else {
    o_fn <- glue::glue('{output_dir}/{tumor_type_simple}{crc_str}{fvf}.pdf')
    pdf(o_fn,
        width = convertUnit(plot_width, 'inch'),
        height = convertUnit(plot_height, 'inch'))
    eval(plot_expr)
    # tryCatch(eval(plot_expr), error = function(e) {
    #   print(e); success <<- F })
    dev.off()
    message(glue::glue('Wrote results to {o_fn}'))
    return(o_fn)
  }
}
## }}}

if (T && sys.nframe() %in% c(0, 4)) {
  expand.grid(
    # compress_rc = c(12, 8, 5, 2, 0),
    compress_rc = c(0, .5, 1),
    # compress_rc = c(12, 8),
    # compress_rc = c(5, 2, 0),
    # compress_rc = c(2),
    # compress_rc = c(0),
    # tumor_type = names(tumor_types)) %>%
    tumor_type = 'Pan*\'-\'*cancer') %>%
    # tumor_type = tumor_types[1]) %>%
    dplyr::mutate(idx = 1:dplyr::n()) %>%
    # { .[9:nrow(.), ] } %>%
    # tumor_types[1] %>%
    # 'melanoma' %>%
    # 'breast_basal' %>%
    plyr::a_ply(1, function(r) {
      idx <- as.integer(r[['idx']])
      messagef('Processing heatmap index %d', idx)
      tryCatch(plot_pan_IE_heatmap(
          tumor_type = as.character(r[['tumor_type']]),
          # restrict_to_powered_analyses = F,
          restrict_to_powered_analyses = restrict_to_powered_analyses,
          # debug_func = T,
          debug_func = F,
          fill_var = 'rc',
          # fill_var = 'rc_CYT',
          p_val_bound = .25,
          # ncores = 32,
          ncores = 40,
          # ncores = 1,
          adaptive_p_val_bound = T,
          do_test_grobs = do_test_grobs,
          return_res = F,
          compress_rc = as.numeric(r[['compress_rc']])),
      error = function(e) {
        print(sprintf('Failed for index %d', idx)); print(e)
      })
    })
  # purrr::map(tumor_types, plot_pan_IE_heatmap, restrict_to_powered_analyses = T)

  # y = NULL
  # for (sp in c(.5, 1)) {
  #   ht <- plot_pan_IE_heatmap(tumor_types[1], restrict_to_powered_analyses = T,
  #     subset_perc = sp, return_res = T, do_test_grobs = do_test_grobs)
  #
  #   op
  #   ht_height = sum(component_height(ht)) + unit(4, "mm")
  #   ht_height = convertHeight(ht_height, "inch", valueOnly = TRUE)
  #   y = c(y, ht_height)
  # }
}
