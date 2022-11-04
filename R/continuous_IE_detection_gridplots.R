#!/usr/bin/env Rscript

# devtools::load_all(file.path('~/libs', 'maartenutils'))
## Combine stats objects into regression coefficient grid plots
# devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
# source(file.path('~/antigenic_space/libs/fasanalysis/R/plot_rates.R'))
ma_dir <- '~/antigenic_space/maarten-analyses'
source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_init.R'))
source('~/libs/maartenutils/R/system.R')
source('~/libs/fasanalysis/R/plot_rates.R')
# source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_helpers.R'))
# head(plot_IE_test_heatmap)

ncores <- 1
ncores <- 8
ncores <- 36
doParallel::registerDoParallel(cores = ncores)
per_allele <- F
options(error = recover)
options(error = traceback)

# dtf <- compile_all_coef_overview(redo = T)

if (per_allele) {
  dplyr::select(all_cont_IE_settings, -focus_allele) %>%
  unique %>%
  plyr::a_ply(1, function(x) {
    analysis_name <- as.character(unlist(x['analysis_name']))
    # if (analysis_name != 'twoD_sens_analysis') return(NULL)
    overlap_var <- as.character(unlist(x['overlap_var']))
    LOH_HLA <- as.logical(x['LOH_HLA'])
    patient_inclusion_crit <-
      as.character(unlist(x['patient_inclusion_crit']))

    if (F) {
      for (focus_allele in focus_hlas) {
        ## Determine proj order based on A0201
        if (focus_allele == 'A0201') {
          proj_order <- NULL
        }
        param_grid <- gen_continuous_param_grid(
          focus_allele = focus_allele,
          method = method,
          LOH_HLA = LOH_HLA,
          patient_inclusion_crit = patient_inclusion_crit,
          stats_idx = 1,
          ncores = 1,
          analysis_name = analysis_name,
          proj_order = proj_order,
          overlap_var = overlap_var,
          plot_sig_var = plot_sig_var,
          project_ranking = 'sum_fill_var')

        if (focus_allele == 'A0201') {
          proj_order <- attr(param_grid, 'proj_order')
        }
      }
    } else {
      param_hist <- gen_continuous_param_hist(
        focus_allele = focus_allele,
        patient_inclusion_crit = patient_inclusion_crit,
        method = method,
        stats_idx = 1,
        ncores = 1,
        analysis_name = analysis_name,
        overlap_var = overlap_var,
        plot_sig_var = plot_sig_var,
        project_ranking = 'sum_fill_var')
    }
  }, .parallel = (T && ncores > 1))
} else if (F) {
  expand.grid(overlap_var = overlap_vars,
              pssi = 1:6,
              analysis_name = analysis_names,
              hla_sim_range = list(NULL, c(0, 1))) %>%
  plyr::a_ply(1, function(x) {
    gen_continuous_param_grid_plot_w_power(
      hla_sim_range = unlist(x[['hla_sim_range']]),
      pssi = as.integer(unlist(x[['pssi']])),
      overlap_var = as.character(unlist(x[['overlap_var']])),
      ncores = ncores,
      plot_ranges = NULL,
      legend_pos = 'bottom',
      redo = F,
      return_res = T,
      analysis_name = as.character(unlist(x[['analysis_name']])))
  })
} else {
  IE_tumor_types <- c("Pan*'-'*cancer", 
    donor_summary[, levels(project_extended)])
  IE_tumor_types %>%
    purrr::map(complex_continuous_param_grid_plot_w_power, redo = F)
}
