library(RColorBrewer)
library(magrittr)
library(grDevices)
library(purrr)
IE_r_dir <- file.path(IE_root, 'R')

## Only load these packages when needed, loading can be slow
if (!any(grepl('maartenutils', search()))) {
  devtools::load_all(file.path('~/libs', 'maartenutils'))
}
if (!any(grepl('fasanalysis', search())) || T) {
  devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
}
# devtools::install_github('jokergoo/ComplexHeatmap')

source('~/libs/result_cacher.R')
source(file.path(IE_r_dir, 'continuous_IE_detection_helpers.R'))
source(file.path(IE_r_dir, 'model_fit.R'))
source(file.path(IE_r_dir, 'plotting.R'))
source(file.path(IE_r_dir, 'cached_functions.R'))
source(file.path(IE_r_dir, 'continuous_IE_checks.R'))
source(file.path(IE_r_dir, 'HLA_presentation_scores_helpers.R'))
source(file.path(IE_r_dir, 'pan_IE_settings_heatmap.R'))
source(file.path(IE_r_dir, 'IE_analysis.R'))


# all_cont_IE_settings <- tidyr::expand_grid(
#   focus_allele = focus_hlas,
#   overlap_var = c('mean_score', 'mean_score_AB'),
#   patient_inclusion_crit = list(c('TR'), c('')),
#   LOH_HLA = c('no_LOHHLA', 'LOHHLA', 'strict_LOHHLA'),
#   analysis_name = analysis_names
# )

all_cont_IE_settings <- tidyr::expand_grid(
  focus_allele = focus_hlas,
  overlap_var = c('mean_score', 'mean_score_AB'),
  # patient_inclusion_crit = c('strict_TR', 'TR', ''),
  patient_CYT = c('all', '>75', '<=75'),
  patient_inclusion_crit = c('none', 'FDR0.01', 'FDR1', 'FDR10'),
  LOH_HLA = c('no_LOHHLA', 'LOHHLA', 'strict_LOHHLA'),
  analysis_name = c('SNV_param_titration', 'clon_param_titration',
    'driv_ess_param_titration', 'marty_param_titration',
    'indel_param_titration')
)

if (exists('setMKLthreads')) setMKLthreads(1)
if (exists('setDTthreads')) {
  setDTthreads(1)
  stopifnot(getDTthreads() == 1)
}
stopifnot(test_dplyr())

datef <- format(Sys.time(), '%Y%m%d')

if (F) {
  datef <- format(Sys.time(), '%Y%m%d')
  pacman::p_load('lme4')
  method <- 'rlm_one_group'
  plot_sig_var = 'signif_effect.adj'
  mod_time <- '2019-08-08 08:27 CET'
  mod_time <- '2019-08-20 08:27 CET'
  mod_time <- '2020-01-09 12:27 CET'

  requires_computation <- maartenutils::gen_time_comparator(
    minimum_mod_time = mod_time,
    verbose = F)

  patient_selection_settings <-
    all_cont_IE_settings[, c('LOH_HLA', 'patient_inclusion_crit')] %>%
    unique %>%
    as.data.table

  basic_settings <-
    all_cont_IE_settings[, c('analysis_name', 'focus_allele', 'overlap_var')] %>%
    unique %>%
    as.data.table
}
