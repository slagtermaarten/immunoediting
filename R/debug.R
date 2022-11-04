# devtools::load_all(file.path('~/libs', 'maartenutils'))
# devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
# source('~/antigenic_space/bin/init.R')
source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_helpers.R'))
source(file.path(ma_dir, 'immune_editing', 'IE_analysis.R'))

base_settings <- list(
  z_normalize = F,
  include_non_ds_tallies = F,
  ncores = 1,
  filter_depletion_outliers = F,
  p_var = 'p_val'
)

thresholds <- list(.5, .25, .1, .05, NA)


cat('Expression prot\n')
list(affinity_expression_prot_idxs) %>% 
  walk(function(idxs) {
    analysis_settings <- list(
      hl_lgd_name = 'Gene expression threshold',
      x_var = 'expression_threshold',
      hl_name = 'expression_threshold',
      hl_func = function(x) {
        x %>%
          dplyr::mutate(hl_var = ifelse(is.na(expression_threshold), 'none',
                                        as.character(expression_threshold))) %>%
          dplyr::mutate(hl_var = factor(hl_var,
                                        levels = c('none', '0', 
                                                   '50', '5000'))) %>%
          { . }
      })
    map(thresholds, 
        ~exec(gen_depletion_tally_tables, 
              idxs = idxs,
              !!!base_settings,
              !!!analysis_settings,
              filter_group_cov = .x))
    map(thresholds, 
        ~exec(IE_lineplot, 
              idxs = idxs,
              !!!base_settings,
              !!!analysis_settings,
              filter_group_cov = .x))

    analysis_settings <- list(
      hl_lgd_name = 'Binding affinity rank percentile',
      x_var = 'percentile_rank',
      hl_name = 'affinity_threshold',
      hl_func = function(x) {
        x %>%
          dplyr::mutate(hl_var = ifelse(is.na(percentile_rank), 'none',
                                        as.character(percentile_rank))) %>%
          dplyr::mutate(hl_var = factor(hl_var,
                                        levels = c('none', '4', '3',
                                                   '1.9', '1'))) %>%
          { . }
      }
    )
    map(thresholds, 
        ~exec(gen_depletion_tally_tables, 
              idxs = idxs,
              !!!base_settings, 
              !!!analysis_settings, 
              filter_group_cov = .x))
    map(thresholds, 
        ~exec(IE_lineplot, 
              idxs = idxs,
              !!!base_settings,
              !!!analysis_settings,
              filter_group_cov = .x))

  })

cat('Variant expression prot\n')
list(
  # main_analysis_idxs,
  # main_plus_affinity_na_idxs,
  # affinity_VE_prot_idxs,
  # affinity_idxs,
  affinity_VE_prot_idxs
  ) %>% 
  walk(function(idxs) {
    analysis_settings <- list(
      hl_lgd_name = 'Variant expression threshold',
      x_var = 'VE_threshold',
      hl_name = 'variant_expression_threshold',
      hl_func = function(x) {
        x %>%
          dplyr::mutate(hl_var = ifelse(is.na(VE_threshold), 'none',
                                        as.character(VE_threshold))) %>%
          dplyr::mutate(hl_var = factor(hl_var,
                                        levels = c('none', '0', 
                                                   '1', '5'))) %>%
          { . }
      })
    map(thresholds, 
        ~exec(gen_depletion_tally_tables, 
              idxs = idxs,
              !!!base_settings,
              !!!analysis_settings,
              filter_group_cov = .x))
    map(thresholds, 
        ~exec(IE_lineplot, 
              idxs = idxs,
              !!!base_settings,
              !!!analysis_settings,
              filter_group_cov = .x))

    analysis_settings <- list(
      hl_lgd_name = 'Binding affinity rank percentile',
      x_var = 'percentile_rank',
      hl_name = 'affinity_threshold',
      hl_func = function(x) {
        x %>%
          dplyr::mutate(hl_var = ifelse(is.na(percentile_rank), 'none',
                                        as.character(percentile_rank))) %>%
          dplyr::mutate(hl_var = factor(hl_var,
                                        levels = c('none', '4', '3',
                                                   '1.9', '1'))) %>%
          { . }
      }
    )
    map(thresholds, 
        ~exec(gen_depletion_tally_tables, 
              idxs = idxs,
              !!!base_settings, 
              !!!analysis_settings, 
              filter_group_cov = .x))
    map(thresholds, 
        ~exec(IE_lineplot, 
              idxs = idxs,
              !!!base_settings,
              !!!analysis_settings,
              filter_group_cov = .x))
  })

