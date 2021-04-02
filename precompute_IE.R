#!/usr/bin/env Rscript

#' This file computes all continuous IE stats objects we require

source('~/antigenic_space/bin/init.R')
source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_init.R'))
source('~/libs/maartenutils/R/0-general_utils.R')

#' SETTINGS
ncores <- 8
ncores <- 8
ncores <- 6
ncores <- 4
ncores <- 38
ncores <- 40
ncores <- 12
ncores <- 24
ncores <- 1
ncores <- 36
options(error = recover)
options(error = traceback)
# test <- function() { stop() }
# test()
redo <- T
redo <- F
reg_method = 'rlm_one_group'
debug_master <- T
debug_master <- F
check_res <- F
fn_suffix = '-im_union'
# ro_fn <- file.path(rds_dir, 'focus_allele_avg_HLA_presentation.rds')
ro_fn <- file.path(rds_dir, glue('focus_allele_avg_HLA_presentation-hla_sim\\
                                 _method_-cutoff-integration_method_\\
                                 -union.rds'))
repertoire_overlap_dat <- readRDS(ro_fn)
repertoire_overlap_dat <- NULL

hla_alleles_to_test = focus_hlas
hla_alleles_to_test = 'A3301'

# class(requires_computation)
# print(requires_computation)
# environment(requires_computation)$minimum_mod_time

doParallel::registerDoParallel(cores = ncores)

projects_that_need_recomputation <-
  donor_summary[project %in% virus_candidate_projects] %>%
  .[, unique(as.character(project_extended))] %>%
  c('Pan*\'-\'*cancer') %>%
  as.character %>%
  grep('Stomach.*EBV.*MSI', ., value = T, invert = T)

all_projects <-
  donor_summary[, unique(as.character(project_extended))] %>%
  c('Pan*\'-\'*cancer') %>%
  as.character %>%
  grep('Stomach.*EBV.*MSI', ., value = T, invert = T)

if (F) {
  settings <-
    as.list(all_cont_IE_settings) %>%
    lapply(unique) %>%
    {
      do.call(expand.grid,
              c(., list('project_extended' = projects_that_need_recomputation,
                        stringsAsFactors = F)))
    }
} else if (F) {
  settings <- all_cont_IE_settings
} else {
  settings <-
    as.list(all_cont_IE_settings) %>%
    lapply(unique) %>%
    { .[names(.) != 'focus_allele'] } %>%
    {
      do.call(expand.grid,
              c(., list(focus_allele = 'A3301', stringsAsFactors = F)))
    }
}


if (T) {
  ## phase 1: perform all regression analyses
  ## iterate over all ie settings
  plyr::l_ply(1:nrow(settings), function(i) {
  # plyr::l_ply(23, function(i) {
    focus_allele <- settings[i, 'focus_allele']
    analysis_name <- settings[i, 'analysis_name']
    if (F && analysis_name != 'twoD_sens_analysis') return(NULL)
    overlap_var <- settings[i, 'overlap_var']
    # hla_sim_range <- settings[[i, 'hla_sim_range']]
    LOH_HLA <- settings[i, 'LOH_HLA']
    patient_inclusion_crit <- unlist(settings[i, 'patient_inclusion_crit'])
    pe <- settings[i, 'project_extended']
    mymessage(glue::glue('Processing {i}/{nrow(settings)}'), 'main loop')
    if (check_res) {
      it_func <- plyr::llply
    } else {
      it_func <- plyr::l_ply
    }
    res <- it_func(main_analysis_idxs, function(j) {
      test_continuous_IE(
        analysis_name = analysis_name,
        idx = j,
        reg_method = 'rlm_one_group',
        tumor_type = pe,
        patient_inclusion_crit = patient_inclusion_crit,
        repertoire_overlap_dat = repertoire_overlap_dat,
        LOH_HLA = LOH_HLA,
        check_res = check_res,
        fn_suffix = fn_suffix,
        # hla_sim_range = hla_sim_range,
        focus_allele = focus_allele,
        ncores = 1,
        return_res = F,
        redo = redo,
        overlap_var = overlap_var)
    }, .parallel = (ncores > 1))
    if (check_res) {
      print(i)
      res <- setNames(res, analysis_idxs)
      print(res)
    }
  }, .parallel = F)
}


## Phase 1.5
## Summarize the analyses by representing each analysis with the most
## representative allele
if (T) {
  source(file.path(ma_dir, 'immune_editing',
                   'continuous_IE_detection_init.R'))
  res <- compile_all_coef_overview(redo = redo,
                                   ncores = ncores,
                                   hla_alleles = hla_alleles_to_test,
                                   repertoire_overlap_dat = repertoire_overlap_dat,
                                   fn_suffix = fn_suffix)
  if (null_dat(res))
    stop('Run and complete compile_all_coef_overview first')
}


# prep_pan_IE_heatmap(tumor_type = 'Pan*\'-\'*cancer',
#                     redo = redo,
#                     fill_var = 'rc',
#                     ncores = ncores,
#                     return_res = F,
#                     fn_suffix = fn_suffix,
#                     p_val_bound = .25,
#                     PPV = .45,
#                     adaptive_p_val_bound = T)

source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_init.R'))
# p_projects <- 'Prostate'
p_projects <- all_projects
p_projects <- 'Pan*\'-\'*cancer'
redo = T
if (T) {
  ## Phase 2: compile all regression analyses into tumor-specific overviews,
  ## summarize the analyses by representing each analysis with the most
  ## representative allele, compute statistical power for the combination of
  ## this allele and analysis settings
  p_projects %>% plyr::l_ply(function(pe) {
    prep_pan_IE_heatmap(tumor_type = pe,
                        redo = redo,
                        fill_var = 'rc',
                        ncores = ncores,
                        return_res = F,
                        fn_suffix = fn_suffix,
                        repertoire_overlap_dat = repertoire_overlap_dat,
                        p_val_bound = .25,
                        PPV = .45,
                        adaptive_p_val_bound = T)
    }, .parallel = F)
}

if (F) {
  p_projects %>%
    plyr::l_ply(function(pe) {
    prep_pan_IE_heatmap(tumor_type = pe,
                        redo = redo,
                        fill_var = 'rc_CYT',
                        ncores = ncores,
                        return_res = F,
                        fn_suffix = fn_suffix,
                        repertoire_overlap_dat = repertoire_overlap_dat,
                        p_val_bound = .25,
                        PPV = .45,
                        adaptive_p_val_bound = T)
    }, .parallel = F)
}

# mail_notify(subject = 'IE analyses', msg = 'Finished IE precompute')
