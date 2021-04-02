#!/usr/bin/env Rscript

#' Started in May 2019
#'
#' Compute the likelihood of presentation for each focus allele peptide by each
#' of a donor's other alleles. Arrive at a mean score per patient.

#' SETTINGS
check_cached <- F
min_coverage <- 0
max_lohhla_cov <- .25
ncores <- 32
ncores <- 16
ncores <- 1
ncores <- 38
ncores <- 42
ncores <- 24
ncores <- 8
debug_func <- F
redo <- F
hla_sim_method <- 'cutoff'
integration_method <- 'union'

## LOAD LIBRARIES
## Clean up BA tables first
# devtools::load_all(file.path('~/libs', 'quickMHC'))
# devtools::load_all(file.path('~/libs', 'maartenutils'))
# devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
source('~/antigenic_space/bin/init.R')
source(file.path(ma_dir, 'immune_editing', 'continuous_IE_detection_init.R'))
# source(file.path(ma_dir, 'immune_editing', 'HLA_presentation_scores_helpers.R'))
library(RPostgreSQL)
library(dbplyr)

doParallel::stopImplicitCluster()
doParallel::registerDoParallel(cores = ncores)

options(warn = 1)

# devtools::load_all(file.path('~/libs', 'quickMHC'))
# ncores = 1
if (F) {
  ## Make sure all alleles are indexed
  plyr::l_ply(c(all_AB_hlas, all_C_hlas), function(hla_allele) {
    message('Indexing ', hla_allele, ' if needed')
    bp <- quickMHC::BindingPredictor$new(hla_allele = hla_allele)
    bp$index_SQL()
    rm(bp)
  }, .parallel = (ncores > 1))
}


if (F) {
  ## Precompute all pairwise relationships
  for (focus_allele in focus_hlas) {
    plyr::l_ply(c(all_AB_hlas, all_C_hlas), function(ancillary_allele) {
      compute_pairwise_corroboration(hla_sim_method = hla_sim_method,
                                     focus_allele = focus_allele,
                                     ancillary_allele = ancillary_allele,
                                     pr = 1.9, pr_ancillary = 1.9)
    }, .parallel = (ncores > 1))
  }
}


# kill_all_connections()
hla_alleles = focus_hlas
hla_alleles = c('B1302', sample(all_AB_hlas))
for (hla_allele in hla_alleles) {
  new_PS <- compute_PS_overview(hla_alleles = hla_allele, mail_notify = T,
                                verbose = F,
                                # ncores = 1,
                                ncores = ncores,
                                redo = F) 
}

# saveRDS(new_PS, paste0(attr(new_PS, 'ro_fn'), '.bak'))
if (F) {
  ro_fn <- file.path(rds_dir,
                     glue::glue('focus_allele_avg_HLA_presentation-hla_sim\\
                                _method_-cutoff-integration_method_-union.rds'))
  # readRDS(ro_fn)
  new_ro <- controlled_merge(readRDS(ro_fn),
                             new_PS[, .(donor_id, hla, mean_score_AB)],
                             by_cols = c('donor_id', 'hla'))
  new_ro[mean_score == mean_score_AB & anchor_present == F & LOH_HLA == 'no_LOHHLA']
  new_ro[, LOH_HLA]
  # compute_PS_overview(hla_alleles = 'A3301', mail_notify = T, verbose = T,
  #                     ncores = ncores)
}
