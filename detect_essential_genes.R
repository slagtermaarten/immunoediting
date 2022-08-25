devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
library(maartenutils)


parse_restifo_data <- function(FDR_thresh = .10) {
  fn <- file.path(ma_dir, 'immune_editing', 'data-raw', 
    'restifo_essential_genes.csv')
  fh <- fread(fn)
  fh[, sigFDR := as.numeric(gsub(',', '.', sigFDR))]
  fh <- fh[sigFDR <= FDR_thresh]
  return(fh)
}
parse_restifo_data(FDR_thresh = .00001)

# if (T) {
#   fh <- parse_restifo_data()
#   immunotherapy_ess_genes <- fh[, gene_id] %>% { .[. != ""] }
# } else {
#   immunotherapy_ess_genes <- c('JAK1', 'STAT1', 'TAP1', 'COL17A1', 'CD58',
#     'HLA-A', 'TAPBP', 'TAF3', 'SRP54', 'B2M', 'RPL23', 'SOX10')
# }

merge_cols <- c('variant_id', 'chromosome', 'start_position', 'ref_allele',
                'alt_allele')


#' Compile an overview of mutations in genes essential to T-cell functioning
#'
#'
detect_essential_genes <- function(donor_id,
  cols_of_interest = c('impact', 'hugo_symbol', 'variant_classification', 'vaf',
                       'AT_MLE_CCF', 'nmd_remark'),
  FDR_thresh = .1) {

  message(donor_id)

  fh <- read_input_variants(donor_id)
  if (null_dat(fh)) return(NULL)
  restifo_fh <- parse_restifo_data(FDR_thresh = FDR_thresh)
  immunotherapy_ess_genes <- restifo_fh[, gene_id]
  # grep('hsa', immunotherapy_ess_genes, invert = T, value = T)
  fh <- fh[hugo_symbol %in% immunotherapy_ess_genes &
           variant_classification != 'silent_mutation']

  if (F) {
    fh_a <- read_raw_maf(donor_id)
    if (!null_dat(fh_a)) {
      fh_a <- fh_a[hugo_symbol %in% immunotherapy_ess_genes]
      fh_a <- fh_a[!grepl('silent', variant_classification, ignore.case = T)]
    } else {
      warning('no raw maf available ', donor_id)
    }

    if (!null_dat(fh_a) && !null_dat(fh)) {
      if (nrow(fh_a) > nrow(fh)) {
        setnames(fh_a,
                 c('chromosome', 'start_position', 'reference_allele',
                   'tumor_seq_allele1'),
                 merge_cols[2:5])

        fh_a[, chromosome := as.character(chromosome)]
        fh[, chromosome := as.character(chromosome)]
        setkeyv(fh_a, merge_cols[2:5])
        setkeyv(fh, merge_cols[2:5])
        # classes <- fh[, lapply(.SD, class), .SDcols = merge_cols]
        # classes <- fh_a[, lapply(.SD, class), .SDcols = merge_cols]
        # classes <- unlist(classes)
        # fh <- maartenutils::set_dt_types(fh, classes)
        # fh_a <- maartenutils::set_dt_types(fh_a, classes)
        dtf <- fh_a
        ## Deduplicate columns
        dtf <- dtf[, -which(duplicated(colnames(dtf))), with = F]
        fh_m <- merge(dtf, fh, all.x = T)
        dtf[, merge_cols[2:5], with = F]
        fh[, merge_cols[2:5], with = F]

        if (fh_a[fh, .N] != row(fh)) {
          browser()
          merge(fh, fh_a, by = c('chromosome', 'start_position', 'ref_allele',
                                 'alt_allele'))
          fh_a[, variant_type]
          fh_a[, variant_classification]
          missing_vars <- setdiff(fh_a[, variant_id], fh[, variant_id])
          fh_a[variant_id %in% missing_vars]
        }
      }
    }
  }
  # browser(expr = donor_id == 'TCGA-CJ-4888')
  cols_of_interest <- intersect(cols_of_interest, colnames(fh))
  ret_val <- fh[, cols_of_interest, with = F]
  ret_val[, 'donor_id' := donor_id]
  return(ret_val)
}

ncores <- 32
doParallel::registerDoParallel(cores = ncores)
donor_ids <- intersect(setdiff(core_donors, 'TCGA-CJ-4888'),
  pipeline_out_fns[project %in% tcga_projects, donor_id])

if (T) {
} else {
  extra_extra_extra_lenient_escape_prev_ext <- rbindlist(
    plyr::llply(donor_ids, detect_essential_genes,
      .parallel = !maartenutils::local_run, FDR_thresh = .00001), 
    fill = T)
  w_saveRDS('extra_extra_extra_lenient_escape_prev_ext')

  extra_extra_lenient_escape_prev_ext <- rbindlist(
    plyr::llply(donor_ids, detect_essential_genes,
      .parallel = !maartenutils::local_run, FDR_thresh = .0001), 
    fill = T)
  w_saveRDS('extra_extra_lenient_escape_prev_ext')

  extra_lenient_escape_prev_ext <- rbindlist(
    plyr::llply(donor_ids, detect_essential_genes,
      .parallel = !maartenutils::local_run, FDR_thresh = .001), 
    fill = T)
  w_saveRDS('extra_lenient_escape_prev_ext')

  escape_prev_ext <- rbindlist(
    plyr::llply(donor_ids, detect_essential_genes,
      .parallel = !maartenutils::local_run), fill = T)
  w_saveRDS('escape_prev_ext')

  lenient_escape_prev_ext <- rbindlist(
    plyr::llply(donor_ids, detect_essential_genes,
      .parallel = !maartenutils::local_run, FDR_thresh = .01), 
    fill = T)
  w_saveRDS('lenient_escape_prev_ext')
}
