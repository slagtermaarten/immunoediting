# focus_hlas <- c('A0201', 'A1101', 'B0702', 'B2705', 'B4001')
if (!exists('all_hlas')) {
  invisible(w_readRDS('optitype_tcga'))
  all_AB_hlas <-
    optitype_tcga[, unique(c(`hla_a1`, `hla_a2`, `hla_b1`, `hla_b2`))] %>%
    sort %>% quickMHC::shortenHLA()
  all_C_hlas <- optitype_tcga[, unique(c(`hla_c1`, `hla_c2`))] %>%
    sort %>% quickMHC::shortenHLA()
  all_hlas <- c(all_AB_hlas, all_C_hlas)
  hla_class_dat <- copy(w_readRDS('hla_class_dat'))
  setkey(hla_class_dat, donor_id)
}

