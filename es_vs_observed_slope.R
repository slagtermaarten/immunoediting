#' SETTINGS
output_dir <- '~/antigenic_space/maarten-analyses/img/cont_IE_heatmaps_V2'
ro_fn <- file.path(rds_dir, glue::glue('focus_allele_avg_HLA_presentation-hla_sim\\
                                 _method_-cutoff-integration_method_\\
                                 -union.rds'))
repertoire_overlap_dat <- readRDS(ro_fn)
fn_suffix = '-im_union'

es_vs_observed_slope <- function(tumor_type = "Pan*'-'*cancer",
                                 restrict_to_powered_analyses = T,
                                 return_res = F,
                                 p_val_bound = 0.25,
                                 PPV = .45,
                                 subset_perc = 1,
                                 power_var = '0.8',
                                 adaptive_p_val_bound = T,
                                 img_dir = output_dir,
                                 do_test_grobs = F) {
  messagef('Starting %s', tumor_type)
  # obj_fn <- glue::glue('~/Projects/antigenic_space/maarten-analyses\\
  #     /rds/ch_test_{tumor_type}.rds')
  # if (!file.exists(obj_fn)) return(NULL)
  # dtf <- readRDS(obj_fn)

  dtf <- prep_pan_IE_heatmap(
    tumor_type = tumor_type,
    fn_suffix = fn_suffix,
    repertoire_overlap_dat = repertoire_overlap_dat,
    p_val_bound = p_val_bound,
    PPV = PPV,
    adaptive_p_val_bound = T
    # redo = redo
    # PPV = PPV,
    # fill_var = fill_var,
    # ncores = ncores,
    # adaptive_p_val_bound = adaptive_p_val_bound
  )
  browser()

  if (null_dat(dtf)) return(NULL)
  stopifnot(is.numeric(subset_perc) && subset_perc <= 1 && subset_perc > 0)
  if (!is.null(p_val_bound)) {
    dtf[p_val.adj > p_val_bound, signif_effect.adj := NA]
  }

  library(ggplot2)
  levs <- c('-', 'L', 'G')
  p_dat <- dtf %>%
    mutate(shape_var = ifelse(is.na(signif_effect.adj),
        '-', signif_effect.adj)) %>%
    mutate(shape_var = factor(shape_var, levels = levs)) %>%
    setDT
  text_labels <- c("'Insignificant'", "'Significantly'<0", "'Significantly'>0")

  shape_labels <- p_dat[, .N, shape_var] %>%
    setkey(shape_var) %>%
    .[c('-', 'L', 'G'), ] %>%
    .[, N := ifelse(is.na(N), 0, N)] %>%
    .[, sprintf('%s~(italic(n)==%d)', text_labels, N)] %>%
    sapply(function(x) parse(text = x)) %>%
    setNames(NULL)
  p_dat %<>%
    # mutate(shape_var = factor(shape_var, levels = levs, labels = shape_labels)) %>%
    mutate(shape_var = factor(shape_var, levels = levs)) %>%
    arrange(shape_var) %>%
    setDT
  cols <- c('grey80', '#A30D1D', '#3F9CB5')
  point_name <- 'Test result (FDR-corrected)'

  p_dat[, 'es' := get(power_var)]
  p <- ggplot(p_dat, aes_string(x = 'es', y = 'rc', shape = 'shape_var',
                                colour = 'shape_var', fill = 'shape_var')) +
    geom_hline(yintercept = 0, col = 'grey80', linetype = 'dashed') +
    geom_point(alpha = .5) +
    scale_shape_manual(name = point_name,
      breaks = p_dat[, levels(shape_var)],
      labels = shape_labels, values = c(16, 25, 24)) +
    scale_colour_manual(name = point_name,
      breaks = p_dat[, levels(shape_var)],
      labels = shape_labels, values = cols) +
    scale_fill_manual(name = point_name,
      breaks = p_dat[, levels(shape_var)],
      labels = shape_labels, values = cols) +
    # scale_shape_manual(name = point_name,
    #   breaks = shape_labels, values = c(16, 25, 24)) +
    # scale_colour_manual(name = point_name,
    #   breaks = shape_labels, values = cols) +
    # scale_fill_manual(name = point_name,
    #   breaks = shape_labels, values = cols) +
    theme_ms(legend.position = c(.5, 0.95),
      legend.direction = 'vertical',
      legend.justification = c(.5, 1)) +
    scale_x_continuous(name = 'Required fractional decrease\nbetween PS=0 and PS=1\n for statistical significance',
      limits = c(0, 1),
      # trans = 'log10', 
      # labels = fancy_scientific) +
      labels = function(x) ifelse(x == 0, x, paste0('-', x))) +
    scale_y_continuous(name = 'Observed slope')

  if (return_res) {
    return(p)
  } else {
    o_fn <- sprintf('%s/es_vs_observed_slope-%s.pdf', img_dir, tumor_type)
    ggsave(o_fn, p, width = 8.7, height = 8, unit = 'cm')
    return(o_fn)
  }
}

if (sys.nframe() %in% c(0, 4)) {
  expand.grid(tumor_type = tumor_types) %>%
  # expand.grid(tumor_type = 'melanoma') %>%
    mutate(idx = 1:n()) %>%
    plyr::a_ply(1, function(r) {
      idx <- as.integer(r[['idx']])
      messagef('Processing index %d', idx)
      tryCatch(es_vs_observed_slope(tumor_type = as.character(r[['tumor_type']])),
               error = function(e) { print(e) })
    })
}
