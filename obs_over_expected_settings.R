compute_sig_score <- function(v) {
  tab <- prespecified_table(v, exp_levels = c('', 'G', 'L'))
  tab['L'] - tab['G']
}


get_unobs_obs_settings <- function(tt) {
  process_settings <- function(dtf) {
    dplyr::select(dtf, analysis_name, overlap_var,
                  patient_inclusion_crit, LOH_HLA, analysis_idx) %>%
    unique %>%
    { cbind(., ds_param_grid[.$analysis_idx, ]) } %>%
    mutate(analysis_idx = NULL) %>%
    mutate(exp_filter = ifelse(!is.na(expression_threshold),
                                      paste0('GExp=', expression_threshold),
                                      ifelse(!is.na(VE_threshold),
                                             paste0('VExp=', VE_threshold),
                                             'none'))) %>%
    mutate(expression_threshold = NULL) %>%
    mutate(VE_threshold = NULL)
  }

  if (!is.null(tt)) {
    dtf <- all_preps %>% filter(project_extended == tt)
  } else {
    dtf <- all_preps
  }
  ## All combinations of tumor types and settings for which at least
  ## one tumor
  ## type is significantly edited
  all_settings <- process_settings(dtf)
  obs_settings <- dtf %>%
    dplyr::filter(p_val_signif_effect.adj == 'L') %>%
    process_settings
  setDT(all_settings)
  setDT(obs_settings)
  setkeyv(all_settings, colnames(all_settings))
  setkeyv(obs_settings, colnames(all_settings))
  ## Don't know how to do easily do this in the tidyverse
  unobs_settings <- all_settings[!obs_settings]
  stopifnot(nrow(unobs_settings) + nrow(obs_settings) == nrow(all_settings))
  return(list('obs_settings' = obs_settings, 'unobs_settings' = unobs_settings))
}


gen_forest_plot <- function(tt = NULL, fontsize = 6, return_gtab = F) {
  ## List of dataframes
  lod <- get_unobs_obs_settings(tt)
  obs_settings <- lod$obs_settings
  unobs_settings <- lod$unobs_settings
  if (null_dat(obs_settings)) return(NULL)

  if (F) {
    trans <- sqrt
    trans_back <- function(x) 2^x
  } else {
    trans <- identity
    trans_back <- identity
  }

  enrich_dtf <- map_dfr(auto_name(colnames(obs_settings)),
                        function(varn = 'overlap_var') {
    obs_f <- dplyr::pull(obs_settings, !!varn)
    unobs_f <- dplyr::pull(unobs_settings, !!varn)
    all_levs <- c(obs_f, unobs_f) %>% unique
    # if (varn == 'patient_inclusion_crit') browser()
    if (FALSE && length(all_levs) == 2) {
      all_levs <- all_levs[-1]
    }
    res <- all_levs %>%
      auto_name %>%
      map_dfr(function(lev) {
        ## Current level vs. other levels
        d <- data.frame(c_lev=c(sum(obs_f == lev, na.rm = T),
                                sum(unobs_f == lev, na.rm = T)),
                        o_lev=c(sum(obs_f != lev, na.rm = T),
                                sum(unobs_f != lev, na.rm = T)))
        row.names(d) <- c('obs', 'unobs')
        test <- fisher.test(d, alternative = 'two.sided')
        stars_s <- gen_p_val_stars(test$p.value, max_stars = 5,
                                   logs_per_star = 3)
        list('mean' = test$estimate,
             'lower' = test$conf.int[1],
             'upper' = test$conf.int[2],
             'OR_string' = glue('{d$c_lev[1]}/{d$o_lev[1]} vs. ',
                                '{d$c_lev[2]}/{d$o_lev[2]} {stars_s}'),
             'p' = test$p.value)
      }, .id = 'level')
    res$level <- as.character(all_levs)
    if (varn %in% names(display_settings)) {
      res %<>%
        mutate(level = display_settings[[varn]]$labels %>%
               { .[match(as.character(level), display_settings[[varn]]$breaks)] })
    }
    return(res)
  }, .id = 'variable') %>%
  dplyr::group_by(variable) %>%
  dplyr::arrange(-mean, .by_group = T) %>%
  mutate(variable = var_names[variable])

  row_height <- unit(1, 'cm')
  s_size <- unit(.25, 'cm')
  s_size <- unit(.1, 'cm')
  s_size <- unit(.5, 'mm')

  CI_range <- c(min(enrich_dtf$lower %>% { .[is.finite(.)] }, na.rm = T),
                max(enrich_dtf$upper %>% { .[is.finite(.)] }, na.rm = T)) %>%
    trans
  pos_within_scale <- function(x, margin = .1) {
    unit(margin + (1 - margin) * (trans(x) - CI_range[1]) / diff(CI_range),
         'npc')
  }
  plot_vp <- viewport(height = row_height)

  idx <- 1
  plot_grobs <- enrich_dtf %>%
    ungroup %>%
    select(lower, mean, upper) %>%
    pmap(function(lower, mean, upper) {
      if (any(!is.finite(c(lower, mean, upper)))) {
        textG <- textGrob('Inf',
                          x = unit(.5, 'npc'),
                          y = unit(.5, 'npc'),
                          gp = gpar(fontsize = fontsize),
                          vjust = .5, hjust = .5)
        vline <- linesGrob(x = pos_within_scale(c(1, 1)),
                           y = unit(c(0, 1), 'npc'),
                           gp = gpar(fill = 'grey80'))
        grobs <- list(textG, vline)
      } else {
        myrect <- rectGrob(x = pos_within_scale(mean),
                           width = sqrt(s_size),
                           height = sqrt(s_size),
                           gp = gpar(fill = 'black'))
        myline <- linesGrob(x = pos_within_scale(c(lower, upper)),
                            y = unit(.5, 'npc'),
                            gp = gpar(fill = 'black'))
        vline <- linesGrob(x = pos_within_scale(c(1, 1)),
                           y = unit(c(0, 1), 'npc'),
                            gp = gpar(fill = 'grey80'))
        grobs <- list(myrect, myline, vline)
      }
      if (idx == nrow(enrich_dtf)) {
        breaks <- seq(CI_range[1], CI_range[2], length.out = 3)
        xaxis <- tryCatch(xaxisGrob(at = breaks,
                                    # label = trans_back(breaks),
                                    gp = gpar(fontsize = fontsize - 2)),
                          error = function(e) { print(e); browser() })
        # grobs <- c(grobs, list(xaxis))
      }
      idx <<- idx + 1
      do.call(grobTree, grobs)
    })
  # grid.draw(plot_grobs[[3]])
  # grid.draw(plot_grobs[[length(plot_grobs) - 2]])
    # browser()
  # grid.draw(plot_grobs[[length(plot_grobs)]])

  library(gtable)
  format_label <- function(label, align = 'left') {
    if (align == 'left') {
      textGrob(label,
               x = unit(0, 'npc') + unit(1, 'mm'),
               y = unit(1, 'npc') - unit(1, 'mm'),
               gp = gpar(fontsize = fontsize),
               vjust = 1, hjust = 0)
    } else if (align == 'right') {
      textGrob(label,
               x = unit(1, 'npc') - unit(1, 'mm'),
               y = unit(1, 'npc') - unit(1, 'mm'),
               gp = gpar(fontsize = fontsize),
               vjust = 1, hjust = 1)
    }
  }
  grob_matrix <-
    list(map(sub_repetitions(unlist(enrich_dtf[, 1])),
             ~format_label(.x)),
         map(sub_repetitions(unlist(enrich_dtf[, 2])),
             ~format_label(.x, align = 'right')),
         map(sub_repetitions(unlist(enrich_dtf$OR_string)),
             ~format_label(.x, align = 'right')),
         plot_grobs) %>%
    unlist(recursive = F) %>%
    { matrix(., ncol = 4, byrow = F) }
  row_heights <- unit(rep(.35, length(plot_grobs)), 'cm')
  # grid::heightDetails(grob_matrix[[4, 1]])
  # browser()
  column_widths <- unit(c(1.9, 2.6, 2.4, 1.8), 'cm')
  sum(as.numeric(column_widths))
  gtab <- gtable_matrix(name = 'forest',
                        grobs = grob_matrix,
                        clip = "off",
                        widths = column_widths, heights = row_heights)
  grey_bg <- T
  class_labels <- map_chr(1:nrow(grob_matrix), ~grob_matrix[[.x]]$label)
  change_locs <- which(class_labels != '')
  for (i in 1:nrow(enrich_dtf)) {
     if (i %in% change_locs)
       grey_bg <- !grey_bg
     bg_col <- ifelse(grey_bg, 'grey95', 'white')
     gtab %<>% gtable_add_grob(
       rectGrob(gp = gpar(fill = bg_col, col = NA, lwd = 0)),
       t = i, b = i, l = 1, r = 4,
       z = 0, clip = "on", name = 'bg')
  }
  column_titles <- map(c('Variable', 'Level', 'OR', ''),
    ~grobTree(rectGrob(gp = gpar(fill = 'grey40', lwd = 0, col = 'grey40')),
              textGrob(label = .x,
                       gp = gpar(color = 'white',
                                 col = 'white',
                                 fontsize = 6,
                                 fontface = 'italic')))) %>%
    { matrix(., ncol = 4, byrow = F) } %>%
    { gtable_matrix(name = 'col_titles',
                    grobs = .,
                    clip = "off",
                    widths = column_widths, heights = unit(.5, 'cm')) }
  gtab <- rbind(column_titles, gtab)
  ## Make 3rd label span 3rd and 4th columns
  gtab$layout[3, ]$r <- 4
  gtab$grobs[[4]] <- nullGrob()

  if (!return_gtab) {
    plot_width <- unit(8.7, 'cm')
    plot_height <- unit(12, 'cm')
    tt_simple <- tumor_types[tt]
    cat(tt_simple, '\n')
    o_fn <- glue::glue('{output_dir}/{tt_simple}.pdf')
    pdf(o_fn,
        width = as.numeric(convertUnit(plot_width, 'inch')),
        height = as.numeric(convertUnit(plot_height, 'inch')))
    plot(gtab)
    dev.off()
  } else {
    return(gtab)
  }

  if (F) {
    f_labels <- map(enrich_dtf[, c(1, 2)], sub_repetitions) %>%
      c(list('p' = map(fancy_scientific(enrich_dtf$p), ~parse(text = .x))))
    f_labels[[1]] <- c('Variable name', f_labels[[1]])
    f_labels[[2]] <- c('Value', f_labels[[2]])
    f_labels[[3]] <- c('p-value (Fisher-exact)', f_labels[[3]])

    if (!require('forestplot')) install.packages('forestplot')
    library('forestplot')
    # debugonce(forestplot:::forestplot.default)
    source(file.path(ma_dir, 'immune_editing', 'forest_plot.R'))
    forestplot(labeltext = f_labels,
               mean = c(NA, enrich_dtf$mean),
               upper = c(NA, enrich_dtf$upper),
               lower = c(NA, enrich_dtf$lower),
               colgap = unit(4, 'mm'),
               zero = 1,
               # txt_gp = fpTxtGp(label = gpar(),
               #                  summary = gpar(),
               #                  xlab = gpar(),
               #                  title = gpar(),
               #                  ticks = gpar(),
               #                  legend = gpar(),
               #                  legend.title = gpar(),
               #                  cex = 2),
               new_page = TRUE,
               is.summary = c(T, rep(F, length(f_labels[[1]]) - 1)),
               xlog = F,
               col = fpColors(box = 'royalblue', line = 'darkblue',
                              summary='royalblue'))
  }
}



#' Prep variables for plotting in one of multiple downstream functions
#'
#'
prep_IE_plot <- function(hl_name = '',
                         idxs = main_analysis_idxs,
                         hl_lgd_name = '',
                         x_var = 'percentile_rank',
                         p_var = 'p_val.adj',
                         filter_qcd = NULL,
                         filter_group_qcd = NULL,
                         filter_cov = NULL,
                         filter_group_cov = NULL,
                         filter_p = NULL,
                         filter_group_p = NULL,
                         include_non_ds_tallies = F,
                         filter_depletion_outliers = F,
                         return_mode = 'preps',
                         # tumor_type = "Pan*'-'*cancer",
                         z_normalize = F,
                         ncores = 1,
                         tumor_type = NULL,
                         filter_single_x_var_projects = T,
                         hl_func = function(x)
                           dplyr::mutate(x, hl_var = factor(''))) {
  if (return_mode == 'preps') {
    all_preps <- get_all_preps(idxs = idxs,
                               ncores = ncores,
                               include_non_ds_tallies = include_non_ds_tallies,
                               z_normalize = z_normalize)
    t_dat <- all_preps %>%
      ## Restrict to analyses where depletion is located in the range [-1, 1]
      # dplyr::filter(depletion_full >= -1 & depletion_full <= 1) %>%
      { cbind(., ds_param_grid[.$analysis_idx, ]) } %>%
      hl_func %>%
      { . }
  } else if (return_mode == 'overview') {
    ## For the lineplot, we want all alleles rather than the median allele only,
    ## which is why we don't use *preps, but *overviews
    t_dat <- compile_all_coef_overview(
      redo = F,
      ncores = ncores,
      z_normalize = z_normalize,
      include_non_ds_tallies = include_non_ds_tallies,
      idxs = idxs,
      repertoire_overlap_dat = repertoire_overlap_dat) %>%
      { cbind(., ds_param_grid[.$analysis_idx, ]) } %>%
      hl_func
    t_dat$hl_var
  }
  orig_row_count <- nrow(t_dat)
  # t_dat %>% filter(rc < 0) %>% pull(depletion_full) %>% summary %>% print
  # t_dat %>% filter(rc > 0) %>% pull(depletion_full) %>% summary %>% print

  ## Group analyses by all relevant factors except for the x_var
  id_vars <- c('project_extended', 'overlap_var',
               'patient_inclusion_crit', 'LOH_HLA',
               'analysis_name', 'focus_allele',
               colnames(ds_param_grid)) %>% setdiff(x_var)
  t_dat <- t_dat %>% {
      dplyr::mutate(.,
                    group_var = apply(select(., all_of(id_vars)), 1,
                                      function(r) paste0(r, collapse = '-')))
    } %>%
    dplyr::mutate(color_var = factor(focus_allele)) %>%
    dplyr::mutate(x_var = factor(as.character(.data[[x_var]]))) %>%
    # dplyr::filter(depletion_full >= -1 & depletion_full <= 1) %>%
    dplyr::filter(analysis_idx %in% idxs) %>%
    hl_func %>%
    { . }

  if (!is.null(tumor_type)) {
    t_dat <- t_dat %>% dplyr::filter(project_extended == tumor_type)
  }
  ## Filter down to robust analyses in a variety of ways
  if (!is.null(filter_qcd) && !is.na(filter_qcd)) {
    t_dat <- t_dat %>% dplyr::filter(abs(MCMC_qcd) <= filter_qcd)
    if (nrow(t_dat) == 0)
      mystop(msg = glue('No analyses remaining after filter_qcd'))
  }
  if (!is.null(filter_cov) && !is.na(filter_cov)) {
    t_dat <- t_dat %>% dplyr::filter(abs(MCMC_cov) <= filter_cov)
    if (nrow(t_dat) == 0)
      mystop(msg = glue('No analyses remaining after filter_cov'))
  }
  if (!is.null(filter_p) && !is.na(filter_p)) {
    t_dat <- t_dat <- dplyr::filter(p_val <= filter_p)
    if (nrow(t_dat) == 0)
      mystop(msg = glue('No analyses remaining after filter_p'))
  }

  if (!is.null(filter_group_cov) && !is.na(filter_group_cov)) {
    t_dat <- t_dat %>%
      dplyr::filter(abs(MCMC_cov) <= filter_group_cov)
    ## Select groups for which all analyses are robust
    clean_groups <- t_dat %>%
      dplyr::group_by(group_var, color_var) %>%
      dplyr::summarize(N = n()) %>%
      dplyr::filter(N == nlevels(t_dat$hl_var)) %>%
      dplyr::pull(group_var)
    if (length(clean_groups) == 0) {
      mywarning('No groups remaining, try a lower value of filter_group_cov',
           call. = F)
      return(NULL)
    } else {
      t_dat <- t_dat %>% filter(group_var %in% clean_groups)
    }
  }

  if (!is.null(filter_group_qcd) && !is.na(filter_group_qcd)) {
    t_dat <- t_dat %>%
      dplyr::filter(abs(MCMC_qcd) <= filter_group_qcd)
    ## Select groups for which all analyses are robust
    clean_groups <- t_dat %>%
      dplyr::group_by(group_var, color_var) %>%
      dplyr::summarize(N = n()) %>%
      dplyr::filter(N == nlevels(t_dat$x_var)) %>%
      dplyr::pull(group_var)
    if (length(clean_groups) == 0) {
      stop('No groups remaining, try a lower value of filter_group_qcd',
           call. = F)
    } else {
      t_dat <- t_dat %>% filter(group_var %in% clean_groups)
    }
  }

  if (!is.null(filter_group_p) && !is.na(filter_group_p)) {
    group_p <- t_dat %>%
      dplyr::group_by(group_var) %>%
      dplyr::summarize(min_p = min(p_val, na.rm = T)) %>%
      dplyr::filter(min_p <= filter_group_p) %>%
      dplyr::pull(group_var) %>%
      unique
    if (length(group_p) == 0) {
      stop('No groups remaining, try a higher value of filter_group_p',
           call. = F)
    } else {
      t_dat <- t_dat %>% filter(group_var %in% group_p)
    }
  }

  if (filter_depletion_outliers) {
    t_dat_f <- remove_outliers(t_dat, by_cols = 'project_extended',
                               test_cols = 'depletion_full',
                               probs = c(.05, .95))
    t_dat <- t_dat_f
  }

  ## Remove tumor types for which only one level of x_var/hl_var remains
  if ((is.na(tumor_type) || is.null(tumor_type)) &&
      filter_single_x_var_projects) {
    inventorization <- t_dat %>%
      group_by(project_extended) %>%
      dplyr::summarize('pass' = length(unique(x_var)) > 1)
    allowed_projects <-
      inventorization %>% filter(pass == T) %>% pull(project_extended)
    t_dat <- t_dat %>% filter(project_extended %in% allowed_projects)
  }

  cat(paste0(format(100 * nrow(t_dat) / orig_row_count, digits = 8), '%'),
      'of all analyses retained \n')

  if (is.null(tumor_type)) {
    cat('Overview of medians per project: \n')
    t_dat %>% group_by(project_extended, hl_var) %>%
      dplyr::summarise('med' = median(depletion_full)) %>%
      group_by(project_extended) %>%
      dplyr::summarise('med' = median(med)) %>%
      arrange(-med) %>%
      print(n = 50)
  }

  project_ordering <- t_dat %>%
    dplyr::group_by(project_extended) %>%
    # dplyr::summarize(order_val = sum(.data[[p_var]] < 1, na.rm = T)) %>%
    # dplyr::filter(order_val > 0) %>%
    dplyr::summarize(order_val = mean(.data[[p_var]], na.rm = T)) %>%
    # dplyr::summarize(order_val = mean(!!p_var, na.rm = T)) %>%
    dplyr::arrange(-order_val) %>%
    pull(project_extended)
  t_dat <- t_dat %>%
    dplyr::filter(project_extended %in% project_ordering) %>%
    dplyr::mutate(project_extended = factor(project_extended,
                                            levels = project_ordering))

  man_colors <- fas_discrete_colors(levels(t_dat$color_var))
  other_levs <- c('other', 'none')
  man_colors[names(man_colors) %in% other_levs] <- 'grey80'

  # t_dat$hl_var %>% fct_relevel(other_levs, after = 0)
  p_val_mod_string <- ifelse(grepl('adj', p_var), 'adj. ', '')

  fn <- glue('{{plot_basename}}',
             '-hl_name={hl_name}',
             '{make_flag(filter_group_cov)}',
             '{make_flag(filter_group_qcd)}',
             '{make_flag(filter_group_p)}',
             '{make_flag(filter_cov)}',
             '{make_flag(filter_qcd)}',
             '{make_flag(filter_p)}',
             '{make_flag(filter_depletion_outliers)}',
             '{make_flag(z_normalize)}',
             '-idxs={attr(idxs, \'name\')}.pdf') %>%
    { file.path(img_loc, .) }
  ## Conventiently return everything in the current scope ('environment')
  return(as.list(environment()))
}


IE_volcano <- function() {
  calling_args <- as.list(environment())
  with(do.call(prep_IE_plot, calling_args), {
    if (null_dat(t_dat)) return(NULL)
    p <- t_dat %>%
      dplyr::arrange(hl_var) %>%
      dplyr::mutate(y_var = -log10(.data[[p_var]])) %>%
      ggplot(aes(x = depletion_full, y = y_var)) +
        # geom_point() +
        geom_hex(bins = 100) +
        scale_y_continuous(name = glue('-log10({p_val_mod_string}p-value)')) +
        scale_x_continuous(name = glue('Neo-antigen depletion between\n',
                                       'PS = 0 and PS = 1')) +
        scale_colour_manual(name = hl_lgd_name, values = man_colors) +
        theme_fas() +
        facet_wrap(~project_extended, labeller = label_parsed, scale = 'free')
    plot_basename <- 'volcano'
    fn <- glue(fn)
    ggsave(filename = fn, plot = p, height = 18, width = 17.4, unit = 'cm')
    return(invisible())
  })
}
formals(IE_volcano) <- formals(prep_IE_plot)


IE_boxplot <- function() {
  calling_args <- as.list(environment())
  with(do.call(prep_IE_plot, calling_args), {
    if (null_dat(t_dat)) return(NULL)
    p <- t_dat %>%
      dplyr::arrange(hl_var) %>%
      ggplot(aes(x = depletion_full, y = hl_var, fill = hl_var)) +
        fas_geom_boxplot() +
        scale_x_continuous(name = glue('Neo-antigen depletion between\n',
                                       'PS = 0 and PS = 1')) +
        scale_y_discrete(name = glue('')) +
        scale_fill_manual(name = hl_lgd_name, values = man_colors) +
        theme_fas(legend.pos = 'none') +
        facet_wrap(~project_extended, labeller = label_parsed, scale = 'free')
    plot_basename <- 'boxplot'
    fn <- glue(fn)
    ggsave(filename = fn, plot = p, height = 12, width = 17.4, unit = 'cm')
    return(invisible())
  })
}
formals(IE_boxplot) <- formals(prep_IE_plot)


IE_lineplot <- function() {
  calling_args <- as.list(environment()) %>% replace('return_mode', 'overview')
  with(do.call(prep_IE_plot, calling_args), {
    if (null_dat(t_dat)) return(NULL)
    # if (T) {
    #   median_dat <- t_dat %>%
    #     group_by(project_extended, .data[[hl_var]]) %>%
    #     summarise(depletion_full = median(depletion_full, na.rm = T)) %>%
    #     dplyr::mutate(hl_var = factor(as.character(.data[[hl_var]]))) %>%
    #     { . }
    # } else {
    #   median_dat <- t_dat %>%
    #     group_by(project_extended, .data[[hl_var]]) %>%
    #     filter(depletion_full == median(depletion_full, na.rm = T))
    #     dplyr::mutate(hl_var = factor(as.character(.data[[hl_var]]))) %>%
    #     { . }
    # }

    p <- ggplot(t_dat, aes(x = hl_var, y = depletion_full,
                           color = color_var, group = group_var)) +
      geom_point(alpha = .1) +
      geom_line(alpha = .1) +
      # geom_point(data = median_dat,
      #            mapping = aes(y = depletion_full, x = x_var),
      #            alpha = .8, color = 'grey40', size = 2,
      #            inherit.aes = F) +
      fas_geom_boxplot(mapping = aes(x = hl_var, y = depletion_full,
                                     group = hl_var),
                       inherit.aes = F, fill = NA) +
      scale_y_continuous(name = glue('Neo-antigen depletion between\n',
                                     'PS = 0 and PS = 1')) +
      scale_x_discrete(name = glue('')) +
      scale_colour_discrete(name = glue('')) +
      theme_fas() +
      facet_wrap(~project_extended, labeller = label_parsed, scale = 'free') +
      gg_legend_alpha_cancel
      # coord_cartesian(ylim = c(0, 2))

    plot_basename <- 'lineplot'
    fn <- glue(fn)
    ggsave(filename = fn, plot = p, height = 18, width = 17.4, unit = 'cm')
    return(invisible())
  })
}
formals(IE_lineplot) <- formals(prep_IE_plot)


gen_depletion_tally_tables <- function() {
  calling_args <- as.list(environment()) %>% replace('return_mode', 'overview')
  with(do.call(prep_IE_plot, calling_args), {
    if (null_dat(t_dat)) return(NULL)
    ## Focus allele based tally. Are some alleles biased towards
    ## enrichment/depletion?
    allele_count <- t_dat %>%
      dplyr::mutate(depletion_binned = cut(depletion_full,
                                           c(-Inf, 1, 1+.Machine$double.eps, Inf),
                                           right = F)) %>%
      dplyr::group_by(project_extended, focus_allele, depletion_binned) %>%
      dplyr::summarise(N = n()) %>%
      tidyr::pivot_wider(id_cols = c(project_extended, focus_allele),
                         names_from = depletion_binned,
                         values_from = N,
                         values_fill = 0) %>%
      dplyr::rename(E = Range_1, D = Range_3) %>%
      dplyr::mutate(total = D + E, score = D - E)

    tabA <- allele_count %>%
      dplyr::ungroup() %>%
      dplyr::summarize(total = sum(total), score = sum(score)) %>%
      gen_table_grob()

    tabB <- allele_count %>%
      dplyr::group_by(focus_allele) %>%
      dplyr::summarize(total = sum(total), score = sum(score)) %>%
      dplyr::arrange(-score) %>%
      gen_table_grob()

    tabC <- allele_count %>%
      dplyr::rename('Tumor type' = project_extended) %>%
      dplyr::arrange(-score) %>%
      gen_table_grob()

    # library(ComplexHeatmap)

    # heightDetails(tabC)
    # grobHeight(tabC)

    # tabA$widths
    # grobWidth(tabA)
    plot_panel_layout(plots = list(tabA, tabB, tabC), nrow = 1,
                      height = .4 * nrow(allele_count) + 1,
                      widths = c(.2, .2, .6),
                      width = 12,
                      filename = glue(fn, plot_basename = 'IE_tables'))
  })
}
formals(gen_depletion_tally_tables) <- formals(prep_IE_plot)
