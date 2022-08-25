source(file.path(IE_root, 'load_optitype.R'))
library(RPostgreSQL)
library(glue)

PS_requires_computation <- maartenutils::gen_time_comparator(
  minimum_mod_time = '2020-07-09 11:12 CET',
  verbose = F)


kill_all_connections <- function() {
  pg <- DBI::dbDriver('PostgreSQL')
  all_cons <- dbListConnections(pg)
  for (con in all_cons) {
    dbDisconnect(con)
  }
  # dbDisconnect(all_cons[[1]])
  return(invisible())
}


#' Execute a GET query in an SQL database
#'
#'
gq <- function(query, con = NULL) {
  if (is.null(con)) {
    pg <- DBI::dbDriver('PostgreSQL')
    con <- RPostgreSQL::dbConnect(pg,
                                  user = 'm.slagter', 
                                  dbname = 'binding_affinity')
    on.exit(RPostgreSQL::dbDisconnect(con))
  }
  tryCatch(RPostgreSQL::dbGetQuery(con, glue(query)),
    error = function(e) {
      message(glue('Could not execute: {query}')); print(e)
  })
}


#' Store all peps presented by a particular allele in a SQL table to minimize
#' subsequent I/O with R
#'
#'
prep_focal_peps_table <- function(focus_allele, pr = 1.9, aff = NULL,
                                  debug = F, redo = F) {
  pg <- DBI::dbDriver('PostgreSQL')
  conn <- RPostgreSQL::dbConnect(pg, user = 'm.slagter',
                                 dbname = 'binding_affinity')
  on.exit(RPostgreSQL::dbDisconnect(conn))

  if (is.null(pr) && !is.null(aff)) {
    name_appendix <- glue('aff_{aff}')
    sel_statement <- glue('where "affinity" <= {aff}')
  } else if (!is.null(pr) && is.null(aff)) {
    name_appendix <- glue('percentile_rank_{pr}')
    sel_statement <- glue('where "percentile_rank" <= {pr}')
  } else {
    stop('Exactly one of pr and aff needs to be non-NULL')
  }

  table_name <- glue('all_peps_focal_{focus_allele}_{name_appendix}')

  if (RPostgreSQL::dbExistsTable(conn, table_name) && !redo) {
    focus_peps <- dbGetQuery(conn, sprintf('select count(*) from "%s"', 
                                           table_name)) %>% unlist()
    stopifnot(focus_peps > 1e5)
    return(invisible(NULL))
  }

  if (debug) {
    maartenutils::mail_notify(subject = 'prep_focal_peps_table',
                              msg = 'browser invoked')
    browser()
  }

  if (F) {
    ## Bulldoze away the previous attempt at this table
    tryCatch(RPostgreSQL::dbSendQuery(conn, 
                                      glue('drop table "{table_name}";')),
             error = function(e) { print(e) })
  }

  create_table <- glue('create table "{table_name}" ("peptide" character(9));')
  tryCatch(RPostgreSQL::dbSendQuery(conn, create_table),
           error = function(e) { print(e); browser() })

  select_test_peps <- glue('\\
    insert into "{table_name}" ( \\
        "peptide" \\
    ) \\
    select \\
        "peptide" \\
    from \\
        "all_peps_focal" \\
            natural inner join \\
        (select "peptide" from "{focus_allele}" \\
         {sel_statement}) as "{focus_allele}_focus_peps"')

  tryCatch(RPostgreSQL::dbSendQuery(conn, select_test_peps),
    error = function(e) { print(e); browser() })
  rm(pg)
  return(invisible())
}
# prep_focal_peps_table('A3301', pr = 1.9, debug = F, redo = F)


#' Generate SQL statement to compute PS (presentation score)
#'
#' @param focus_allele Peptides from this allele are used for PS
#' @param hla_alleles Alleles to check for presentation of focus allele peptides
#' @param method. If set to 'mean', compute the average amount of alleles with
#' which focus allele peptides are presented. If set to 'union', compute
#' the fraction of focus allele peptides to be presented by at least one HLA
#' allele
#'
gen_master_statement <- function(focus_allele = 'A0201',
                                 pr = 1.9,
                                 pr_ancillary = 1.9,
                                 aff = NULL,
                                 aff_ancillary = NULL,
                                 hla_alleles = c('A0201', 'B2705'),
                                 integration_method = 'mean') {
  hla_alleles <- unique(hla_alleles)
  if (is.null(pr) && !is.null(aff)) {
    bool_statement <- purrr::map(hla_alleles, function(h)
      glue('("{h}_aff" <= {aff_ancillary})')) %>%
      unlist %>%
      paste(collapse = ' or ')

    merge_prs <- purrr::imap(hla_alleles, function(h, i)
      glue('natural inner join \\
        (select "peptide", "affinity" as "{h}_aff" from "{h}") as "t{i}"')) %>%
      unlist %>%
      paste(collapse = ' ')

    res <- glue('select avg("pep_pres_bool") from (\n\\
        select cast({bool_statement} as integer) as "pep_pres_bool" from (\n\\
          (select "peptide" from "all_peps_focal_{focus_allele}") as "focus_peps" \n\\
            {merge_prs}
          ) as "coerced"
        ) as "computed";')
  } else if (!is.null(pr) && is.null(aff)) {
    ## This is an example of the SQL query we're creating here for a patient
    ## with HLA alleles A2501 and B3503 for reference allele A0201
    # select avg("pep_pres_bool") from (
    #   select cast(("A2501_pr" <= 1.9) or ("B5201_pr" <= 1.9) as integer) as 
    #   "pep_pres_bool" from (
    #     (select "peptide" from "all_peps_focal_A0201") as "focus_peps"
    #       natural inner join
    #     (select "peptide", "percentile_rank" as "A2501_pr" from "A2501") 
    #     as "t1"
    #       natural inner join
    #     (select "peptide", "percentile_rank" as "B5201_pr" from "B5201") 
    #     as "t2"
    #   ) as "coerced"
    # ) as "computed";
    # select avg("pep_pres_bool") from (
    #  select cast(("A1101_pr" <= 1.9) or ("B0702_pr" <= 1.9) as integer) 
    #  as "pep_pres_bool", "A1101_pr", "B0702_pr" from (
    #    (select "peptide" from "all_peps_focal_A0201") as "focus_peps" 
    #      natural inner join
    #    (select "peptide", "percentile_rank" as "A1101_pr" from "A1101") 
    #    as "t1"
    #      natural inner join
    #    (select "peptide", "percentile_rank" as "B0702_pr" from "B0702") 
    #    as "t2"
    #  ) as "coerced"
    # ) as "computed";
    bool_statement <- purrr::map(hla_alleles, function(h)
      glue('("{h}_pr" <= {pr_ancillary})')) %>%
      unlist %>%
      paste(collapse = ' or ')

    merge_prs <- purrr::imap(hla_alleles, function(h, i)
      glue('natural inner join \\
        (select "peptide", "percentile_rank" as "{h}_pr" from "{h}") as "t{i}"')) %>%
      unlist %>%
      paste(collapse = ' ')

    res <- glue('select avg("pep_pres_bool") from (\n\\
        select cast({bool_statement} as integer) as "pep_pres_bool" from (\n\\
          (select "peptide" from "all_peps_focal_{focus_allele}") as "focus_peps" \n\\
            {merge_prs}
          ) as "coerced"
        ) as "computed";')
  } else {
    stop('Exactly one of pr and aff needs to be non-NULL')
  }
  return(res)
}


#' Master function to compute the presentation score
#'
#' @param integration_method If 'mean' or 'union', compute the overlap
#' between entire repertoire and focus allele, if 'sum', sum the pairwise
#' overlaps between individual alleles and the focus allele
#'
compute_presentation_score <- function(focus_allele, donor_hlas,
                                       pr = 1.9, pr_ancillary = 1.9,
                                       aff = NULL, aff_ancillary = NULL,
                                       hla_sim_method = '',
                                       redo = F,
                                       limit_hla_sim = F,
                                       rds_sub_dir = file.path(rds_dir,
                                         'pres_scores'),
                                       check_cached = F,
                                       integration_method = 'union',
                                       omit_mean_score = T,
                                       verbose = F,
                                       ncores = 1) {
  fn_sel <- gen_affinity_cutoff_flags(pr = pr, pr_ancillary = pr_ancillary,
    aff_ancillary = aff_ancillary, aff = aff)
  donor_hlas_f <- paste(donor_hlas, collapse = '-')
  res_fn <- file.path(rds_sub_dir, glue('{focus_allele}_{donor_hlas_f}\\
  -{fn_sel}\\
  -hla_sim_method_{hla_sim_method}\\
  -integration_method_{integration_method}.rds'))
  
  ## Determine (look-up) filename to write results to
  if (!redo && file.exists(res_fn) && 
      !check_cached && !PS_requires_computation(res_fn)) {
    if (verbose) mymessage(glue('Reading from {res_fn} - {file.mtime(res_fn)}'))
    res <- readRDS(res_fn)
    return(res)
  } else {
    compute_date <- Sys.time()
  }

  if (length(donor_hlas) == 0) {
    res <- list('hla' = focus_allele, 'donor_hlas' = NA, 'ps'= 0)
    return(res)
  }

  ## Not needed if this function is called from compute_PS_overview, 
  ## already one in compute_PS_overview
  prep_focal_peps_table(focus_allele = focus_allele, pr = pr, aff = aff)

  if (integration_method == 'union' && focus_allele %in% donor_hlas) {
    v <- 1
  } else {
    sql_query <- gen_master_statement(focus_allele = focus_allele,
                                      pr = pr,
                                      pr_ancillary = pr_ancillary,
                                      aff = aff,
                                      aff_ancillary = aff_ancillary,
                                      hla_alleles = donor_hlas,
                                      integration_method =
                                        integration_method)
    v <- gq(sql_query) %>% unlist %>% unname
  }

  res <- list('hla' = focus_allele, 'donor_hlas' = donor_hlas, 'ps'= v)

  if (check_cached && file.exists(res_fn)) {
    res_stored <- readRDS(res_fn)
    if (compare_nums(res$ps, res_stored$ps) == F) {
      res_fn <- maartenutils::append_date_to_fn(res_fn)
      mywarning(sprintf('Detected mismatch with cached result: %s', res_fn))
    }
  }

  if (!dir.exists(rds_sub_dir)) {
    dir.create(rds_sub_dir, showWarnings = F)
  }
  if (verbose) { mymessage(glue('Writing to {res_fn}')) }
  saveRDS(res, res_fn)
  return(res)
}


identify_lost_alleles <- function(donor_id, min_coverage = 0, 
  max_lohhla_cov = Inf, debug_func = F) {
  w_readRDS('lohhla_overview', force_reload = F)
  l_min_coverage <- min_coverage
  l_donor_id <- donor_id
  ## Select relevant subset
  subs <- lohhla_overview %>%
    .[min_coverage == l_min_coverage & donor_id == l_donor_id]
  if (debug_func) {
    print(subs)
  }
  if (null_dat(subs) ||
      subs[, any(grepl('no_filtered_bam_present', message))]) {
    return(NULL)
  }
  subs_analysis <- subs[
    grepl('analysis_completed', message) &
    hla_allele_1_cov <= max_lohhla_cov &
    hla_allele_2_cov <= max_lohhla_cov]
  LOHHLA_complete <- subs_analysis[, .N] == 3
  LOHHLA_AB_complete <-
    subs_analysis[grepl('_(a|b)_', hla_allele_1), .N] == 2

  loss_indicators <- 'no_pileup_for_tumor'
  escape_indicators <- 'normal|homozygous'
  # grepl(escape_indicators, 'tumor;normal;homozygous_alleles_not_implemented')
  # grepl(loss_indicators, 'tumor;normal;homozygous_alleles_not_implemented')
  lost_alleles <- c()
  for (i in 1:nrow(subs)) {
    messages <- strsplit(subs[i, message], ';')[[1]]
    for (message in messages) {
      if (grepl(loss_indicators, message) &&
          !grepl(escape_indicators, message)) {
        gene_lost <- gsub('.*_(\\w)_\\d', '\\1', message)
        allele_lost <- as.integer(gsub('.*_(\\d)', '\\1', message))
        lost_alleles <- c(lost_alleles,
          subs[i, get(glue('hla_allele_{allele_lost}'))])
      }
    }
    if (subs[i, !is.na(hla_allele_1_cn_ci_upper)] &&
        subs[i, hla_allele_1_cn_ci_upper < 0] &&
        subs[i, !is.na(hla_allele_1_cov)] &&
        subs[i, hla_allele_1_cov <= max_lohhla_cov]) {
      lost_alleles <- c(lost_alleles, subs[i, hla_allele_1])
    }
    if (subs[i, !is.na(hla_allele_2_cn_ci_upper)] &&
        subs[i, hla_allele_2_cn_ci_upper < 0] &&
        subs[i, !is.na(hla_allele_2_cov)] &&
        subs[i, hla_allele_2_cov <= max_lohhla_cov]) {
      lost_alleles <- c(lost_alleles, subs[i, hla_allele_2])
    }
  }
  lost_alleles <- unique(quickMHC::shortenHLA(lost_alleles))
  return(list('alleles' = lost_alleles,
              'LOHHLA_complete' = LOHHLA_complete,
              'LOHHLA_AB_complete' = LOHHLA_AB_complete))
}


gen_affinity_cutoff_flags <- function(
  aff = NULL, aff_ancillary = NULL, pr = 1.9, pr_ancillary = 1.9) {

  if (is.null(pr) && !is.null(aff)) {
    fn_sel <- glue('aff={aff}-aff_ancillary={aff_ancillary}')
  } else if (!is.null(pr) && is.null(aff)) {
    fn_sel <- glue('pr={pr}-pr_anc={pr_ancillary}')
  } else {
    stop('Exactly one of pr and affinity needs to be non-NULL')
  }
  ## if (fn_sel == 'pr=1.9-pr_anc=1.9') fn_sel <- ''
  ## 2021-02-10 11:03 I should have used = at the time, but don't want to
  ## recompute or change file names of already computed files
  fn_sel <- gsub('=', '_', fn_sel)
  return(fn_sel)
}


gen_PS_overview_fn <- function(hla_sim_method = 'cutoff',
                               integration_method = 'union',
                               hla_alleles = 'A0201',
                               pr = NULL, pr_ancillary = NULL,
                               aff = NULL, aff_ancillary = NULL) {
  fn_sel <- gen_affinity_cutoff_flags(pr = pr, pr_ancillary = pr_ancillary,
    aff_ancillary = aff_ancillary, aff = aff)
  hla_flag <- paste0('-hla_alleles=', paste(hla_alleles, collapse = "_"))
  file.path(rds_dir, glue::glue('focus_allele_avg_HLA_presentation\\
      {hla_flag}\\
      {make_flag(hla_sim_method)}\\
      -{fn_sel}\\
      {make_flag(integration_method)}.rds')
    )
}


#' Compute presentation scores for all OptiType'd TCGA patients
#'
#' @param hla_alleles 'Focus' alleles to compute presentation scores for
#'
#'
compute_PS_overview <- function(hla_alleles = focus_hlas,
                                hla_sim_method = 'cutoff',
                                integration_method = 'union',
                                pr = 1.9, pr_ancillary = 1.9,
                                # pr = NULL, pr_ancillary = NULL,
                                # aff = 100, aff_ancillary = 100,
                                aff = NULL, aff_ancillary = NULL,
                                mail_notify = T,
                                ncores = 1,
                                verbose = T,
                                min_coverage = 0,
                                redo = F) {
  fn <- gen_PS_overview_fn(integration_method = integration_method,
                           hla_sim_method = hla_sim_method,
                           hla_alleles = hla_alleles,
                           pr = pr, pr_ancillary = pr_ancillary,
                           aff = aff, aff_ancillary = aff_ancillary)
  if (!PS_requires_computation(fn) && !redo) {
    return(readRDS(fn))
  }

  ## Ensure focal pep tables exist before we start parallelizing
  mymessage(paste0('Prepping all peptide predictions for ',
                   paste(hla_alleles, collapse = ', ')))
  for (hla_allele in hla_alleles) {
    prep_focal_peps_table(hla_allele, pr = pr, aff = aff, debug = F)
  }

  mymessage(paste0('Computing PS_overview for ',
                   paste(hla_alleles, collapse = ', ')))

  param_table_arg <- list(
    'focus_allele' = hla_alleles,
    'donor_id' = hla_class_dat[, donor_id],
    'aff' = aff,
    'aff_ancillary' = aff_ancillary,
    'pr_ancillary' = pr_ancillary,
    LOH_HLA = c('LOHHLA', 'no_LOHHLA'),
    ## 2019-08-14 13:01 The 'strict_LOHHLA' will result in identical scores but
    ## will merely require some patients to be excluded, we thus don't need to
    ## recompute scores but can suffice by flagging whether all HLA alleles were
    ## sufficiently covered or not
    'pr' = pr,
    stringsAsFactors = F
  )

  if (is.null(pr) && !is.null(aff)) {
    param_table_arg <- 
      param_table_arg[setdiff(names(param_table_arg), 
                              grep('pr', names(param_table_arg), value = T))]
  } else if (!is.null(pr) && is.null(aff)) {
    param_table_arg <- 
      param_table_arg[setdiff(names(param_table_arg), 
                              grep('aff', names(param_table_arg), value = T))]
  } else {
    stop('Exactly one of pr and affinity needs to be non-NULL')
  }
  param_table <- do.call(expand.grid, param_table_arg)

  kill_all_connections()

  doParallel::registerDoParallel(cores = ncores)
  ROs <- plyr::llply(1:nrow(param_table), function(i) {
    if (verbose) mymessage(sprintf('batch %d/%d', i, nrow(param_table)))
    donor_id <- as.character(param_table[i, 'donor_id'])
    focus_allele <- as.character(param_table[i, 'focus_allele'])
    pr <- as.numeric(param_table[i, 'pr'])
    pr_ancillary <- as.numeric(param_table[i, 'pr_ancillary'])
    aff <- as.numeric(param_table[i, 'aff'])
    aff_ancillary <- as.numeric(param_table[i, 'aff_ancillary'])
    HLA_LOH <- as.character(param_table[i, 'LOH_HLA'])

    donor_AB_hlas <-
      hla_class_dat[donor_id, c(hla_a1, hla_a2, hla_b1, hla_b2)]
    donor_C_hlas <-
      hla_class_dat[donor_id, c(hla_c1, hla_c2)]

    if (HLA_LOH == 'LOHHLA') {
      lost_alleles <- identify_lost_alleles(donor_id,
                                            min_coverage = min_coverage)
      donor_AB_hlas <- setdiff(donor_AB_hlas, lost_alleles$alleles)
      donor_C_hlas <- setdiff(donor_C_hlas, lost_alleles$alleles)
    }
    donor_hlas <- c(donor_AB_hlas, donor_C_hlas)

    scores <- list('mean_score' = donor_hlas, 
                   'mean_score_AB' = donor_AB_hlas) %>%
      purrr::map(function(l_hla_repertoire) {
        res <- compute_presentation_score(focus_allele = focus_allele,
                                          donor_hlas = l_hla_repertoire,
                                          redo = F,
                                          hla_sim_method = hla_sim_method,
                                          check_cached = F,
                                          integration_method = integration_method,
                                          verbose = verbose,
                                          pr = pr,
                                          pr_ancillary = pr_ancillary)
      })

    res <- list('donor_id' = donor_id,
                'hla' = focus_allele,
                'mean_score' = scores[['mean_score']][['ps']],
                'mean_score_AB' = scores[['mean_score_AB']][['ps']],
                'LOH_HLA' = HLA_LOH,
                'N_AB_alleles' = length(donor_AB_hlas),
                'N_C_alleles' = length(donor_C_hlas),
                'hla_sim_method' = hla_sim_method)

    if (HLA_LOH == 'LOHHLA') {
      res[['LOHHLA_complete']] <- lost_alleles$LOHHLA_complete
      res[['LOHHLA_AB_complete']] <- lost_alleles$LOHHLA_AB_complete
    }
    return(res)
  }, .parallel = (ncores > 1), .progress = 'text') %>% rbindlist(fill = T)

  setkey(ROs, donor_id, hla)
  ROs <- maartenutils::controlled_merge(ROs, hla_class_dat,
                                        by_cols = 'donor_id')
  ROs[, 'anchor_present' :=
      hla %in% c(hla_a1, hla_a2, hla_b1, hla_b2, hla_c1, hla_c2),
      by = 1:nrow(ROs)]
  ROs <- clean_columns('', ROs, 'hla_allele_status')
  attr(ROs, 'ro_fn') <- fn
  saveRDS(ROs, fn)
  if (mail_notify) {
    mail_notify(subject = 'Presentation score overview completed',
                msg = fn)
  }
  return(invisible(ROs))
}


#' Compute the binding overlap between two HLA alleles
#'
#' Used in Fig 3A of ms
#'
compute_pairwise_corroboration <- function(
  focus_allele = 'A0201',
  ancillary_allele = 'B2705',
  pr = 1.9, pr_ancillary = 1.9,
  hla_sim_method = 'cutoff',
  limit_hla_sim = F, redo = F,
  rds_sub_dir = file.path(rds_dir, 'focus_allele_avg_HLA_presentation')) {

  focus_allele <- quickMHC::shortenHLA(focus_allele)
  ancillary_allele <- quickMHC::shortenHLA(ancillary_allele)
  if (focus_allele == ancillary_allele) return(1)
  hla_sim_method <- match.arg(hla_sim_method,
    choices = c('', 'cutoff', 'cosine', 'scalar_projection',
                'ols_coef', 'trimmed_ols_coef', 'ols_coef_pr'))
  if (hla_sim_method == 'cutoff')
    hla_sim_method <- ''
  method_string <- prepend_hyphen(hla_sim_method)

  l_fn <- file.path(rds_sub_dir,
    glue::glue('pairwise_corroboration-{focus_allele}-{ancillary_allele}',
         '-{pr}-{pr_ancillary}{method_string}.rds'))

  if (length(l_fn) > 1) browser()
  if (file.exists(l_fn) && !redo) {
    message(glue::glue('Returning precomputed {basename(l_fn)}'))
    if (limit_hla_sim) {
      return(min(1, readRDS(l_fn)))
    } else {
      return(readRDS(l_fn))
    }
  }

  pg <- DBI::dbDriver('PostgreSQL')
  conn <- RPostgreSQL::dbConnect(pg, user = 'm.slagter',
                                 dbname = 'binding_affinity')
  on.exit(RPostgreSQL::dbDisconnect(conn))

  prep_focal_peps_table(focus_allele, pr = pr)

  ## Perform predictions for all alleles for which these aren't available yet
  for (hla_allele in c(focus_allele, ancillary_allele)) {
    if (!RPostgreSQL::dbExistsTable(conn, hla_allele)) {
      message(glue::glue('Creating table {hla_allele}'))
      quickMHC::precompute_peps(peptides = all_peps, hla_allele = hla_allele,
        perform_STS = F, percentile_rank = pr, ncores = 10, batch_size = 1e6,
        index = T)
    }
  }
  h <- ancillary_allele

  if (hla_sim_method %in% c('', 'cutoff')) {
    ## This is an example of the SQL we're creating here for ancillary allele
    ## A2501 and reference allele A0201
    ## select avg("psum") from (
    ##   select cast ("A2501_pr" <= 1.9 as integer) as "psum" from (
    ##     (select "peptide" from "all_peps_focal_A0201") as "focus_peps"
    ##       natural inner join
    ##     (select "peptide", "percentile_rank" as "A2501_pr" from "A2501") as "t1"
    ##   )
    ##   as "merged")
    ## as "computed";
    conversion_statement <- glue::glue('cast ("{h}_pr" <= {pr_ancillary} as integer)')

    merge_prs <- glue::glue('natural inner join \\
        (select "peptide", "percentile_rank" as "{h}_pr" from "{h}") as "t1"')

    combined_sql <- glue::glue('select avg("psum") from (
        select {conversion_statement} as "psum" from (
          (select "peptide" from "all_peps_focal_{focus_allele}") as "focus_peps"
            {merge_prs}
          )
        as "merged")
      as "computed";')

    res <- as.numeric(unlist(gq(combined_sql, con = conn)))
  } else if (hla_sim_method == 'cosine') {
    sql <- glue::glue('select "affinity", "ancillary_affinity" from (
      "all_peps_focal"
        natural inner join
      (select "peptide", "affinity" from "{focus_allele}") as "focus_peps"
        natural inner join
      (select "peptide", "affinity" as "ancillary_affinity" from "{h}")
        as "merged") where "affinity" <= 500;')
    sql_table <- gq(sql, con = conn)
    res <- cosine_dist(sql_table$affinity, sql_table$ancillary_affinity)
  } else if (hla_sim_method == 'scalar_projection') {
    ## Chosen as such:
    ## select min(peptide_score_log50k) from "A0201" where "percentile_rank" = 1.9;
    min_peptide_score <- .47922
    sql <- glue('select "peptide_score_log50k", "ancillary_affinity" from (
      "all_peps_focal"
        natural inner join
      (select "peptide", "peptide_score_log50k" from "{focus_allele}") as "focus_peps"
        natural inner join
      (select "peptide", "peptide_score_log50k" as "ancillary_affinity" from "{h}")
        as "merged") where "peptide_score_log50k" >= {min_peptide_score};')

    sql_table <- gq(sql, con = conn)
    ## Normalize to 500 nM in order to make the absolute scores for all alleles
    ## comparable (the direction is what truly matters here)
    # sql_table$peptide_score_log50k <- sql_table$peptide_score_log50k *
    #   (500 / max(sql_table$affinity, na.rm = T))
    # ## Inverse the affinities in order to make small values more influential and
    # ## vice versa
    # sql_table$affinity <- 1 / sql_table$affinity
    # sql_table$ancillary_affinity <- 1 / sql_table$ancillary_affinity
    res <- scalar_projection(
      sql_table$peptide_score_log50k, sql_table$ancillary_affinity)
  } else if (hla_sim_method == 'ols_coef') {
    sql <- glue('select "affinity", "ancillary_affinity" from (
      "all_peps_focal"
        natural inner join
      (select "peptide", "affinity" from "{focus_allele}") as "focus_peps"
        natural inner join
      (select "peptide", "affinity" as "ancillary_affinity" from "{h}")
        as "merged") where "affinity" <= 255;')

    sql_table <- gq(sql, con = conn)
    res <- ols_coef(sql_table$ancillary_affinity, sql_table$affinity)
  } else if (hla_sim_method == 'ols_coef_pr') {
    sql <- glue('select "percentile_rank", "a_percentile_rank" from (
      "all_peps_focal"
        natural inner join
      (select "peptide", "percentile_rank" from "{focus_allele}") as "focus_peps"
        natural inner join
      (select "peptide", "percentile_rank" as "a_percentile_rank" from "{h}")
        as "merged") where "percentile_rank" <= {pr};')

    sql_table <- gq(sql, con = conn)
    # res <- ols_coef(sql_table$percentile_rank, sql_table$a_percentile_rank)
    res <- ols_coef(sql_table$a_percentile_rank, sql_table$percentile_rank)
  } else if (hla_sim_method == 'trimmed_ols_coef') {
    sql <- glue('select "peptide", "affinity", "ancillary_affinity" from (
      "all_peps_focal"
        natural inner join
      (select "peptide", "affinity" from "{focus_allele}") as "focus_peps"
        natural inner join
      (select "peptide", "affinity" as "ancillary_affinity" from "{h}")
        as "merged") where "affinity" <= 255;')

    sql_table <- gq(sql, con = conn)
    sql_table <- trim_table(dtf = sql_table, cn = 'ancillary_affinity', .1)
    res <- ols_coef(sql_table$ancillary_affinity, sql_table$affinity)
    sql_table$predicted <- sql_table$ancillary_affinity * res
    sql_table$residuals <- sql_table$predicted - sql_table$affinity
    sql_table$residuals_n <- sql_table$residuals / sql_table$affinity
    setDT(sql_table)
    print(sql_table[, .(focus_allele, ancillary_allele, .N),
      keyby = cut(residuals_n, breaks = 10)] %>%
      .[, 'perc' := 100 * N / sum(N)])
    # if (res >= 5 && ncores == 1) browser()
  }
  saveRDS(res, l_fn)
  message(glue('Writing {basename(l_fn)} -> {res}'))

  if (limit_hla_sim) {
    return(min(1, readRDS(l_fn)))
  } else {
    return(res)
  }
}


pairwise_HLA_groups <- function(dtf, 
                                focus_allele = attr(dtf, 'ds_opts')$hla_allele, 
                                ncores = 32) {
  doParallel::registerDoParallel(cores = ncores)
  source('~/antigenic_space/maarten-analyses/immune_editing/load_optitype.R')

  CS <- plyr::llply(all_hlas, function(x) {
    compute_pairwise_corroboration(focus_allele = focus_allele, 
                                   ancillary_allele = x) > .2
  }, .parallel = (ncores > 1)) %>% unlist

  hla_levs <- extract_hla_levs(dtf, 'hla_allele_status')

  patient_categories <- plyr::adply(optitype_tcga, 1, function(r) {
    if (any(unlist(r) == gsub('\\*', '', ppHLA(focus_allele)))) {
      return(data.frame('hla_allele_status_b' = hla_levs[['test']]))
    } else if (any(unlist(r) %in% names(which(CS)))) {
      return(data.frame('hla_allele_status_b' = hla_levs[['like']]))
    } else {
      return(data.frame('hla_allele_status_b' = hla_levs[['ctrl']]))
    }
  }, .parallel = (ncores > 1))

  patient_categories$hla_allele_status_b <- 
    factor(patient_categories$hla_allele_status_b, levels = hla_levs)

  if ('hla_allele_status_b' %in% colnames(dtf)) {
    dtf[, hla_allele_status_b := NULL]
  }

  dtf <- controlled_merge(dtf, 
                          patient_categories[, .(donor_id, hla_allele_status_b)])
  return(dtf)
}
