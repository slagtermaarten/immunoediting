#!/usr/bin/env Rscript

# devtools::install('~/libs/maartenutils')

if (!require('devtools'))
  try(install.packages('devtools'))
# if (!require('pacman'))
#   try(install.packages('pacman'))

hostname <- Sys.info()[["nodename"]]
servers <- c('paranoid', 'steroid', 'medoid', 'void', 'coley')
local_run <- !(hostname %in% servers)
maarten_libs <- ifelse(local_run, '/Users/maartenslagter/libs',
                       '/home/m.slagter/libs')
fas_libs <- ifelse(local_run, '/Users/maartenslagter/antigenic_space/libs',
                   '/home/m.slagter/antigenic_space/libs')
root_folder <- ifelse(local_run, '/Users/maartenslagter/antigenic_space',
                      '/DATA/research_projects/antigenic_space')
username <- Sys.info()['user']
pkg_dir <- .libPaths()[1]
message('Installing packages in ', pkg_dir)

if (!exists('force_installation'))
  force_installation <- F


#' Fully deteach package
#'
#' @param pkg \code{character} name of package to unload
detach_package <- function(pkg, character.only = FALSE) {
  if (!character.only) {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste('package', pkg, sep = ':')
  while(search_item %in% search()) {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}
# detach_package('fasanalysis')


#' Update installed version of developmental package
#'
#' @description Only install package if development version has undergone
#' changes since last install
#'
#' @param pkg_name name of package, assumed to equal directory containing it
#' @param libs_dir libs directory where devel versions are stored
#' @param devel_dir full path of development version of library
#' @param install_loc where to install package
#'
#' @return NULL invisibly
update_devel_pkg <- function(pkg_name = 'maartenutils',
  force_installation = F,
  libs_dir = '~/libs',
  devel_dir = file.path(libs_dir, pkg_name),
  install_loc = .libPaths()[1],
  installed_pkg_dir = list.files(install_loc, full.names = T, 
    pattern = pkg_name)) {

  dev_pkg_files <- 
    c(list.files(file.path(devel_dir, 'R'), recursive = T, full.names = T),
      list.files(file.path(devel_dir, 'src'), recursive = T, full.names = T),
      file.path(devel_dir, 'DESCRIPTION'), 
      file.path(devel_dir, 'NAMESPACE')) 
  latest_dev_time <- max(file.mtime(dev_pkg_files), na.rm = T)
  installed_pkg_files <-
    list.files(installed_pkg_dir, recursive = T, full.names = T)
  ## NA value will be returned if file does not exist in target directory
  latest_install_time <- tryCatch(max(file.mtime(installed_pkg_files)),
      error = function(e) { NA })
  if (is.na(latest_install_time) || latest_install_time < latest_dev_time || 
      force_installation) {
    devtools::install(devel_dir, lib = install_loc)
  }
  invisible()
}
# update_devel_pkg(pkg_name = 'genesets')


#' Remove locks in order to install packages
#'
#'
remove_locks <- function(pkg_dir = .libPaths()[1]) {
  locks <- list.files(pkg_dir, pattern = '00LOCK.*', full.names = T)
  if (length(locks) > 0) {
    print(locks)
    unlink(locks, recursive = T)
  }
  invisible()
}
remove_locks()

if (T) {
  ## Sometimes dependencies need to be manually and individually installed in
  ## order to not cause problems...
  pacman::p_load('gridExtra', lib = pkg_dir)
  pacman::p_load('gtable', lib = pkg_dir)
  remove_locks(pkg_dir = .libPaths()[1])
  # remove.pacakges('gdata')
  pacman::p_load('gdata', lib = pkg_dir)
  pacman::p_load('ggpubr')
  pacman::p_load('ggplot2', lib = pkg_dir)
  pacman::p_load('ggbeeswarm')
  pacman::p_load('tidyverse', lib = pkg_dir)
  pacman::p_load('munsell', lib = pkg_dir)
  pacman::p_load('RPostgreSQL', lib = pkg_dir)
  pacman::p_load('naturalsort', lib = pkg_dir)
  pacman::p_load('gridExtra', lib = pkg_dir)
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
  biocLite("edgeR")
  devtools::install_github('slowkow/ggrepel')
  install.packages('cowplot', lib = pkg_dir)
  # devtools::install_url('https://github.com/wilkelab/cowplot/archive/0.6.2.zip')
  # install.packages('operator.tools', lib = pkg_dir)
  # file.remove(file.path(pkg_dir, 'data.table', 'R', 'data.table.rdb'))
  install.packages('data.table')
  install.packages('Rcpp', lib = pkg_dir)
  # install.packages('Rcpp', type='source', lib = pkg_dir)
  # install.packages('git2r', type='source', lib = pkg_dir)
  # devtools::install_github('git2r', 'ropensci')
  # devtools::install_github('ropensci/git2r')
  install.packages('git2r')
  # install.packages('git2r', type = 'mac.binary.mavericks')
  # install.packages('data.table', repos = 'https://rdatatable.github.io/data.table')
}


if (F) {
  update.packages(ask = F, lib.loc = pkg_dir, instlib = pkg_dir)
}


if (T) {
  for (p in c('serializer', 'maartenutils', 'quickMHC')) {
    message(p)
    detach_package(p)
    update_devel_pkg(pkg_name = p, libs_dir = maarten_libs,
                     install_loc = pkg_dir,
                     force_installation = force_installation)
    detach_package(p)
  }
}

if (F) {
  fas_packages <- c('firehosedownload', 'fasanalysis')
  # fas_packages <- c('firehosedownload')
  for (p in fas_packages) {
    message(p)
    detach_package(p)
    update_devel_pkg(pkg_name = p, libs_dir = fas_libs, install_loc = pkg_dir,
                     force_installation = force_installation)
    detach_package(p)
  }
}

if (T) {
  for (p in c('genesets')) {
    message(p)
    detach_package(p)
    update_devel_pkg(pkg_name = p, libs_dir = maarten_libs,
                     install_loc = pkg_dir,
                     force_installation = force_installation)
    detach_package(p)
  }
}

if (F) {
  if (T) {
    devtools::load_all(file.path(root_folder, 'libs', 'maartenutils'))
    devtools::load_all(file.path(root_folder, 'libs', 'fasanalysis'))
  }

  if (exists("twoD_sens_analysis")) rm(twoD_sens_analysis)
  if (exists("donor_summary")) rm(donor_summary)
  if (exists("dtf")) rm(dtf)
  if (exists("ds")) rm(ds)
  if (exists("purity_f")) rm(purity_f)
}
