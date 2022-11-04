#' ---
#' output:
#'   md_document:
#'     pandoc_args:
#'       - '--from=markdown-implicit_figures'
#'       - '--to=commonmark'
#'       - '--wrap=preserve'
#' ---



#+ reprex-setup, include = FALSE
options(tidyverse.quiet = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = TRUE)
knitr::opts_knit$set(upload.fun = knitr::imgur_upload)

#+ reprex-body
# install.packages('reprex')
library(reprex)

test_dplyr <- function() {
  library(dplyr)
  res <- tryCatch({
    by_cyl <- mtcars %>% dplyr::group_by(cyl)
    by_cyl %>% dplyr::filter(disp == max(disp))
  }, error = function(e) NULL)
  !is.null(res)
}
test_dplyr()

## Newer version of ComplexHeatmap does not yield problems
library(ComplexHeatmap)
test_dplyr()

## Newer version of ComplexHeatmap does not yield problems
devtools::load_all('~/libs/ComplexHeatmap')
test_dplyr()
## The session is ruined, restarting R is the only way I know of to reset it
library(ComplexHeatmap)
test_dplyr()



#' <sup>Created on `r Sys.Date()` by the [reprex package](https://reprex.tidyverse.org) (v`r utils::packageVersion("reprex")`)</sup>


