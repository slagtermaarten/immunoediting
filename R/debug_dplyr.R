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
library(testthat)

## Newer version of ComplexHeatmap does not yield problems
library(ComplexHeatmap)
test_dplyr()

## Newer version of ComplexHeatmap does not yield problems
devtools::load_all('~/libs/ComplexHeatmap')
test_dplyr()
## The session is ruined, restarting R is the only way I know of to reset it
library(ComplexHeatmap)
test_dplyr()
