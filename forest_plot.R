unlockBinding("forestplot.default", as.environment("package:forestplot"))
lockBinding("forestplot.default", as.environment("package:forestplot"))
assign('forestplot.default', 
       envir = as.environment('package:forestplot'), 
       value = function(
    labeltext, mean, 
    lower, upper, align, 
    is.summary = FALSE,
    graph.pos = "right", hrzl_lines, clip = c(-Inf, Inf), xlab = "",
    zero = ifelse(xlog, 1, 0), graphwidth = "auto", colgap, lineheight = "auto",
    line.margin, col = fpColors(), txt_gp = fpTxtGp(), xlog = FALSE,
    xticks, xticks.digits = 2, grid = FALSE, lwd.xaxis, lwd.zero,
    lwd.ci, lty.ci = 1, ci.vertices, ci.vertices.height = 0.1,
    boxsize, mar = unit(rep(5, times = 4), "mm"), title, legend,
    legend_args = fpLegend(), 
    new_page = getOption("forestplot_new_page", TRUE), 
    fn.ci_norm = fpDrawNormalCI, fn.ci_sum = fpDrawSummaryCI,
    fn.legend, shapes_gp = fpShapesGp(), ...)
{
  if (missing(colgap)) {
      colgap <- convertUnit(unit(6, "mm"), "npc", valueOnly = TRUE)
      if (colgap < 0.1)
          colgap <- unit(0.05, "npc")
      else colgap <- unit(colgap, "npc")
  }
  else if (!grid::is.unit(colgap)) {
      colgap <- as.numeric(colgap)
      if (is.na(colgap))
          stop("Invalid colgap argument")
  }
  colgap <- convertUnit(colgap, "mm")
  assert_class(txt_gp, "fpTxtGp")
  assert_class(col, "fpColors")
  if (missing(lower) && missing(upper) && missing(mean)) {
      if (missing(labeltext))
          stop("You need to provide the labeltext or", " the mean/lower/upper arguments")
      mean <- labeltext
      labeltext <- rownames(mean)
  }
  if (missing(lower) && missing(upper)) {
      assert(check_matrix(mean, ncols = 3), check_array(mean,
          d = 3), check_integer(dim(mean)[2], lower = 3, upper = 3))
  }
  assert_vector(zero, max.len = 2)
  if (missing(labeltext))
      labeltext <- rownames(mean)
  if (is.null(labeltext))
      stop("You must provide labeltext either in the direct form as an argument",
          " or as rownames for the mean argument.")
  if (missing(lower) && missing(upper)) {
      if (NCOL(mean) != 3)
          stop("If you do not provide lower/upper arguments your mean needs to have 3 columns")
      all <- prFpConvertMultidimArray(mean)
      mean <- all$mean
      lower <- all$lower
      upper <- all$upper
  }
  if (NCOL(mean) != NCOL(lower) || NCOL(lower) != NCOL(upper) ||
      NCOL(mean) == 0)
      stop("Mean, lower and upper contain invalid number of columns",
          " Mean columns:", ncol(mean), " Lower bound columns:",
          ncol(lower), " Upper bound columns:", ncol(upper))
  if (NCOL(mean) != length(col$box)) {
      col$box <- rep(col$box, length.out = NCOL(mean))
      col$line <- rep(col$lines, length.out = NCOL(mean))
  }
  if (!missing(legend)) {
      fn.legend <- prFpPrepareLegendMarker(fn.legend = fn.legend,
          col_no = NCOL(mean), row_no = NROW(mean), fn.ci_norm = fn.ci_norm)
  }
  if (!is.unit(lineheight) && !lineheight %in% c("auto", "lines"))
      stop("The argument lineheight must either be of type unit or set to 'auto',",
          " you have provided a '", class(lineheight), "' class")
  if (!missing(legend)) {
      if (length(legend) != ncol(mean))
          stop("If you want a legend you need to provide the same number of",
              " legend descriptors as you have boxes per line, currently you have ",
              ncol(mean), " boxes and ", length(legend), " legends.")
      if (is.list(legend_args$pos)) {
          legend_args$pos <- prFpGetLegendBoxPosition(legend_args$pos)
      }
      else if (!legend_args$pos %in% c("top", "right")) {
          stop("The legend is either a list positioning it inside the main plot or at the 'top' or 'right' side,",
              " the position '", legend_args$pos, "' is not valid.")
      }
      if (inherits(legend_args$gp, "gpar")) {
          if (!"col" %in% names(legend_args$gp)) {
              if (any(c("lwd", "lwd") %in% names(legend_args$gp))) {
                legend_args$gp[["col"]] = "black"
              }
              else {
                legend_args$gp[["col"]] = NA
              }
          }
      }
  }
  if (is.data.frame(mean))
      mean <- as.matrix(mean)
  if (is.data.frame(lower))
      lower <- as.matrix(lower)
  if (is.data.frame(upper))
      upper <- as.matrix(upper)
  if (new_page || dev.cur() == 1)
      grid.newpage()
  if (xlog) {
      if (any(mean < 0, na.rm = TRUE) || any(lower < 0, na.rm = TRUE) ||
          any(upper < 0, na.rm = TRUE) || (!is.na(zero) &&
          zero <= 0) || (!missing(clip) && any(clip <= 0, na.rm = TRUE)) ||
          (!missing(grid) && any(grid <= 0, na.rm = TRUE))) {
          stop("All argument values (mean, lower, upper, zero, grid and clip)",
              " should be provided as exponentials when using the log scale.",
              " This is an intentional break with the original forestplot function in order",
              " to simplify other arguments such as ticks, clips, and more.")
      }
      org_mean <- log(mean)
      org_lower <- log(lower)
      org_upper <- log(upper)
  }
  else {
      org_mean <- mean
      org_lower <- lower
      org_upper <- upper
  }
  if (NCOL(mean) > 1) {
      mean <- as.vector(mean)
      lower <- as.vector(lower)
      upper <- as.vector(upper)
  }
  nr <- NROW(org_mean)
  if (is.expression(labeltext)) {
      widthcolumn <- c(TRUE)
      nc <- 1
      label_type = "expression"
      label_nr <- length(labeltext)
  }
  else if (is.list(labeltext)) {
      if (all(sapply(labeltext, function(x) {
          length(x) == 1 && !is.list(x)
      }))) {
          labeltext <- list(labeltext)
      }
      if (!prFpValidateLabelList(labeltext))
          stop("Invalid labellist, it has to be formed as a matrix m x n elements")
      nc <- length(labeltext)
      widthcolumn = c()
      for (col.no in seq(along = labeltext)) {
          empty_row <- TRUE
          for (row.no in seq(along = labeltext[[col.no]])) {
              if (is.expression(labeltext[[col.no]][[row.no]]) ||
                !is.na(labeltext[[col.no]][[row.no]])) {
                empty_row <- FALSE
                break
              }
          }
          widthcolumn <- append(widthcolumn, empty_row)
      }
      label_type = "list"
      label_nr <- length(labeltext[[1]])
  }
  else if (is.vector(labeltext)) {
      widthcolumn <- c(FALSE)
      nc = 1
      labeltext <- matrix(labeltext, ncol = 1)
      label_type = "matrix"
      label_nr <- NROW(labeltext)
  }
  else {
      widthcolumn <- !apply(is.na(labeltext), 1, any)
      nc <- NCOL(labeltext)
      label_type = "matrix"
      label_nr <- NROW(labeltext)
  }
  if (nr != label_nr) {
      stop("You have provided ", nr, " rows in your", " mean arguement while the labels have ",
          label_nr, " rows")
  }
  if (is.character(graph.pos)) {
      graph.pos <- switch(graph.pos, right = nc + 1, last = nc +
          1, left = 1, first = 1, stop("The graph.pos argument has an invalid text argument.",
          " The only values accepted are 'left'/'right' or 'first'/'last'.",
          " You have provided the value '", graph.pos, "'"))
  }
  else if (is.numeric(graph.pos)) {
      if (!graph.pos %in% 1:(nc + 1))
          stop("The graph position must be between 1 and ",
              (nc + 1), ".", " You have provided the value '",
              graph.pos, "'.")
  }
  else {
      stop("The graph pos must either be a string consisting of 'left'/'right' (alt. 'first'/'last')",
          ", or an integer value between 1 and ", (nc + 1))
  }
  if (missing(align)) {
      if (graph.pos == 1)
          align <- rep("l", nc)
      else if (graph.pos == nc + 1)
          align <- c("l", rep("r", nc - 1))
      else align <- c("l", rep("c", nc - 1))
  }
  else {
      align <- rep(align, length.out = nc)
  }
  is.summary <- rep(is.summary, length = nr)
  if (is.matrix(mean)) {
      missing_rows <- apply(mean, 2, function(row) all(is.na(row)))
  }
  else {
      missing_rows <- sapply(mean, is.na)
  }
  fn.ci_norm <- prFpGetConfintFnList(fn = fn.ci_norm, no_rows = NROW(org_mean),
      no_cols = NCOL(org_mean), missing_rows = missing_rows,
      is.summary = is.summary, summary = FALSE)
  fn.ci_sum <- prFpGetConfintFnList(fn = fn.ci_sum, no_rows = NROW(org_mean),
      no_cols = NCOL(org_mean), missing_rows = missing_rows,
      is.summary = is.summary, summary = TRUE)
  lty.ci <- prPopulateList(lty.ci, no_rows = NROW(org_mean),
      no_cols = NCOL(org_mean))
  hrzl_lines <- prFpGetLines(hrzl_lines = hrzl_lines, is.summary = is.summary,
      total_columns = nc + 1, col = col, shapes_gp = shapes_gp)
  labels <- prFpGetLabels(label_type = label_type, labeltext = labeltext,
      align = align, nc = nc, nr = nr, is.summary = is.summary,
      txt_gp = txt_gp, col = col)
  colwidths <- unit.c(prFpFindWidestGrob(labels[[1]]))
  if (nc > 1) {
      for (i in 2:nc) {
          colwidths <- unit.c(colwidths, colgap, prFpFindWidestGrob(labels[[i]]))
      }
  }
  axisList <- prFpGetGraphTicksAndClips(xticks = xticks, xticks.digits = xticks.digits,
      grid = grid, xlog = xlog, xlab = xlab, lwd.xaxis = lwd.xaxis,
      txt_gp = txt_gp, col = col, clip = clip, zero = zero,
      x_range = prFpXrange(upper = upper, lower = lower, clip = clip,
          zero = zero, xticks = xticks, xlog = xlog), mean = org_mean,
      graph.pos = graph.pos, shapes_gp = shapes_gp)
  clip <- axisList$clip
  marList <- list()
  marList$bottom <- convertY(mar[1], "npc")
  marList$left <- convertX(mar[2], "npc")
  marList$top <- convertY(mar[3], "npc")
  marList$right <- convertX(mar[4], "npc")
  prPushMarginViewport(bottom = marList$bottom, left = marList$left,
      top = marList$top, right = marList$right, name = "forestplot_margins")
  if (!missing(title)) {
      prGridPlotTitle(title = title, gp = txt_gp$title)
  }
  if (!missing(legend)) {
      lGrobs <- prFpGetLegendGrobs(legend = legend, txt_gp = txt_gp,
          title = legend_args$title)
      legend_colgap <- colgap
      if (convertUnit(legend_colgap, unitTo = "mm", valueOnly = TRUE) >
          convertUnit(attr(lGrobs, "max_height"), unitTo = "mm",
              valueOnly = TRUE)) {
          legend_colgap <- attr(lGrobs, "max_height")
      }
      legend_horizontal_height <- sum(legend_args$padding,
          attr(lGrobs, "max_height"), legend_args$padding)
      if (!is.null(attr(lGrobs, "title"))) {
          legend_horizontal_height <- sum(attr(lGrobs, "titleHeight"),
              attr(lGrobs, "line_height_and_spacing")[2], legend_horizontal_height)
      }
      legend_vertical_width <- sum(unit.c(legend_args$padding,
          attr(lGrobs, "max_height"), legend_colgap, attr(lGrobs,
              "max_width"), legend_args$padding))
      if ((!is.list(legend_args$pos) && legend_args$pos ==
          "top") || ("align" %in% names(legend_args$pos) &&
          legend_args$pos[["align"]] == "horizontal")) {
          legend_layout <- grid.layout(nrow = 3, ncol = 1,
              heights = unit.c(legend_horizontal_height, legend_colgap +
                legend_colgap, unit(1, "npc") - legend_horizontal_height -
                legend_colgap - legend_colgap))
          legend_pos <- list(row = 1, col = 1)
          main_pos <- list(row = 3, col = 1)
      }
      else {
          legend_layout <- grid.layout(nrow = 1, ncol = 3,
              widths = unit.c(unit(1, "npc") - legend_colgap -
                legend_vertical_width, legend_colgap, legend_vertical_width))
          legend_pos <- list(row = 1, col = 3)
          main_pos <- list(row = 1, col = 1)
      }
  }
  if (!missing(legend) > 0 && !is.list(legend_args$pos)) {
      pushViewport(prFpGetLayoutVP(lineheight = lineheight,
          labels = labels, nr = nr, legend_layout = legend_layout))
      vp <- viewport(layout.pos.row = legend_pos$row, layout.pos.col = legend_pos$col,
          name = "legend")
      pushViewport(vp)
      prFpDrawLegend(lGrobs = lGrobs, col = col, colgap = convertUnit(legend_colgap,
          unitTo = "mm"), pos = legend_args$pos, gp = legend_args$gp,
          r = legend_args$r, padding = legend_args$padding,
          fn.legend = fn.legend, ...)
      upViewport()
      vp <- viewport(layout.pos.row = main_pos$row, layout.pos.col = main_pos$col,
          name = "main")
      pushViewport(vp)
  }
  else {
      pushViewport(prFpGetLayoutVP(lineheight = lineheight,
          labels = labels, nr = nr))
  }
  if (!is.unit(graphwidth) && graphwidth == "auto") {
      npc_colwidths <- convertUnit(unit.c(colwidths, colgap),
          "npc", valueOnly = TRUE)
      graphwidth <- unit(max(0.05, 1 - sum(npc_colwidths)),
          "npc")
  }
  else if (!is.unit(graphwidth)) {
      stop("You have to provide graph width either as a unit() object or as 'auto'.",
          " Auto sizes the graph to maximally use the available space.",
          " If you want to have exact mm width then use graphwidth = unit(34, 'mm').")
  }
  if (graph.pos == 1) {
      colwidths <- unit.c(graphwidth, colgap, colwidths)
  }
  else if (graph.pos == nc + 1) {
      colwidths <- unit.c(colwidths, colgap, graphwidth)
  }
  else {
      spl_position <- ((graph.pos - 1) * 2 - 1)
      colwidths <- unit.c(colwidths[1:spl_position], colgap,
          graphwidth, colwidths[(spl_position + 1):length(colwidths)])
  }
  axis_height <- unit(0, "npc")
  if (is.grob(axisList$axisGrob))
      axis_height <- axis_height + grobHeight(axisList$axisGrob)
  if (is.grob(axisList$labGrob)) {
      gp_lab_cex <- prGetTextGrobCex(axisList$labGrob)
      axis_height <- axis_height + unit(gp_lab_cex + 0.5, "line")
  }
  axis_layout <- grid.layout(nrow = 2, ncol = 1, heights = unit.c(unit(1,
      "npc") - axis_height, axis_height))
  pushViewport(viewport(layout = axis_layout, name = "axis_margin"))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  main_grid_layout <- grid.layout(nrow = nr, ncol = length(colwidths),
      widths = colwidths, heights = unit(rep(1/nr, nr), "npc"),
      respect = TRUE)
  pushViewport(viewport(layout = main_grid_layout, name = "BaseGrid"))
  browser()
  if (!missing(boxsize)) {
      info <- rep(boxsize, length = length(mean))
  }
  else {
      cwidth <- (upper - lower)
      cwidth[cwidth <= 0 | is.na(cwidth)] <- min(cwidth[cwidth >
          0])
      textHeight <- convertUnit(grobHeight(textGrob("A", gp = do.call(gpar,
          txt_gp$label))), unitTo = "npc", valueOnly = TRUE)
      info <- 1/cwidth * 0.75
      info <- info/max(info[!is.summary], na.rm = TRUE)
      if (any(textHeight * (nr + 0.5) * 1.5 < info))
          info <- textHeight * (nr + 0.5) * 1.5 * info/max(info,
              na.rm = TRUE) + textHeight * (nr + 0.5) * 1.5/4
      info[is.summary] <- 1/NCOL(org_mean)
  }
  prFpPrintLabels(labels = labels, nc = nc, nr = nr, graph.pos = graph.pos)
  prFpDrawLines(hrzl_lines = hrzl_lines, nr = nr, colwidths = colwidths,
      graph.pos = graph.pos)
  prFpPrintXaxis(axisList = axisList, col = col, lwd.zero = lwd.zero,
      shapes_gp = shapes_gp)
  for (i in 1:nr) {
      if (is.matrix(org_mean)) {
          low_values <- org_lower[i, ]
          mean_values <- org_mean[i, ]
          up_values <- org_upper[i, ]
          info_values <- matrix(info, ncol = length(low_values))[i,
              ]
      }
      else {
          low_values <- org_lower[i]
          mean_values <- org_mean[i]
          up_values <- org_upper[i]
          info_values <- info[i]
      }
      clr.line <- rep(col$line, length.out = length(low_values))
      clr.marker <- rep(col$box, length.out = length(low_values))
      clr.summary <- rep(col$summary, length.out = length(low_values))
      line_vp <- viewport(layout.pos.row = i, layout.pos.col = graph.pos *
          2 - 1, xscale = axisList$x_range, name = sprintf("Line_%d_%d",
          i, graph.pos * 2 - 1))
      pushViewport(line_vp)
      if (length(low_values) > 1) {
          b_height <- max(info_values)
          if (is.unit(b_height))
              b_height <- convertUnit(b_height, unitTo = "npc",
                valueOnly = TRUE)
          if (missing(line.margin)) {
              line.margin <- 0.1 + 0.2/(length(low_values) -
                1)
          }
          else if (is.unit(line.margin)) {
              line.margin <- convertUnit(line.margin, unitTo = "npc",
                valueOnly = TRUE)
          }
          y.offset_base <- b_height/2 + line.margin
          y.offset_increase <- (1 - line.margin * 2 - b_height)/(length(low_values) -
              1)
          for (j in length(low_values):1) {
              current_y.offset <- y.offset_base + (length(low_values) -
                j) * y.offset_increase
              if (is.na(mean_values[j]))
                next
              shape_coordinates <- c(i, j)
              attr(shape_coordinates, "max.coords") <- c(nr,
                length(low_values))
              if (is.summary[i]) {
                call_list <- list(fn.ci_sum[[i]][[j]], lower_limit = low_values[j],
                  estimate = mean_values[j], upper_limit = up_values[j],
                  size = info_values[j], y.offset = current_y.offset,
                  col = clr.summary[j], shapes_gp = shapes_gp,
                  shape_coordinates = shape_coordinates)
              }
              else {
                call_list <- list(fn.ci_norm[[i]][[j]], lower_limit = low_values[j],
                  estimate = mean_values[j], upper_limit = up_values[j],
                  size = info_values[j], y.offset = current_y.offset,
                  clr.line = clr.line[j], clr.marker = clr.marker[j],
                  lty = lty.ci[[i]][[j]], vertices.height = ci.vertices.height,
                  shapes_gp = shapes_gp, shape_coordinates = shape_coordinates)
                if (!missing(ci.vertices))
                  call_list$vertices = ci.vertices
                if (!missing(lwd.ci))
                  call_list$lwd <- lwd.ci
              }
              if (length(list(...)) > 0) {
                ll <- list(...)
                for (name in names(ll)) {
                  call_list[[name]] <- ll[[name]]
                }
              }
              tryCatch(eval(as.call(call_list)), error = function(e) {
                stop("On row ", i, " the print of the estimate failed: ",
                  e$message)
              })
          }
      }
      else {
          shape_coordinates <- c(i, 1)
          attr(shape_coordinates, "max.coords") <- c(nr, 1)
          if (is.summary[i]) {
              call_list <- list(fn.ci_sum[[i]], lower_limit = low_values,
                estimate = mean_values, upper_limit = up_values,
                size = info_values, col = clr.summary, shapes_gp = shapes_gp,
                shape_coordinates = shape_coordinates)
          }
          else {
              call_list <- list(fn.ci_norm[[i]], lower_limit = low_values,
                estimate = mean_values, upper_limit = up_values,
                size = info_values, clr.line = clr.line, clr.marker = clr.marker,
                lty = lty.ci[[i]], vertices.height = ci.vertices.height,
                shapes_gp = shapes_gp, shape_coordinates = shape_coordinates)
              if (!missing(ci.vertices))
                call_list$vertices = ci.vertices
              if (!missing(lwd.ci))
                call_list$lwd <- lwd.ci
          }
          if (length(list(...)) > 0) {
              ll <- list(...)
              for (name in names(ll)) {
                call_list[[name]] <- ll[[name]]
              }
          }
          if (!is.na(mean_values)) {
              tryCatch(eval(as.call(call_list)), error = function(e) {
                stop("On row ", i, " the print of the estimate failed: ",
                  e$message)
              })
          }
      }
      upViewport()
  }
  if (!missing(legend) && is.list(legend_args$pos)) {
      plot_vp <- viewport(layout.pos.row = 1:nr, layout.pos.col = 2 *
          graph.pos - 1, name = "main_plot_area")
      pushViewport(plot_vp)
      if ("align" %in% names(legend_args$pos) && legend_args$pos[["align"]] ==
          "horizontal") {
          height <- legend_horizontal_height
          width <- 0
          for (i in 1:length(lGrobs)) {
              if (width > 0) {
                width <- width + convertUnit(legend_colgap,
                  unitTo = "npc", valueOnly = TRUE)
              }
              width <- width + convertUnit(attr(lGrobs, "max_height") +
                legend_colgap + attr(lGrobs[[i]], "width"),
                unitTo = "npc", valueOnly = TRUE)
          }
          width <- unit(width + convertUnit(legend_args$padding,
              unitTo = "npc", valueOnly = TRUE) * 2, "npc")
      }
      else {
          legend_height <- attr(lGrobs, "line_height_and_spacing")[rep(1:2,
              length.out = length(legend) * 2 - 1)]
          if (!is.null(attr(lGrobs, "title"))) {
              legend_height <- unit.c(attr(lGrobs, "titleHeight"),
                attr(lGrobs, "line_height_and_spacing")[2],
                legend_height)
          }
          height <- sum(legend_args$padding, legend_height,
              legend_args$padding)
          width <- legend_vertical_width
      }
      pushViewport(viewport(x = legend_args$pos[["x"]], y = legend_args$pos[["y"]],
          width = width, height = height, just = legend_args$pos[["just"]]))
      prFpDrawLegend(lGrobs = lGrobs, col = col, colgap = legend_colgap,
          pos = legend_args$pos, gp = legend_args$gp, r = legend_args$r,
          padding = legend_args$padding, fn.legend = fn.legend,
          ...)
      upViewport(2)
  }
  seekViewport("forestplot_margins")
  upViewport(2)
})

