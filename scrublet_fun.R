#' @title Doublet score
#' @description Calculate a doublet score and predict which cells are doublets.
#' @param edata Expression data (raw, unnormalized)
#' @param mdata Metadata, Default: NULL
#' @param by Library column, Default: NULL
#' @param score_threshold Score to call a cell as a doublet, Default: 0.2
#' @param sim_doublet_ratio Expected doublet ratio, Default: 2
#' @param n_neighbors Number of neighbors used to construct the KNN graph,
#' Default: NA (round(0.5 * sqrt(colnames(edata))).
#' @param expected_doublet_rate Based on 10x guides, Default: 0.1
#' @param do_parallel Do it in parallel.
#' @param verbose Show progress, Default: 2.
#' @param ... More parameters for scrub_doublets (you need to specify the type
#' variable, ie., as.integer)
#' @return Metadata with extra doublet information.
#' @examples
#' \dontrun{
#'  if(interactive()) mscrublet = doublets_scrublet(edata = counts)
#' }
#' @references \url{https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb}
#' @seealso
#'  \code{\link[reticulate]{import}}
#'  \code{\link[Matrix]{t}}
#' @rdname doublets
#' @export
#' @importFrom reticulate import
#' @importFrom Matrix t
#' @importFrom foreach %dopar% foreach
doublets_scrublet = function(
  edata,
  mdata = NULL,
  by = NULL,
  score_threshold = 0.2,
  sim_doublet_ratio = 2.0, # Run
  n_neighbors = NA,
  expected_doublet_rate = 0.06,
  do_parallel = TRUE,
  verbose = TRUE,
  ...
) {
  scr = reticulate::import("scrublet")

  if(is.null(mdata))
    mdata = data.frame(AllData_ = rep("Data", ncol(edata)), row.names = colnames(edata))
  if(is.null(by)){ mdata$AllData_ = "Data"; by = "AllData_" }
  if(verbose) cat("By:", by, "\n")
  library_id = unique(mdata[, by])
  edata = edata[, rownames(mdata)]
  if(length(expected_doublet_rate) == 1){
    expected_doublet_rate = rep(expected_doublet_rate, length(library_id))
  }
  if(!is.null(names(expected_doublet_rate)))
    expected_doublet_rate = expected_doublet_rate[library_id]
  # if(length(n_neighbors) == 1) n_neighbors = rep(n_neighbors, 3)
  # if(!is.null(names(n_neighbors))) n_neighbors = n_neighbors[library_id]

  do_par = all(requireNamespace(c("doParallel", "parallel"), quietly = TRUE))
  if(do_par && length(library_id) > 2 && isTRUE(do_parallel)){
    doParallel::registerDoParallel(cores = min(c(4, parallel::detectCores())))
  }else{ foreach::registerDoSEQ() }
  doublets = foreach(i = library_id, .combine = rbind) %dopar% {
    edata_i = Matrix::t(edata[, mdata[, by] == i])
    # n_nei = n_neighbors[which(library_id %in% i)]
    # if(is.na(n_nei)) n_nei = 0.5 * sqrt(ncol(edata_i))
    scrub_object = scr$Scrublet(
      counts_matrix = edata_i,
      sim_doublet_ratio = sim_doublet_ratio,
      # n_neighbors = as.integer(round(n_nei)),
      expected_doublet_rate = expected_doublet_rate[which(library_id %in% i)]
    )
    doublet_scores = scrub_object$scrub_doublets(...)
    if(is.null(doublet_scores[[2]]))
      doublet_scores[[2]] = scrub_object$call_doublets(threshold=score_threshold)
    # if(length(dev.list())) scrub_object$plot_histogram();
    names(doublet_scores) = c("doublet_scores", "predicted_doublets_auto")
    doublet_scores = as.data.frame(doublet_scores)
    doublet_scores$predicted_doublets = doublet_scores$doublet_scores > score_threshold
    rownames(doublet_scores) = rownames(edata_i)
    str(doublet_scores)
    doublet_scores
  }
  mdata = cbind(mdata, doublets[rownames(mdata), ])
  mdata = mdata[, !colnames(mdata) %in% "AllData_"]
  return(mdata)
}

#' @title Plot doublet score
#' @description Plot histograms of scores.
#' @param mdata Metadata.
#' @param by Library column, Default: NULL.
#' @param score_threshold Threshold to call a doublet or it can also
#' be a named vector, Default: 0.2.
#' @param verbose Show progress, Default: 2.
#' @param ... Extra qplot arguments.
#' @return ggplot from cowplot
#' @examples
#' \dontrun{
#'  if(interactive()) mscrublet = doublets_scrublet(edata = counts)
#' }
#' @seealso
#'  \code{\link[ggplot2]{qplot}}
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname doublets
#' @export
#' @importFrom ggplot2 qplot
#' @importFrom cowplot plot_grid
doublets_scrublet_plot = function(
  mdata,
  by = NULL,
  score_threshold = NULL,
  verbose = TRUE,
  ...
) {
  if(is.null(by)){ mdata$AllData_ = "Data"; by = "AllData_" }
  if(verbose) cat("By:", by, "\n")
  p_list = lapply(X = unique(mdata[, by]), function(i) {
    if(verbose) cat(" -", i, "\n")
    mdata_i = mdata[mdata[, by] == i, ]
    ythr = if(is.null(score_threshold)){
      min(mdata_i$doublet_scores[mdata_i$predicted_doublets_auto])
    }else if(isTRUE(i %in% names(score_threshold))){
      score_threshold[i][1]
    }else{ score_threshold }
    if(grepl("^p", ythr))
      ythr <- quantile(mdata_i$doublet_scores, 1 - as.numeric(sub("p", "", ythr)))
    n_prop <- mdata_i$doublet_scores > ythr
    mdata_stitle = paste0(
      "Threshold: ", round(ythr, 3), " | N = ",
      sum(n_prop), " (", round(mean(n_prop)*100, 2), "%)"
    )
    if(verbose) cat(mdata_stitle, "\n")

    ggplot2::qplot(mdata_i$doublet_scores, data = mdata_i, geom = "histogram", bins = 50, ...) +
      ggplot2::geom_vline(xintercept = ythr, linetype = "dashed", color = "#007385") +
      ggplot2::theme_minimal() + ggplot2::scale_y_log10() + ggplot2::xlim(-0.05, 1.05) +
      ggplot2::labs(title = i, subtitle = mdata_stitle, x = "Score")
  })
  cowplot::plot_grid(plotlist = p_list)
}
