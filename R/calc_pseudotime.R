#' Calculate pseudotime
#'
#' Calculate pseudotemporal order via principal curves for an arbitrary low dimensional embedding for scRNAseq data.
#'
#' @param embedding data.frame
#' @param clusters factor
#'
#' @return data.frame
#' @export

calc_pseudotime <- function(embedding,clusters) {

  if (class(embedding)=="data.frame") {
    embedding <- as.matrix(embedding)
  }

  row.names <- rownames(embedding)

  #Calculate pseudotime
  prin.fit <- princurve::principal_curve(embedding,smoother="periodic_lowess")

  #Return data.frame of pseudotime results
  fit.frame <- data.frame(cell_order=prin.fit$ord,pseudotime=prin.fit$lambda)
  rownames(fit.frame) <- row.names

  return(fit.frame)

}
