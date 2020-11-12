#' Order clusters temporally
#'
#' Derive new cluster identities by progression through pseudotime. This function returns a new data.frame with clusters ordered by pseudotemporal progression.
#'
#' @param embedding data.frame
#' @param clusters factor
#' @param pseudotime numeric
#'
#' @importFrom magrittr "%>%"
#' @return data.frame
#' @export

order_clusters_by_time <- function(embedding,clusters,pseudotime) {

  fit.frame <- data.frame(embedding,clusters,pseudotime)

  #Reorder clusters by progression through pseudotime
  mean.pseudotime.cluster <- fit.frame %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(mean_lambda=mean(pseudotime))

  mean.pseudotime.cluster <- mean.pseudotime.cluster %>%
    dplyr::arrange(mean_lambda) %>% dplyr::mutate(new_cluster_id=seq_along(1:nrow(mean.pseudotime.cluster)))

  return(mean.pseudotime.cluster)

}
