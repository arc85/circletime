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
    dplyr::group_by(clusters) %>%
    dplyr::summarize(mean_lambda=mean(pseudotime))

  mean.pseudotime.cluster <- mean.pseudotime.cluster %>%
    dplyr::arrange(mean_lambda) %>% dplyr::mutate(new_cluster_id=seq_along(1:nrow(mean.pseudotime.cluster)))

  #Calculate embedding means of each cluster
  embeds.calc.mean <- grep("clusters|pseudotime",colnames(fit.frame),invert=TRUE,value=TRUE)

  cluster.split <- split(fit.frame,fit.frame$clusters)

  cluster.means <- t(sapply(cluster.split,function(x) {
    apply(x[,embeds.calc.mean],2,mean)
  }))

  colnames(cluster.means) <- paste(colnames(cluster.means),"mean",sep="_")

  cluster.means <- tibble::as_tibble(cluster.means,rownames="clusters")
  cluster.means$clusters <- as.factor(cluster.means$clusters)

  mean.pseudotime.cluster <- dplyr::left_join(mean.pseudotime.cluster,cluster.means,by="clusters")

  mean.pseudotime.cluster <- rbind(mean.pseudotime.cluster,mean.pseudotime.cluster[1,])

  return(data.frame(mean.pseudotime.cluster))

}
