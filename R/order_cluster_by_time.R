#' Order clusters temporally
#'
#' Derive new cluster identities by progression through pseudotime
#'
#' @param embedding data.frame
#' @param clusters factor
#' @param pseudotime numeric
#'
#' @return data.frame
#' @export

order_clusters_by_time <- function(embedding,clusters,pseudotime) {

  #Reorder clusters by progression through pseudotime
  mean.pseudotime.cluster <- fit.frame %>%
    group_by(cluster) %>%
    summarize(mean_lambda=mean(pseudotime))

  mean.pseudotime.cluster <- mean.pseudotime.cluster %>%
    arrange(mean_lambda) %>% mutate(new_cluster_id=seq_along(1:nrow(mean.pseudotime.cluster)))

  level.key <- c(rd3.means$new_cluster_id)
  names(level.key) <- rd3.means$clusters
  pca.plot2 <- pca.plot %>% mutate(new_clusters=recode(clusters,!!!level.key))

  cluster.means <- left_join(cluster.means,rd3.means,by="clusters")
  cluster.means <- cluster.means %>% arrange(new_cluster_id)
  cluster.means <- rbind(cluster.means,cluster.means[1,])

}
