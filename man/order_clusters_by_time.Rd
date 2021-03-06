% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/order_clusters_by_time.R
\name{order_clusters_by_time}
\alias{order_clusters_by_time}
\title{Order clusters temporally}
\usage{
order_clusters_by_time(embedding, clusters, pseudotime)
}
\arguments{
\item{embedding}{data.frame}

\item{clusters}{factor}

\item{pseudotime}{numeric}
}
\value{
data.frame
}
\description{
Derive new cluster identities by progression through pseudotime. This function returns a new data.frame with clusters ordered by pseudotemporal progression.
}
