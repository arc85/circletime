#' Fit genes as a function of pseudotime
#'
#' Utilize generalized additive models to fit gene expression as a function of pseudotime.
#'
#' @param expr.matrix data.frame
#' @param pseudotime numeric
#'
#' @return list
#' @export

genes_over_time <- function(expr.matrix,pseudotime) {

expr.list <- split(expr.matrix,rownames(expr.matrix))
t <- pseudotime

p.val <- lapply(expr.list,function(z) {
  d <- data.frame(z=z,t=t)
  tmp <- gam::gam(z~gam::lo(t),data=d)
  p <- summary(tmp)[4][[1]][1,5]
})

fits <- lapply(expr.list,function(z) {
  d <- data.frame(z=z,t=t)
  tmp <- gam::gam(z~gam::lo(t),data=d)
  stats::predict(tmp,data=d)
})

names(fits) <- names(p.val) <- rownames(expr.matrix)

list(p_values=p.val,expression_fits=fits)

}
