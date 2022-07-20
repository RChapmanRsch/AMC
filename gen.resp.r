gen.resp <- function(ip, x = NULL, D = 1.7) {
  p <- lab.irf(ip, x, D)$f
  if (is.null(dim(p))) dim(p) <- c(1, length(p))
  np <- dim(p)[1]
  ni <- dim(p)[2]
  u <- matrix(runif(np * ni), ncol = ni)
  re <- ifelse(p > u, 1, 0)
  colnames(re) <- paste0("item", 1:ncol(re))
  return(re)
} # end function.
