createRatio <- function(exprs, x) {
  g1 <- x[1];
  g2 <- x[2];
  g1g2_ratio <- 2^(exprs[g1,]-exprs[g2,])
  return(g1g2_ratio)
}
