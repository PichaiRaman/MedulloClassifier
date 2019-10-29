createRatio <- function(x) {
  g1 <- x[1];
  g2 <- x[2];
  g1g2_ratio <- 2^(exprs_37418[g1,]-exprs_37418[g2,])
  return(g1g2_ratio)
}