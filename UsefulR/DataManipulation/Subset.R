# Conditionally subset the rownames
sub.row <- function(x, z.r=0.001, z.a=0.01){
  x.prop <- x/colSums(x)
  id.cat <- c()
  p.sub <- function(k){
    record1 <- sum(k>z.r)
    record2 <- sum(k>z.a)
    return( (record1==length(k))&(record2>0)&(record2<length(k)) )
  }
  sub.row.vec <- apply(x.prop, 1, p.sub)
  return(rownames(x)[sub.row.vec==TRUE])
}
