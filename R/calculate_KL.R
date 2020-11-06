#' Calculate a modified KL score that uses absolute value of the ratios
#' It will detect changes in the relative proportion of enzymes connected to a metabolite
#' @param x expression vector for the reference/baseline condition, in FPKM, or TPM value - normalized for depth and gene length.
#' @param y expression vector for the perturbed/treated condition, in FPKM, or TPM value - normalized for depth and gene length.
#' @return KL value denoting the difference between the distribution of values in x and y
calculateKL_abs<-function(x,y){
  x<-x+0.1#to avoid zeros
  y<-y+0.1
  x<-x/sum(x)
  y<-y/sum(y)
  KL<-sum(abs(y*log2(y/x)))
  return(KL)
}


#' Calculate a modified KL score that uses absolute value of the ratios
#' It will detect cases where the relative proportion did not change, but the total amount of expressed enzymes changed.
#' @param x expression vector for the reference/baseline condition, in FPKM, or TPM value - normalized for depth and gene length.
#' @param y expression vector for the perturbed/treated condition, in FPKM, or TPM value - normalized for depth and gene length.
#' @return KL value denoting the difference between the distribution of values in x and y

calculateKL_Dirichlet<-function(x,y){
  x<-x+1
  y<-y+1
  sum_x<-sum(x)
  sum_y<-sum(y)
  max_sum<-max(c(sum_x,sum_y))
  a<- (sum_x/max_sum)*(x/sum(x))
  b<- (sum_y/max_sum)*(y/sum(y))
  KL<-lgamma(sum(a)) - lgamma(sum(b)) - sum(lgamma(a)) + sum(lgamma(b)) + sum((a-b)*(digamma(a) - digamma(sum(a))))
  return(KL)
}
