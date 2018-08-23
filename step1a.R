# Marginal sieve model likelihood function under H0: SNP_1 effect = 0
# para include beta, phi (bernstein coefficients)
# x1, x2 are design matrices of two margins
# bl1, br1: Bernstein basis polynomials for Left and Right endpoints of oberved intervals for the 1st margin; similar to bl2 and br2
# return -loglik

step1a <- function(para, x1, x2, bl1, br1, bl2, br2) { #total 6 parameters = 2non-gene_beta + 4polynomial
  
  beta<-c(para[1:(p-1)],0) # set SNP as 0
  phi<-para[((p-1)+1):((p-1)+1+m)]  
  
  ep<-cumsum(exp(phi)) # reparameterization back to Bernstein coefficients
  
  gLl1<-G(exp(x1%*%beta)*(bl1%*%ep)) 
  gLr1<-G(exp(x1%*%beta)*(br1%*%ep)) 
  gLl2<-G(exp(x2%*%beta)*(bl2%*%ep)) 
  gLr2<-G(exp(x2%*%beta)*(br2%*%ep))
  
  logL1<-sum(log(exp(-gLl1)-exp(-gLr1)))
  logL2<-sum(log(exp(-gLl2)-exp(-gLr2)))
  logL <- -1*(logL1+logL2)
  return(logL)
}