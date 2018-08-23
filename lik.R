# Two-parameter copual sieve model full likelihood function
# para include beta, phi (bernstein coefficients), copula parameters (alpha, kappa)
# x1, x2 are design matrices of two margins
# bl1, br1: Bernstein basis polynomials for Left and Right endpoints of oberved intervals for the 1st margin; similar to bl2 and br2
# status1, status2: censoring status for two margins, with 0 for right censor and 1 for interval censor 
# return -loglik

lik <- function(para, x1, x2, bl1, br1, bl2, br2, status1, status2) { 
  
  beta<-para[1:p]
  phi<-para[(p+1):(p+1+m)] 
  alpha<-para[p+1+m+1] 
  kappa<-para[p+1+m+1+1] 
  
  ep<-cumsum(exp(phi))
  
  gLl1<-G(exp(x1%*%beta)*(bl1%*%ep)) 
  gLr1<-G(exp(x1%*%beta)*(br1%*%ep)) 
  gLl2<-G(exp(x2%*%beta)*(bl2%*%ep)) 
  gLr2<-G(exp(x2%*%beta)*(br2%*%ep))
  
  u1_left = exp(-gLl1)
  u1_right = exp(-gLr1)
  u2_left = exp(-gLl2)
  u2_right = exp(-gLr2)
  
  # Copula function for joint distribution probability
  C_val_1 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2_left^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
  C_val_2 <- (1 + ((u1_left^(-1/kappa)-1)^(1/alpha) + (u2_right^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
  C_val_3 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2_left^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
  C_val_4 <- (1 + ((u1_right^(-1/kappa)-1)^(1/alpha) + (u2_right^(-1/kappa)-1)^(1/alpha))^alpha)^(-kappa)
  
  # Use Copula functions to write each block of likelihood function
  term1 <- log(C_val_1 - C_val_2 - C_val_3 + C_val_4)
  term1 <- ifelse((status1 == 1) & (status2 == 1), term1, 0)
  
  term2 <- log(C_val_1 - C_val_3)
  term2 <- ifelse((status1 == 1) & (status2 == 0), term2, 0)
  
  term3 <- log(C_val_1 - C_val_2)
  term3 <- ifelse((status1 == 0) & (status2 == 1), term3, 0)
  
  term4 <- log(C_val_1)
  term4 <- ifelse((status1 == 0) & (status2 == 0), term4, 0)
  
  logL<-(-1)*sum( term1 + term2 + term3 + term4 )
  
  return(logL)
  
}