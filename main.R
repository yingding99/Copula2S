##############################################################################################
### Estimation and Inference under the two-parameter Copula-Sieve transformation mdoel ###
###############################################################################################

source("model_functions.R") 
source("step1a.R") 
source("step1b.R") 
source("step2.R") 
source("lik.R") 
library("stats") 
library("pracma")
library("corpcor")

# data input
load("data_example.RData") # 1,000 rows of observations
var_list=c("SevScaleBL","SEX","SNP_1") 
p=length(var_list) # number of unknown beta = 3
n = nrow(data_example)/2 # cluster size = 500

# assign a large number to R_ij = Inf
max_end = max(data_example$Left,max(data_example$Right[is.finite(data_example$Right)]))
data_example$Right[data_example$status==0] = max_end + 1

# separate into two margins
indata1 = subset(data_example, ind==1)
indata2 = subset(data_example, ind==2)
dim_m<-dim(as.matrix(indata1[,var_list]))
x1<-matrix(as.numeric(as.character(as.matrix(indata1[,var_list]))),dim_m[1],dim_m[2])
x2<-matrix(as.numeric(as.character(as.matrix(indata2[,var_list]))),dim_m[1],dim_m[2])   
t_left <- data_example[,"Left"]
t_right <- data_example[,"Right"]
status1 <- indata1$status # censoring status for 1st margin (0 for RC, 1 for IC)
status2 <- indata2$status

# specify range [l,u]
l=min(data_example$Left, data_example$Right)
u=max(data_example$Left, data_example$Right)

# specify transformation function and Bernstein degree
m=3 # degrees
r=3 # r = 3 for PO; r=1 for PH
G <- function(para) {
  G_fun(para,3)
}

# setup Bernstein basis polynomials (assume the same bais between two margins here)
bl=br=matrix(0,nrow = 2*n,ncol = m+1)
for (i in 0:m) {
  bl[,(i+1)] = bern(i,m,l,u,t_left)
  br[,(i+1)] = bern(i,m,l,u,t_right)
}
odd_index <- seq(1,2*n,by=2)
even_index <- seq(2,2*n,by=2)
bl1 <- bl[odd_index,] #1st margin, left end
br1 <- br[odd_index,] #1st margin, right end
bl2 <- bl[even_index,] #2nd margin, left end
br2 <- br[even_index,] #2nd margin, right end



##################################################################################################
### Fitting Copula2-Sieve using two covariates SevScaleBL and Sex (set SNP_1 effect = 0)
# number of regression coefs = p - 1
# number of Bernstein coefs = m + 1 (assuming same covariate effects between two margins)
##################################################################################################
#### step 1a: marginal model
model_step1a = nlm(step1a,p = c(rep(0,p-1),rep(0,m+1)),hessian = F, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2) 
phi = model_step1a$estimate[((p-1)+1):((p-1)+1+m)] 
ep <- cumsum(exp(phi)) # reparameterization back to Bernstein coefficients
    
####  step 1b: update alpha and kappa only, others are fixed at previous beta and ep 
model_step1b = nlm(step1b,p = c(1,1),hessian = F, x1, x2, bl1, br1, bl2, br2, status1, status2, 
                   beta = c(model_step1a$estimate[1:(p-1)],0), ep)
    
### step 2, with all starting values of p, optimize all parameters at the same time 
model_step2 = nlm(step2,p = c(model_step1a$estimate,model_step1b$estimate),hessian = T,  x1, x2, bl1, br1, bl2, br2, status1, status2)
est = c(model_step2$estimate[1:(p-1)],model_step2$estimate[(length(model_step2$estimate)-1):length(model_step2$estimate)])
se = sqrt(diag(pseudoinverse(model_step2$hessian))[c(1:(p-1), (length(model_step2$estimate)-1):length(model_step2$estimate))])
p_value_wald <- pchisq((est/se)^2,1, lower.tail = F)
output <- data.frame(estimates = est, se = se, p = p_value_wald)
rownames(output) <- c(var_list[1:(p-1)], "alpha", "kappa")
output 
# estimates         se            p
# SevScaleBL 0.08771905 0.01660505 1.273131e-07
# SEX        0.09142112 0.13313663 4.922898e-01
# alpha      0.96303266 0.07616787 1.214372e-36
# kappa      0.37531307 0.05776231 8.164235e-11


############################################################################
### Generalized score test Under Null (H0: effect of SNP_1 =0)
############################################################################
estimates = c(model_step2$estimate[1:(p-1)],0,model_step2$estimate[((p-1)+1):((p-1)+1+m+1+1)]) # including SNP_1 = 0
score = grad(lik, x0 = estimates, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, status1=status1, status2=status2) # score function by numerical approximation
hes = hessian(lik, x0 = estimates, x1=x1, x2=x2, bl1=bl1, br1=br1, bl2=bl2, br2=br2, status1=status1, status2=status2) # hessian matrix by numeircal approximation
test_stat = t(score) %*% pseudoinverse(hes) %*% score  # score test statistics
test_stat
# 0.5870823
p_value_score <- pchisq(test_stat,1, lower.tail = F) # score test p value
p_value_score
# 0.4435496  

    

