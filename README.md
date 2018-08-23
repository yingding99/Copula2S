# Copula2S contains the key R functions for implementing the two-parameter copula-based sieve semiparametric transformation model for bivariate interval-censored data. An example dataset is included to demonstrate how these functions can be called to analyze the data.
1.	model_functions.R: contains two key functions for constructing the model log-likelihood;
2.	step1a.R: contains the marginal sieve log-likelihood function used in step 1a of the proposed estimation procedure to obtain initial estimates of the regression parameters and Bernstein basis polynomial coefficients;
3.	step1b.R: contains the pseudo-log-likelihood function used in step 1b to obtain pseudo-MLE of the copula dependence parameters;
4.	step2.R: contains the joint log-likelihood function used in step 2 to update all model parameters;
5.	lik.R: returns the full joint log-likelihood function;
6.	main.R demonstrate the entire estimation procedure with a generalized score test step-by-step. It outputs the sieve estimators as well as the p-value(s) from the score test.

The dataset “data_example.RData” contains a bivariate interval-censored dataset for 500 subjects. Each row gives data for one eye, with three covariates (one continuous, one binary and one SNP), the observed time interval (“Left” and “Right”), and a censoring indicator “status” (0 for right-censored and 1 for interval-censored).
