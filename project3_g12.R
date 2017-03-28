# -------------------- #
# Mac install gurobi
# -------------------- #

# load("~/Downloads/data.rdata") # Ann Path
# install.packages('/Library/gurobi702/mac64/R/gurobi_7.0-2.tgz', repos=NULL,type="source")
# change the directory if you are not using MAX OS
# install.packages("glmnet", repos = "http://cran.us.r-project.org")

# -------------------- #
# Windows install gurobi
# -------------------- #

#install.packages('slam')
#install.packages('c:\\gurobi702\\win64\\R\\gurobi_7.0-2.zip', repos=NULL)

#setwd("~/5 - MSBA/2017 01 Optimization/Group Projects/Project 3") # Nicole Path
load("data.rdata")

# -------------------- #
# ---- Question 2 ---- #
# -------------------- #

# Indirect Feature Selection (Lasso Regression)
library(glmnet)
fit = glmnet(X, y, alpha = 1)
plot(fit, xvar = "lambda", label = TRUE)

# The plot shows a cut-off at Log Lambda = -2; we will use e^(-2) as the penalty term of our Lasso model
lasso_beta = coef(fit,s=exp(-2))[-1]
lasso_real_comp = cbind(lasso_beta,beta_real)

# Comparing the real betas and the results returned by Lasso, we can see that Lasso successfully regularized the irrelevant features by shrinking their corresponding coefficients towards zero  
lasso_real_comp

# Remove rows that are both 0
lasso_real_comp = as.matrix(subset(lasso_real_comp,subset=lasso_beta>0 | beta_real>0))

# Direct Selection (MIP)
construct_MIQP = function(X,y,k,M){
  # this function is to formulate the Mixed Interger Programming Problem
  n = dim(X)[2]
  # n = 64 in our case
  
  # we created 128 variables in total:
  # the first 64 (B1, B2,.., B64) are continuous variables representing each independent variable; 
  # the last 64 (Z1, Z2,.., Z64) are binary variables: Zi (i in 1:64) indicates whether Bi is 0
  model <- list()
  model$vtype <- c(rep('C',n),rep('B',n))
  
  # Formulate Constraints
  A = matrix(0,2*n,2*n)
  # 1. for -M*Zi <= Bi <= M*Zi:
  A[1:n,1:n] = -1*diag(n)
  A[1:n,(n+1):(2*n)] = -M*diag(n)
  A[(n+1):(2*n),1:n] = diag(n)
  A[(n+1):(2*n),(n+1):(2*n)] = -M*diag(n)
  # 2. for Z1 + Z2 + ... + Z64 = k:
  A1 = c(rep(0,n),rep(1,n))
  model$A <- rbind(A1,A)
  model$sense <- rep("<=",2*n+1)
  model$rhs <- c(k,rep(0,2*n))
  
  # Formulate Objective (i.e. sum of squared residuals in our case)
  Q = matrix(0,2*n,2*n)
  # 1. quadratic component of the objective function
  Q[0:n,0:n] = t(X) %*% X
  model$Q = Q
  # 2. linear component of the objective function
  model$obj <- c(-2*y %*% X,rep(0,n))
  
  result <- gurobi(model, list(ResultFile='model.mps',OutputFlag=0))
  return(result)
}

library(gurobi)
k=8
M = 0.1
initial_sol = construct_MIQP(X,y,k,M)
max = max(abs(initial_sol$x[1:64]))
max

while (M<=max){
  M = M*2
  current_sol = construct_MIQP(X,y,k,M)
  max = max(abs(current_sol$x[1:64]))
  print(M)
}

MIQP_beta = current_sol$x[1:64]
sol_compare = cbind(lasso_beta,MIQP_beta,beta_real)
sol_compare

# Find all MIQP non-zero coefficients
MIQP_non_zero = as.matrix(cbind(match(MIQP_beta[MIQP_beta>0],MIQP_beta),MIQP_beta[MIQP_beta>0]))
colnames(MIQP_non_zero) = c("Coef. Number","MIQP_beta")

MIQP_real_comp = as.matrix(cbind(MIQP_beta,beta_real))
# Remove rows that are both 0
MIQP_real_comp = as.matrix(subset(MIQP_real_comp,subset=MIQP_beta>0 | beta_real>0))

# -------------------- #
# ---- Question 3 ---- #
# -------------------- #

compute_error = function(sol_beta, real_beta, X){
  error = X%*%sol_beta - X%*%real_beta
  return(sum(error^2)/sum((X%*%real_beta)^2))
}
compute_error(beta_real, beta_real, X)
compute_error(lasso_beta, beta_real, X)
compute_error(MIQP_beta, beta_real, X)


