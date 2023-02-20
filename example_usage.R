# load functions and libraries
source("functions.R")
library(R.utils)

# simulation parameters
n = 400  # sample size
p = 100  # dimension of precision matrix
num_repeat = 100 # number of repetitions
method_id = 0 # which method to run, SGM when 0

# True precision matrix generation parameters
prob = 0.05 
cond = p
n_blocks = 5
single_offdiag_element = c(0.3)
triple_offdiag_elements = c(0.4,0.3,0.2)

# SGM parameters 
alpha1 = sqrt(log(p)/n)
alpha2 = sqrt(log(p)/n)/20

# simulation loop
set.seed(2018)
sgm_estimates = list()
true_precision_matrices = list()
for(i in 1:num_repeat){
  print(i)
  # multiple methods for generating the true precision matrix
  # CHOOSE ONE and comment out the others
  theta_star = genTheta_random_block(p,prob,cond,n_blocks) # option 1
  theta_star = genTheta_random(p,prob,cond = cond) # option 2
  theta_star = genTheta(p, single_offdiag_element) # option 3
  theta_star = genTheta(p, triple_offdiag_elements) # option 4
  # end of true precision matrix generation

  # Generate sample covariance matrices
  svd_out = svd(theta_star)
  X1 = matrix(rnorm(n*p), nrow=n, ncol=p)%*%diag(1/sqrt(svd_out$d))%*%t(svd_out$v)
  X2 = matrix(rnorm(n*p), nrow=n, ncol=p)%*%diag(1/sqrt(svd_out$d))%*%t(svd_out$v)
  Sigma = var(X1)
  Sigma_ref = var(X2)

  # SGM estimate step 1
  out1 = step1estimator(Sigma, lambda.length=100, alpha=alpha1, maxit=1e3, error=1e-2)
  outmin1 = find_min(out1$precisionList, Sigma_ref=Sigma_ref, loss_function="Dtrace")
  Theta1 = outmin1$precision
  # SGM estimate step 2
  weight_sgb = weight_matrix(Theta1, 1/n)
  out2 = step2estimator(Sigma, lambda.length=100, alpha=alpha2, weight=weight_sgb,
                        maxit=1e3, error=1e-2)
  outmin2 = find_min(out2$precisionList, Sigma_ref=Sigma_ref, loss_function="Dtrace")
  Theta2 = outmin2$precision

  # save estimates and true precision matrices
  sgm_estimates[[i]] = Theta2
  true_precision_matrices[[i]] = theta_star
}

