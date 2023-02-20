# system("R CMD SHLIB comb.c")
dyn.load("comb.so")

lamMax_sgb = function(Sigma, alpha, weight){
  p = dim(Sigma)[1]
  sig = diag(Sigma)
  mat = matrix(rep(1,p))%*%matrix(1/(sig+alpha),1,p)
  mat = mat + t(mat)
  mat = abs(mat*Sigma/2/weight)
  diag(mat) = 0
  return(max(mat))
}
lamMax_glasso = function(Sigma, weight){
  mat = abs(Sigma/weight)
  diag(mat) = 0
  return(max(mat))
}
gradient_check1 = function(Sigma, Theta.List, lambda.list, alpha){
  ans = 0
  if(is.vector(lambda.list)) lamLength = length(lambda.list)
  else lamLength = 1
  for(i in 1:lamLength){
    if(is.vector(lambda.list)){
      Theta = Theta.List[,,i]
      lambda = lambda.list[i]
    }
    else{
      Theta = Theta.List
      lambda = lambda.list
    }
    SigmaTheta = Sigma %*% Theta
    tmp = SigmaTheta + alpha*Theta - diag(dim(Sigma)[1])
    zeroIndex = (Theta == 0)
    submax = abs(tmp[zeroIndex]) - lambda
    ans = max(submax, ans)
    submax = abs(tmp[!zeroIndex] + lambda*sign(Theta[!zeroIndex]))
    ans = max(submax, ans)
  }
  return(ans)
}
gradient_check2 = function(Sigma, Theta.List, lambda.list, alpha, weight){
  ans = 0
  if(is.vector(lambda.list)){
    lamLength = length(lambda.list)
  } else {
    lamLength = 1
  }
  for(i in 1:lamLength){
    if(is.vector(lambda.list)){
      Theta = Theta.List[,,i]
      lambda = lambda.list[i]
    }
    else{
      Theta = Theta.List
      lambda = lambda.list
    }
    SigmaTheta = Sigma %*% Theta
    tmp = (SigmaTheta+t(SigmaTheta))/2 + alpha*Theta - diag(dim(Sigma)[1])
    zeroIndex = (Theta == 0)
    submax = abs(tmp[zeroIndex]) - lambda*weight[zeroIndex]
    ans = max(submax, ans)
    submax = abs(tmp[!zeroIndex] + lambda*weight[!zeroIndex]*sign(Theta[!zeroIndex]))
    ans = max(submax, ans)
  }
  return(ans)
}
step1estimator = function(Sigma, lambda.list = NULL, lambda.length = 100, alpha, maxit, error){
  p = dim(Sigma)[1]
  if(missing(lambda.list)){
    lambdaMax = 1
    lambda.list = rev(seq(lambdaMax/lambda.length, lambdaMax, length.out = lambda.length))
  }
  else lambda.length = length(lambda.list)
  Sigma = as.vector(Sigma)
  out = .C("step1matrix",
           Sigma = as.double(Sigma), p = as.integer(p), lamList = as.double(lambda.list),
           lamLength = as.integer(lambda.length), alpha = as.double(alpha),
           maxit = as.integer(maxit), error = as.double(error),
           precisionList = as.double(rep(0,p*p*lambda.length)), niterList = as.integer(rep(0,lambda.length)))
  dim(out$precisionList) = c(p,p,lambda.length)
  return(list(precisionList = out$precisionList, lambdaList = out$lamList, niterList=out$niterList))
}
step2estimator = function(Sigma, lambda.list = NULL, lambda.length = 100, alpha, weight = NULL,
                          maxit, error){
  p = dim(Sigma)[1]
  if(missing(weight)){
    weight = matrix(1,p,p)
    diag(weight) = 0
  }
  if(missing(lambda.list)){
    lambdaMax = lamMax_sgb(Sigma, alpha, weight)
    lambda.list = rev(seq(lambdaMax/lambda.length, lambdaMax, length.out = lambda.length))
  }
  else lambda.length = length(lambda.list)
  Sigma = as.vector(Sigma)
  weight = as.vector(weight)
  out = .C("precision",
           Sigma = as.double(Sigma), p = as.integer(p), lamList = as.double(lambda.list),
           lamLength = as.integer(lambda.length), alpha = as.double(alpha), weight = as.double(weight),
           maxit = as.integer(maxit), error = as.double(error),
           precisionList = as.double(rep(0,p*p*lambda.length)), niterList = as.integer(rep(0,lambda.length)))
  dim(out$precisionList) = c(p,p,lambda.length)
  dim(out$weight) = c(p,p)
  return(list(precisionList = out$precisionList, lambdaList = out$lamList,
              niterList=out$niterList, weight = out$weight))
}
genTheta = function(n, val){
  m = matrix(0, n, n)
  diag(m) = 1
  for(k in 1:length(val)){
    for(i in 1:(n-k)){
      m[i,i+k] = m[i+k,i] = val[k]
    }
  }
  return(m)
}

genTheta_random = function(p,prob, cond){
  B = matrix(rbinom(p*p,1,prob = sqrt(prob)), nrow = p, ncol = p)
  diag(B) = 0
  B = B*t(B)
  eig_out = eigen(B,symmetric = TRUE)
  dlt = (eig_out$values[1] - eig_out$values[p]*cond)/(cond-1)
  M = B/dlt + diag(p)
  return(M)
}
genTheta_random_block = function(p,prob, cond, n_blocks){
  M = matrix(0,nrow = p, ncol = p)
  p_div = as.integer(p/n_blocks)  
  for (i in 1:n_blocks){
    M[(i-1)*p_div+1:p_div,(i-1)*p_div+1:p_div] = genTheta_random(p_div,prob*n_blocks,p_div)
  }
  return(M)
}

lossDtrace = function(Sigma, Theta){
  r = -sum(diag(Theta)) + sum(diag(t(Theta)%*%Sigma%*%Theta))/2
  return(r)
}
lossGlasso = function(Sigma, Theta){
  return(-log(det(Theta)) + sum(diag(Sigma%*%Theta)))
}



# tpr and tnr
tpr_tnr = function(true, est,eps = 1e-7){
  p = dim(true)[1]
  tp = tn = 0
  np = nn = 0
  for(i in 1:p){
    for(j in 1:p){
      if(true[i,j] != 0){
        np = np + 1
        if(abs(est[i,j]) >eps) tp = tp + 1
      }
      else{
        nn = nn + 1
        if(abs(est[i,j])<=eps) tn = tn + 1
      }
    }
  }
  return(list(tpr = tp/np, tnr = tn/nn))
}
# Frobenius norm
f_norm = function(mat1, mat2){
  return(norm(mat1-mat2, type = c("f")))
}
# l2 norm (the largest singular value)
l2_norm = function(mat1, mat2){
  return(norm(mat1-mat2, type = c("2")))
}
# l1 norm
l1_norm = function(mat1, mat2){
  return(norm(mat1-mat2, type = c("1")))
}
# maximum absolute value
max_norm = function(mat1, mat2){
  return(norm(mat1-mat2, type = c("m")))
}
cv_matrixlist = function(precisionList, Sigma_ref, loss_function){
  len = dim(precisionList)[3]
  loss = rep(0, len)
  for(i in 1:len){
    if(loss_function == "Dtrace"){
      loss[i] = lossDtrace(Sigma_ref, precisionList[,,i])
    }
    if(loss_function == "Glasso"){
      loss[i] = lossGlasso(Sigma_ref, precisionList[,,i])
    }
  }
  return(precisionList[,,which.min(loss)])
}
weight_matrix = function(Thetahat, u){
  p = dim(Thetahat)[1]
  weight = matrix(0, p, p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      weight[i,j] = weight[j,i] = 1/(min(abs(Thetahat[i,j]), abs(Thetahat[j,i]))+u)
    }
  }
  return(weight)
}
find_min = function(precisionList, Sigma_ref, loss_function){
  len = dim(precisionList)[3]
  loss = rep(0, len)
  for(i in 1:len){
    if(loss_function == "Dtrace"){
      loss[i] = lossDtrace(Sigma_ref, precisionList[,,i])
    }
    if(loss_function == "Glasso"){
      loss[i] = lossGlasso(Sigma_ref, precisionList[,,i])
    }
  }
  return(list(id = which.min(loss), precision = precisionList[,,which.min(loss)]))
}
