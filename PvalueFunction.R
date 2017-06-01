#Fundamental Simulation Function
##Generate p-value function
compute_pvalue = function(Yk, sigk, N, K){
  n <- N/K
  nk = rep(n, K) # assume evenly distributed number of subgroups nk
  # initialize the parameters
  tao = NULL # causal estimand
  Yob = NULL
  sig = NULL
  V = NULL # neyman estimators of variance
  t = NULL # t-statiscs
  # p-value matrix
  p = matrix(0, ncol = K, nrow = n)
  # gen YK obs 
  for (j in 1:n){
    for (i in 1:K){ 
      # derive the sample mean and sample variance 
      Yob[i] = rnorm(1, Yk[i], sigk[i] * sqrt(K/N))
      sig[i] = sigk[i] * rchisq(1, nk[i] - 1) / (nk[i] - 1)
      Y0 = Yob[1] # non-treat mean
      tao[i] = Yob[i] - Y0 
      V[i] = sig[i]/nk[i] + sig[1]/nk[1]
      # calculate the t-statistics
      t[i] = (tao[i] - 0) / sqrt(V[i])
      # get the two-sided p-value
      p[j,i] = (1 - pt(abs(t[i]), nk[i] - 2)) * 2
    }
  }
  return(p)
}


##Generate p-combine function
compute_pcombine <- function(Yk, sigk, N, K){
  n <- N/K
  p_comb_val <- NULL
  pmatrix <- compute_pvalue(Yk, sigk, N, K)
  for (k in 1:n){
  p_comb_val[k] = min(pmatrix[k,-1]) #minimum
  }
  return (p_comb_val)
}


##Generate power function
compute_power <- function(Yk, sigk, N, K, p_cv){
  n <- N/K
  power = 0
  pcomb = compute_pcombine(Yk, sigk, N, K)
  rejec = 0
  rec = 0
  for (i in 1:n){
    if(pcomb[i] <= p_cv){
      rejec = rejec + 1
    }else{
      rec = rec + 1
    }
  }
  power = rejec / (rejec + rec)
}




