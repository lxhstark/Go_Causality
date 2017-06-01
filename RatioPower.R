# define a function that returns the power value given 
# ratio r1 = y1/sig1 and ratio rk = sigk/sig1
ratio_power = function(Y, r1, rk, N, K, compute_pcombine, compute_power, ktype = "homo", knum = 1){
  # the treatment level y1
  y1 = Y[2]
  # derive the standard deviation sig array 
  if (ktype == "homo"){
    sigk = rep(rk* y1/r1, K)
    sigk[2] = y1/r1 
  }else{
    ex_mean = log(y1 * rk/r1)
    sigk = exp(rnorm(K, ex_mean))
    sigk[2] = y1/r1
  }

  # compute the corresponding p_power 
  power = NULL
  # derive the critical value p_cv for null hypothesis
  p_cv = quantile(compute_pcombine(Y, sigk, N, K), 0.05)
  # get the sparse effect 
  Y_treat = Y
  Y_treat[2] = 1
  for (j in 1:knum){
    power[j] = compute_power(Y_treat, sigk, N, K, p_cv)
  }
  avg_power = mean(power)

  return(avg_power) 
}
