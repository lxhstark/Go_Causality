# define a function that returns the power value given 
# ratio r1 = y1/sig1 and ratio rk = sigk/sig1
ratio_power = function(Y, r1, rk, N, K, compute_pcombine, compute_power, knum = 1, ktype = "homo", sig_general){
  # the treatment level y1
  y1 = Y[2]
  # derive the standard deviation sig array 
  sigk = NULL
  if (ktype == "homo"){
    sigk = rep(rk* y1/r1, K)
  }else{
    if (ktype == "normal"){
    ex_mean = log(y1 * rk/r1)
    sigk = exp(rnorm(K, ex_mean))
    }else{
    if (ktype == "min" | ktype == "max" | ktype == "median")
      {
     sigk[1] = sig_general[1]
     sigk[3:K] = sig_general[2:K-1]
    
     }else{
       print("Error: not choosing the existent type")
      }
     }
    }
  sigk[2] = y1/r1 

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

# minimal sigma case 
sig_stat_gen = function(mean, K, type = "min"){
  sig_general = exp(rnorm(K-1, log(mean)))
  sig_min_item = min(sig_general) 
  sig_max_item = max(sig_general)
  sig_median_item = median(sig_general)
  if (type == "min"){
    sig_item = sig_min_item
  }else{
    if (type == "max" ){
    sig_item = sig_max_item  
  }else{
   if (type == "median"){
    sig_item = sig_median_item
  }else{
    print("Error: not choosing the existent type")
      }
    }
  }
  return( list(item = sig_item,  sigk_rest = sig_general) )
}
