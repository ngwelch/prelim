# mu prior
get_a = function(B, alpha=0.05){
  a = log(alpha/2)/log(B/(1+B))
  return(a)
}

get_mu_dtail = function(B, xmax=630, alpha=0.05){
  a = get_a(B, alpha)
  p = 1/(1+B)
  oneMinusP = 1-p
  
  tailProb = 0
  for(x in 0:(xmax-1)){
    l = lgamma(a+x) - lgamma(x+1) - lgamma(a) + a*log(B/(B+1)) + x*log(1/(B+1))
    el = exp(l)
    tailProb = tailProb+el
  }
  out = tailProb-(1-alpha/2)
  return(-out)
}

b = uniroot(get_mu_dtail, c(0.0001, 10))


# theta prior
get_theta_dtail = function(B, xmax=1920, alpha=0.05){
  a = get_a(B, alpha)
  out = a*log(B) - a*log(B+xmax) - log(0.975)
  return(-out)
}

uniroot(get_theta_dtail, c(0.0001, 100))
