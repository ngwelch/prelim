getLogLambda = function(tau, theta, mu, sigma, dst, Tlast=30){
  
  iLessTlast = which(tau < Tlast)
  sumLog = 0
  for(i in iLessTlast){
    sumTauJLessTauI = 0
    jLessI = which(tau<tau[i])
    for(j in jLessI){
      sumTauJLessTauI = sumTauJLessTauI + dnorm(dst[i,j], mean=mu, sd=sigma)
    }
    sumLog = sumLog + log(mu + theta*sumTauJLessTauI)
  }
  return(sumLog)
}

###############################################################################
# mu Likelihood Function
###############################################################################


llRatio_mu = function(mu, muStar, theta, sigma, tau, dst, Tlast=30, 
                      mu_a=0.7, mu_b=0.004){
  
  tauLessTIndex = which(tau < Tlast)
  
  dmu = mu - muStar
  
  dllTerm1 = sum(tau[tauLessTIndex]*dmu)
  dllTerm2 = length(tauLessTIndex) * Tlast * dmu
  
  dllTerm3a = getLogLambda(tau, theta, mu, sigma, dst, Tlast)
  dllTerm3b = getLogLambda(tau, theta, mu=muStar, sigma, dst, Tlast)
  dllTerm3 = dllTerm3a - dllTerm3b
  
  prior = dgamma(mu, shape=mu_a, rate=mu_b)
  priorStar = dgamma(muStar, shape=mu_a, rate=mu_b)
  
  dll = dllTerm1 + dllTerm2 - dllTerm3 + prior - priorStar

  return(dll)
}#end update_mu



###############################################################################
# sigma Likelihood Function
###############################################################################


llRatio_sigma = function(sigma, sigmaStar, mu, tau, theta, dst, Tlast=30, 
                         sigma_a=0.5, sigma_b=100){
  
  dllTerm1 = getSumLessT_sigma(sigma, sigmaStar, mu, tau, theta, dst)
  
  dllTerm2 = getSumGreaterT_sigma(sigma, sigmaStar, mu, tau, theta, dst)
  
  dllTerm3 = getLogLambda(tau, theta, mu, sigmaStar, dst, Tlast) - 
    getLogLambda(tau, theta, mu, sigma, dst, Tlast)
  
  priorStar = dgamma(sigmaStar, shape=sigma_a, rate=sigma_b)
  prior = dgamma(sigma, shape=sigma_a, rate=sigma_b)
  
  dll = dllTerm1 + dllTerm2 + dllTerm3 + priorStar - prior
  return(dll)
}

# llRatio_sigma helper function 1

getSumLessT_sigma = function(sigma, sigmaStar, mu, tau, theta, dst, Tlast=30){
  
  tauLessTIndex = which(tau <= Tlast)
  
  sumLessT = 0
  for(i in tauLessTIndex){
    sum_ij = 0
    tauJLessTauIIndex = which( tau < tau[i])
    for(j in tauJLessTauIIndex){
      sum_ij = sum_ij + (tau[i] - tau[j])*(
        dnorm(dst[i,j], mean=mu, sd=sigma) - 
          dnorm(dst[i,j], mean=mu, sd=sigmaStar))
    }
    sumLessT = sumLessT + theta*sum_ij
  }
  return(sumLessT)
}

# llRatio_sigma helper function 2

getSumGreaterT_sigma = function(sigma, sigmaStar, mu, tau, theta, dst, Tlast=30){
  
  tauLessTIndex = which(tau <= Tlast)
  tauGreaterTIndex = which(tau > Tlast)
  
  sumGreaterT = 0
  for(i in tauGreaterTIndex){
    sum_ij = 0
    for(j in tauLessTIndex){
      sum_ij = sum_ij + (Tlast - tau[j])*(
        dnorm(dst[i,j], mean=mu, sd=sigma) - 
          dnorm(dst[i,j], mean=mu, sd=sigmaStar))
    }
    sumGreaterT = sumGreaterT + theta*sum_ij
  }
  return(sumGreaterT)
}

###############################################################################
# theta Likelihood Function
###############################################################################


llRatio_theta = function(theta, thetaStar, mu, tau, sigma, dst, 
                         Tlast=30, theta_a=0.8, theta_b=10){
  
  dllTerm1 = getSumLessT_theta(theta, thetaStar, mu, tau, sigma, dst, Tlast)
  
  dllTerm2 = getSumGreaterT_theta(theta, thetaStar, mu, tau, sigma, dst, Tlast)
  
  dllTerm3 = getLogLambda(tau, thetaStar, mu, sigma, dst, Tlast) - 
    getLogLambda(tau, theta, mu, sigma, dst, Tlast)
  
  priorStar = dgamma(thetaStar, shape=theta_a, rate=theta_b)
  prior = dgamma(theta, shape=theta_a, rate=theta_b)
  
  dll = dllTerm1 + dllTerm2 + dllTerm3 + priorStar - prior
  return(dll)
}


# llRatio_theta helper function 1

getSumLessT_theta = function(theta, thetaStar, mu, tau, sigma, dst, Tlast=30){
  
  tauLessTIndex = which(tau <= Tlast)
  
  sumLessT = 0
  for(i in tauLessTIndex){
    sum_ij = 0
    tauJLessTauIIndex = which( tau < tau[i])
    for(j in tauJLessTauIIndex){
      sum_ij = sum_ij + (tau[i] - tau[j])*(theta - thetaStar)*
        dnorm(dst[i,j], mean=mu, sd=sigma)
    }
    sumLessT = sumLessT + theta*sum_ij
  }
  return(sumLessT)
}

# llRatio_theta helper function 2

getSumGreaterT_theta = function(theta, thetaStar, mu, tau, sigma, dst, Tlast=30){
  
  tauLessTIndex = which(tau <= Tlast)
  tauGreaterTIndex = which(tau > Tlast)
  
  sumGreaterT = 0
  for(i in tauGreaterTIndex){
    sum_ij = 0
    for(j in tauLessTIndex){
      sum_ij = sum_ij + (Tlast - tau[j])*(theta - thetaStar)*
        dnorm(dst[i,j], mean=mu, sd=sigma)
    }
    sumGreaterT = sumGreaterT + theta*sum_ij
  }
  return(sumGreaterT)
}


###############################################################################
# tau Likelihood Function
###############################################################################


llRatio_tau = function(tau, tauiStar, iStar, mu, theta, sigma, dst, Tlast=30){
  
  N = length(tau)
  
  dllTerm1 = mu*(tau[iStar] - tauiStar)
  
  dllTerm2 = getMinSum(iStar, tauiStar, mu, tau, theta, sigma, dst)
  
  dllTerm3 = getMaxSum(iStar, tauiStar, mu, tau, theta, sigma, dst)
  
  dllTerm4 = getSumIntervalBelowTauStar(iStar, tauiStar, tau, theta, 
  										sigma, dst)
  
  dllTerm5 = getSumIntervalAboveTauStar(iStar, tauiStar, tau, theta, 
  										sigma, dst)
  
  dllTerm6 = getLogLambda(tau, theta, mu, sigma, dst, Tlast)
  
  tauStar = tau #c(tau[0:(iStar-1)], tauiStar, tau[(iStar+1):N])
  tauStar[iStar] = tauiStar
  dllTerm7 = getLogLambda(tauStar, theta, mu, sigma, dst, Tlast)
  
  dll = dllTerm1 + dllTerm2 + dllTerm3 + dllTerm4 + dllTerm5 + 
  	dllTerm6 - dllTerm7 
  return(dll)
}


# llRatio_tau helper function 1

getMinSum = function(iStar, tauiStar, mu, tau, theta, sigma, dst){
  
  ti = tau[iStar]
  minTau = min(tauiStar, ti)
  
  jLessMin = which(tau < minTau)
  minSum = 0
  for(j in jLessMin){
    minSum = minSum + dnorm(dst[iStar,j], mean=mu, sd=sigma)
  }
  minSum = (ti - tauiStar)*theta*minSum
  return(minSum)
}

# llRatio_tau helper function 2

getMaxSum = function(iStar, tauiStar, mu, tau, theta, sigma, dst){
  
  ti = tau[iStar]
  maxTau = max(tauiStar, ti)
  
  jGreaterMax = which(tau > maxTau)
  maxSum = 0
  for(j in jGreaterMax){
    maxSum = maxSum + dnorm(dst[iStar,j], mean=mu, sd=sigma)
  }
  maxSum = (tauiStar - ti)*theta*maxSum
  return(maxSum)
}

# llRatio_tau helper function 3

getSumIntervalBelowTauStar = function(iStar, tauiStar, mu, tau, theta, sigma, dst){
  
  ti = tau[iStar]
  
  jInBetween = which((tau > ti) & (tau < tauiStar))
  sumBelow = 0
  for(j in jInBetween){
    sumBelow = sumBelow + 
      (2*tau[j] - ti - tauiStar)*dnorm(dst[iStar,j], mean=mu, sd=sigma)
  }
  sumBelow = theta*sumBelow
  return(sumBelow)
}

# llRatio_tau helper function 4

getSumIntervalAboveTauStar = function(iStar, tauiStar, mu, tau, theta, sigma, dst){
  
  ti = tau[iStar]
  
  jInBetween = which((tau > tauiStar) & (tau < ti))
  sumAbove = 0
  for(j in jInBetween){
    sumAbove = sumAbove + 
      (ti + tauiStar - 2*tau[j])*dnorm(dst[iStar,j], mean=mu, sd=sigma)
  }
  sumAbove = theta*sumAbove
  return(sumAbove)
}

