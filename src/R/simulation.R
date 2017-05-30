##############################################################################
########################## Epidemic Simulation ###############################
##############################################################################

### Simulated Epidemic Function
scalarEpiSim = function(mu=0.0028, sigma=1.1, theta=0.11,
                  startTime=0, endTime=30, dxixj, susceptableAtStart=NULL, 
                  simCount=1, returnSummary=TRUE){
  
  N=nrow(dxixj)
  if(length(susceptableAtStart)==0){
    susceptableAtStart=rep(TRUE,nrow(dxixj))
  }
  
  dstSq = dxixj**2
  twoSigmaSq = 2*(sigma**2)
  thetafxixj = (theta/(twoSigmaSq*pi))*exp(-dstSq/twoSigmaSq)
  
  weeks = c(0, seq(0.5, min(endTime-0.5, 69.5)))
  out = matrix(weeks, nrow=length(weeks), ncol=(simCount+1))

  # Force only one sim if the goal is to simulate infection times
  if(!returnSummary){simCount=1}
  getRates = function(r){
    rates = suppressWarnings(rexp(N, r))
    rates
  } 
  for(s in 1:simCount){
    
    origInfTimes = Inf
    if(sum(susceptableAtStart==FALSE)>0){
      ratesAtInfectedPlants = thetafxixj[,!susceptableAtStart, drop=FALSE]
      ratesAtInfectedPlants[ratesAtInfectedPlants<.Machine$double.xmin]=2*.Machine$double.xmin
      sampleNewInfTimes = apply(ratesAtInfectedPlants, 2, getRates)
      origInfTimes = startTime + sampleNewInfTimes
    }
    
    timeToNextInfection = cbind(startTime + rexp(N, mu), origInfTimes, endTime)
    infTime = apply(timeToNextInfection, 1, min)
    susceptibles = which((infTime < endTime) & susceptableAtStart)
    
    while(length(susceptibles)){
      
      currentTime = min(infTime[susceptibles])
      currentPlant = susceptibles[which.min(infTime[susceptibles])]
      currentPlantRates = thetafxixj[,currentPlant]
      
      nextInfectionTime = suppressWarnings(rexp(N, currentPlantRates))
      nextInfectionTime[is.na(nextInfectionTime)] = Inf
      infTime = pmin(infTime, currentTime + nextInfectionTime)
      susceptibles = which((infTime < endTime) & (infTime > currentTime) & 
                             susceptableAtStart)
    }
    infTime = ifelse(infTime>=endTime, Inf, infTime)
    infTime[!susceptableAtStart] = startTime
    
    # create an index and count up number of infections per index interval
    results=c()
    for(w in weeks){
      results = c(results, sum(abs(infTime-w)<0.5))
    }
   out[,s+1] = results 
  }
  
  if(returnSummary){
    colnames(out) = c('week', paste0('sim', 1:simCount))
    return(data.frame(out))
  } else {
    return(infTime)
  }
}


### Sampler that runs over a vector of parameter values
vectorEpiSim = function(muSigmaThetaVector, startTime=0, endTime=30, 
                        susceptableAtStart=NULL, dxixj, 
                        simCount=1, returnSummary=TRUE){
  
  N=nrow(dxixj)
  if(length(susceptableAtStart)==0){
    susceptableAtStart=rep(TRUE,nrow(dxixj))
  }
  
  if(returnSummary){
    result = NULL
    for(s in 1:nrow(muSigmaThetaVector)){
      mu = muSigmaThetaVector[s,1]
      sigma = muSigmaThetaVector[s,2]
      theta = muSigmaThetaVector[s,3]
      sim = scalarEpiSim(mu, sigma, theta, startTime, endTime, dxixj,
                         susceptableAtStart, simCount, returnSummary=TRUE)
      
      if(s>1){
        result = cbind(result, sim[,-1, drop=FALSE])
      } else {
        result = sim
      }
    }
    tmp = apply(result[,-1], 1, function(w) c(mean(w), quantile(w, c(0.025, 0.975))))
    sampleCI = t(tmp)
    colnames(sampleCI)[1] = "mean"
    
    out = cbind(result[,1,drop=FALSE], sampleCI, result[,-1,drop=FALSE])
    
  } else {
    out = NULL
    for(s in 1:nrow(muSigmaThetaVector)){
      mu = muSigmaThetaVector[s,1]
      sigma = muSigmaThetaVector[s,2]
      theta = muSigmaThetaVector[s,3]
      sim = scalarEpiSim(mu, sigma, theta, startTime, endTime, dxixj,  
                         susceptableAtStart, simCount=1, returnSummary=FALSE)
      out=cbind(out, sim)
    }
    colnames(out) = paste0("sim", 1:nrow(muSigmaThetaVector))
  }
  return(data.frame(out))
}

##############################################################################
########################## Posterior Sampler #################################
##############################################################################

postEpiSim = function(tau){
  
  # create an index and count up number of infections per index interval
  weeks = c(0, seq(0.5, 29.5))
  post = NULL
  for(w in weeks){
    post = cbind(post, apply(tau, 1, function(t) sum(abs(t-w)<0.5)))
  }
  rownames(post) = paste0('sample', 1:nrow(tau))
  colnames(post) = paste0('week', weeks)
  
  # find the CI for each of sample
  postCI = apply(post, 2, function(w) c(mean(w), quantile(w, c(0.025, 0.975))))
  rownames(postCI)[1] = "mean"
  
  # return both results in the same matrix
  rbind(weeks, postCI, post)
}
