dst = read.csv(file="/Users/nwelch/prelim/data/dxixjTest.csv")


##############################################################################
########################## Epidemic Simulation ###############################
##############################################################################
### Simulated Epidemic Function
scalarEpiSim = function(sigma=1.1, mu=0.0028, theta=0.11,
                  startTime=0, endTime=30, susceptableAtStart=rep(TRUE,nrow(dst)), 
                  dxixj=dst, N=nrow(dst), simCount=1, returnSummary=TRUE){
  
  dstSq = dxixj**2
  twoSigmaSq = 2*(sigma**2)
  thetafxixj = (theta/(twoSigmaSq*pi))*exp(-dstSq/twoSigmaSq)
  
  weeks = c(0, seq(0.5, min(endTime-0.5, 69.5)))
  out = matrix(weeks, nrow=length(weeks), ncol=(simCount+1))

  # Force only one sim if the goal is to simulate infection times
  if(!returnSummary){simCount=1}
  for(s in 1:simCount){
    
    origInfTimes = Inf
    if(sum(susceptableAtStart==FALSE)>0){
      ratesAtInfectedPlants = thetafxixj[,!susceptableAtStart, drop=FALSE]
      sampleNewInfTimes = apply(ratesAtInfectedPlants, 2, function(r) rexp(N, r))
      origInfTimes = startTime + sampleNewInfTimes
    }
    
    infTime = pmin(startTime + rexp(N, mu), origInfTimes, endTime)
    susceptibles = which((infTime < endTime) & susceptableAtStart)
    
    while(length(susceptibles)){
      
      currentTime = min(infTime[susceptibles])
      currentPlant = susceptibles[which.min(infTime[susceptibles])]
      currentPlantRates = thetafxixj[,currentPlant]
      
      nextInfectionTime = rexp(N, currentPlantRates)
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

oneSim = scalarEpiSim(mu=0.004274672,
                      sigma=1.229166,
                      theta=0.1180905,
                      simCount=1000)
hist(colSums(oneSim[,-1]), breaks=30)
timeSim = scalarEpiSim(returnSummary = FALSE)

### Sampler that runs over a vector of parameter values
vectorEpiSim = function(muSigmaThetaVector, startTime=0, endTime=30, 
                        susceptableAtStart=rep(TRUE,nrow(dst)), dxixj=dst, 
                        N=nrow(dst), simCount=1, returnSummary=TRUE){
  
  if(returnSummary){
    result = NULL
    for(s in 1:nrow(muSigmaThetaVector)){
      mu = muSigmaThetaVector[s,1]
      sigma = muSigmaThetaVector[s,2]
      theta = muSigmaThetaVector[s,3]
      sim = scalarEpiSim(sigma, mu, theta, startTime, endTime, susceptableAtStart,
                         dxixj, N, simCount, returnSummary=TRUE)
      
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
      sim = scalarEpiSim(sigma, mu, theta, startTime, endTime, susceptableAtStart,
                         dxixj, N, simCount=1, returnSummary=FALSE)
      out=cbind(out, sim)
    }
    colnames(out) = paste0("sim", 1:nrow(muSigmaThetaVector))
  }
  return(data.frame(out))
}

u = c(0.0028, 0.0025, 0.003)
s = c(1.1, 1.15, 1.08)
t = c(0.11, 0.08, 0.15)
vecSim = vectorEpiSim(muSigmaThetaVector=cbind(u,s,t), simCount=100)


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

sampleTau = read.csv(file="~/prelim/data/mcmc_chain_julia.csv")[,1:34]
pstSim = postEpiSim(sampleTau)

par(mfcol=c(1,2), mar=c(3.1, 3.1, 1, 1), mgp=c(2,0.5,0))
matplot(x=pstSim['weeks',], y=t(pstSim[5:200,]), type='l', lty=3, col="#50505030",
        ylab="infections per week", xlab="week")
matlines(x=pstSim['weeks',], y=t(pstSim[2:4,]), col='black', lty=c(1,2,2), lwd=3)

matplot(x=vecSim$week, y=vecSim[,5:ncol(vecSim)], type='l', lty=3, col="#50505030",
        ylab="infections per week", xlab="week")
matlines(x=vecSim$week, y=vecSim[,2:4], col='black', lty=c(1,2,2), lwd=3)
par(mfcol=c(1,1), mar=c(5,4,4,2)+0.1, mgp=c(3,1,0))

##############################################################################
############################ Density Plots ###################################
##############################################################################

par(mfcol=c(1,2), mar=c(3.1, 3.1, 1, 1), mgp=c(2,0.5,0))
### First time period
breakPoints = seq(10,14,len=10)
tmp = apply(sampleTau[,1:7], 2, hist, breaks=breakPoints, plot=F)
postSampleDensity = NULL
for(psd in 1:length(tmp))
  postSampleDensity = cbind(postSampleDensity, tmp[[psd]]$density)

matplot(breakPoints[-1]-diff(breakPoints)/2, postSampleDensity, 
        type="l", lty=rep(1:3, c(3,3,3)),
        col=rep(grey(c(0, 0.4, 0.6, 0.8)),3),
        xlab="week", ylab="density")

### Last time period
breakPoints = seq(23,30,len=10)
tmp = apply(sampleTau[,25:34], 2, hist, breaks=breakPoints, plot=F)
postSampleDensity = NULL
for(psd in 1:length(tmp))
  postSampleDensity = cbind(postSampleDensity, tmp[[psd]]$density)

postSampleMean = apply(postSampleDensity, 2, mean)
postSampleMean = matrix(postSampleMean, nrow(postSampleDensity), 
                        ncol(postSampleDensity))
centering = postSampleDensity - postSampleMean
centeringL = apply(centering[1:3,], 2, sum)
centeringH = apply(centering[7:9,], 2, sum)

plotLines=c(order(centeringL, decreasing=T)[c(1,10)],
            order(centeringH, decreasing=T)[c(1,10)])

matplot(breakPoints[-1]-diff(breakPoints)/2, postSampleDensity, 
        type="l",
        col="#50505030",
        xlab="week", ylab="density")
matlines(breakPoints[-1]-diff(breakPoints)/2, 
         postSampleDensity[,toplot], 
         col="white", lty=1, lwd=4)
matlines(breakPoints[-1]-diff(breakPoints)/2, 
         postSampleDensity[,toplot], 
         col="black",lty=c(1,1,2,2))
par(mfcol=c(1,1), mar=c(5,4,4,2)+0.1, mgp=c(3,1,0))

##############################################################################
############################# Prediction #####################################
##############################################################################

t6 = read.csv(file="/Users/nwelch/prelim/data/plantData.csv")$t6
susceptable = ifelse(t6==0, TRUE, FALSE)
onePredSim = scalarEpiSim(startTime=30, endTime=Inf, 
                          susceptableAtStart=susceptable,
                          returnSummary=FALSE)
testMatrix = cbind(u,s,t)
for(i in 1:4){ testMatrix = rbind(testMatrix, testMatrix)}
manyPredSim = vectorEpiSim(muSigmaThetaVector=testMatrix,
                           startTime=30, endTime=Inf, 
                           susceptableAtStart=susceptable,
                           returnSummary=FALSE)

cumTimes = c(35, 40, 50, 60)
colBreaks = c(-0.1, 0.2, 0.8, 0.9999, 1.01)
colLabel = c('zeroTo19', 'twentyTo79', 'eightyTo99', 'knownInfection')
probInfected = probInfLabel = NULL
for(ct in cumTimes){
  probabilities = apply(manyPredSim, 1, function(t) mean(t<ct))
  probInfected = cbind(probInfected, probabilities)
  probLabels = as.character(cut(probabilities, breaks=colBreaks, labels=colLabel))
  probInfLabel = cbind(probInfLabel, probLabels)
}
colnames(probInfected) = colnames(probInfLabel) = paste0("week", cumTimes)
write.csv(probInfected, file="/Users/nwelch/prelim/data/fig5Probabilities.csv", 
          row.names=FALSE)
write.csv(probInfLabel, file="/Users/nwelch/prelim/data/fig5ProbLabels.csv", 
          row.names=FALSE)
