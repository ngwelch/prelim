library(MASS)
library(xtable)

##############################################################################
################################## mu prior ##################################
##############################################################################

get_mu_a = function(B, alpha=0.05){
  a = log(alpha/2)/log(B/(1+B))
  return(a)
}

# xmax = 160 for 25% sample; 120 for simulated 25% sample
get_mu_dtail = function(B, xmax=160, alpha=0.05){
  a = get_mu_a(B, alpha)
  
  tailProb = 0
  for(x in 0:(xmax-1)){
    l = lgamma(a+x) - lgamma(x+1) - lgamma(a) + a*log(B/(B+1)) + x*log(1/(B+1))
    el = exp(l)
    tailProb = tailProb+el
  }
  out = tailProb-(1-alpha/2)
  return(-out)
}

b = uniroot(get_mu_dtail, c(0.0001, 10), tol=.Machine$double.eps^0.5)
mu_shape = get_mu_a(b$root)
mu_rate = b$root*30
curve(dgamma(x, shape=mu_shape, rate=mu_rate), from=0.01, to=80)

q=c(0.025, 0.975)
replicatedQuantile = qgamma(q, shape=mu_shape, mu_rate)
replicatedMean = mu_shape/mu_rate
replicatedVar = mu_rate/(mu_rate**2)

brownQuantile = qgamma(q, shape=0.7, rate=0.004)
brownMean = 0.7/0.004
brownVar = 0.7/(0.004**2)

muPriorTable = data.frame(rbind(c(brownQuantile, brownMean, brownVar), 
                          c(replicatedQuantile, replicatedMean, replicatedVar)))
colnames(muPriorTable) = c("2.5%", "97.5%", "mean", "variance")
rownames(muPriorTable) = c("Brown et. al.", "Replicated")
tmp = xtable(muPriorTable, caption = "")
print(tmp)

##############################################################################
################################# theta prior ################################
##############################################################################

#qgamma(c(0.025, 0.975), shape=0.76, rate=0.0495)
source('~/prelim/src/R/simulation.R')
tmp = read.csv(file="/Users/nwelch/prelim/data/plantData.csv")[,c('x', 'y')]
tmp = sort(unique(tmp$x + 1.0i*tmp$y))
xi = tmp[Re(tmp)<=4.5 & Im(tmp)<=1]
ninedxixj = Mod(outer(xi, xi, FUN="-"))

testPlants = data.frame(xi, plant=1:9)
plot(testPlants$xi, pch="")
text(testPlants$xi, label=testPlants$plant)

########### Simulation when the infected plant can be any of the 9 ###########
it=1
alpha = seq(0.1, 1, 0.1)
beta_rate = seq(0.025, 0.5, 0.05)
sigma = seq(0.5, 1.5, 0.1)
coverage = matrix(0, nrow=length(alpha)*length(beta_rate)*length(sigma), ncol=4)
colnames(coverage) = c("alpha", "beta_rate", "sigma", "cover")
sampleSize=100

for(a in alpha){
  for(b in beta_rate){
    for(s in sigma){
      betweenDay1And480=0
      for(th in rgamma(sampleSize, shape=a, rate=b)){
        for(plant in 1:9){
          seed = rep(TRUE, nrow(ninedxixj))
          seed[plant] = FALSE
          infDay = scalarEpiSim(sigma=s, mu=1e-5, theta=th,
                                startTime=0, endTime=Inf, susceptableAtStart=seed, 
                                dxixj=ninedxixj, N=9, simCount=1, returnSummary=FALSE)*7
          nextInfDay = min(infDay[infDay>0], na.rm=TRUE)
          betweenDay1And480 = betweenDay1And480+(nextInfDay>=1 & nextInfDay<=480)
        }
      }
      coverage[it,] = c(a, b, s, betweenDay1And480/(9*sampleSize))
      it = it+1
    }
  }
}

coverage = as.data.frame(coverage)
coverage = coverage[order(coverage$cover, decreasing=TRUE),]
dci=0.005
tmp = coverage[coverage$cover<(0.95+dci) & coverage$cover>(0.95-dci),]
candidateMeans = tmp[order(tmp$cover, decreasing=TRUE),]
colMeans(candidateMeans)

###### Simulation when the infected plant can only be the center plant #######
it=1
nextit=1
alpha = seq(0.5, 0.95, 0.05)
beta_rate = seq(0.05, 0.2, 0.005)
sigma = seq(4, 5, 0.25)
coverage = matrix(0, nrow=length(alpha)*length(beta_rate)*length(sigma), ncol=4)
colnames(coverage) = c("alpha", "beta_rate", "sigma", "cover")
sampleSize=500

fixedSeed = rep(TRUE, nrow(ninedxixj))
fixedSeed[5] = FALSE

cat("Number of Evaluations to Complete:", nrow(coverage))
for(a in alpha){
  for(b in beta_rate){
    for(s in sigma){
      betweenDay1And480=0
      for(th in rgamma(sampleSize, shape=a, rate=b)){
        infDay = scalarEpiSim(sigma=s, mu=1e-5, theta=th,
                              startTime=0, endTime=Inf, dxixj=ninedxixj, 
                              susceptableAtStart=fixedSeed,
                              simCount=1, returnSummary=FALSE)*7
        nextInfDay = min(infDay[infDay>0], na.rm=TRUE)
        betweenDay1And480 = betweenDay1And480+(nextInfDay>=1 & nextInfDay<=480)
      }
      coverage[it,] = c(a, b, s, betweenDay1And480/sampleSize)
      if((it %% 20)==0){cat("Iteration: ", it, "\n")}
      it = it+1
    }
  }
}

coverage = as.data.frame(coverage)
coverage = coverage[order(coverage$cover, decreasing=TRUE),]
dci=0.005
tmp = coverage[coverage$cover<(0.95+dci) & coverage$cover>(0.95-dci),]
candidates = tmp[order(tmp$cover, decreasing=TRUE),]
colMeans(candidates)
candidates[candidates$beta_rate==min(candidates$beta_rate),]
colMeans(candidates[candidates$beta_rate==min(candidates$beta_rate),])
