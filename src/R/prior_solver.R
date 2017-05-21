library(MASS)
library(xtable)

# mu prior
get_mu_a = function(B, alpha=0.05){
  a = log(alpha/2)/log(B/(1+B))
  return(a)
}

get_mu_dtail = function(B, xmax=630, alpha=0.05){
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

q=c(0.025, 0.975)
replicatedQuantile = qgamma(q, shape=get_mu_a(b$root), rate=b$root)
replicatedMean = get_mu_a(b$root)/b$root
replicatedVar = get_mu_a(b$root)/(b$root**2)

brownQuantile = qgamma(q, shape=0.7, rate=0.004)
brownMean = 0.7/0.004
brownVar = 0.7/(0.004**2)

muPriorTable = data.frame(rbind(c(brownQuantile, brownMean, brownVar), 
                          c(replicatedQuantile, replicatedMean, replicatedVar)))
colnames(muPriorTable) = c("2.5%", "97.5%", "mean", "variance")
rownames(muPriorTable) = c("Brown et. al.", "Replicated")
tmp = xtable(muPriorTable, caption = "")
print(tmp)

# sigma prior
geta = function(a,b=10){
  3.9*(gamma(a) - pgamma(b*0.1,a, lower=F)*gamma(a))-gamma(a)+pgamma(b*50,a, lower=F)*gamma(a)
}

a=seq(0.1, 1,by=0.05)
b=seq(5,200,by=1)
N=10000
eps=1e-3
coverage = matrix(0, nrow=length(a)*length(b), ncol=4)
colnames(coverage) = c("a", "b", "cover", "abCoversIt")
it=1
for(i in a){
  for(j in b){
    smpl = rgamma(N, shape=i, scale=j)
    F500 = sum(smpl<500)/N
    F0.1 = sum(smpl<0.1)/N
    cover_ij = F500 - F0.1
    abCoversIt = F500<0.975 & F0.1<0.025
    #abCoversIt = ifelse(abs(cover_ij-0.95)<=eps, TRUE, FALSE)
    coverage[it,] = c(i, j, cover_ij, abCoversIt)
    it = it+1
  }
}
coverage = as.data.frame(coverage)
coverage = coverage[coverage$abCoversIt==TRUE,]
weakestSigmaPrior = coverage[which.max(coverage$a*(coverage$b**2)),]

# theta prior
source('~/prelim/src/R/simulation.R')
tmp = read.csv(file="/Users/nwelch/prelim/data/plantData.csv")[,c('x', 'y')]
tmp = sort(unique(tmp$x + 1.0i*tmp$y))
xi = tmp[Re(tmp)<=4.5 & Im(tmp)<=1]
dxixj = Mod(outer(xi, xi, FUN="-"))

it=1
theta = seq(0.05, 0.2, by=0.01)
#52.6225#seq(0.5, 2, by=0.05)
sigma = sort(rgamma(10, shape=weakestSigmaPrior$a, scale=weakestSigmaPrior$b))
sampleSize=5
coverage = matrix(0, nrow=length(theta)*length(sigma), ncol=3)
colnames(coverage) = c("theta", "sigma", "cover")
#out = array(dim=c(sampleSize,9,length(theta)))

for(th in theta){
  for(s in sigma){
    betweenDay1And448=0
    for(i in 1:sampleSize){
      for(j in 1:9){
        seed = rep(TRUE, nrow(dxixj))
        seed[j] = FALSE
        infDay = scalarEpiSim(sigma=s, mu=1e-5, theta=th,
                              startTime=0, endTime=Inf, susceptableAtStart=seed, 
                              dxixj=dxixj, N=9, simCount=1, returnSummary=FALSE)*7
        nextInfDay = min(infDay[infDay>0], na.rm=TRUE)
        betweenDay1And448 = betweenDay1And448+(nextInfDay>=1 & nextInfDay<=448)
      }
    }
    coverage[it,] = c(th, s, betweenDay1And448/(9*sampleSize))
    it = it+1
  }
}
coverage = as.data.frame(coverage)
candidateMeans = coverage[round(coverage$cover,2)==0.95,]
colMeans(candidateMeans)
