library(dplyr)
library(MASS)
#library(profr)

source('~/prelim/src/R/likelihood_functions.R')

# TEST!
df = read.csv(file="/Users/nwelch/prelim/data/infectionDataTest.csv")
dst = read.csv(file="/Users/nwelch/prelim/data/infectedDistancesTest.csv")

trials = 100
infectionCount = sum(df$t6)
chain = matrix(NA, nrow=trials, ncol=(3+infectionCount))

#priors
mu_a=0.7; mu_b=0.004
theta_a=0.8; theta_b=10
sigma_a=0.5; sigma_b=100

# initial values
theta = theta_a*theta_b #0.105
sigma = sigma_a*sigma_b #1.15
mu = mu_a*mu_b #0.003
tau = t(df$tau)
chain[1,] = c(tau, theta, sigma, mu)
colnames(chain) = c(paste0("t", 1:infectionCount), "theta", "sigma", "mu")

accept_count = c(rep(0, infectionCount+3), trials, rep(0,3))
names(accept_count) = c(colnames(chain), 'trials', 'user', 'system', 'elapsed')

t = system.time(
for(r in 2:trials){
  
  tau = chain[r-1, 1:infectionCount]
  theta = chain[r-1, "theta"]
  mu = chain[r-1, "mu"]
  sigma = chain[r-1, "sigma"]
  
  logU = log( runif(infectionCount+3) )
  
  thetaStar = rnorm(1, mean=theta, sd=0.005)
  if(thetaStar>0){
    logLRatio_theta = llRatio_theta(theta=theta, thetaStar=thetaStar, 
                                    mu=mu, tau=tau, sigma=sigma,
                                    dst=dst)
    
    if(logU[1] < logLRatio_theta){
      theta = thetaStar
      accept_count['theta'] = accept_count['theta'] + 1
    }
  }
  
  muStar = rnorm(1, mean=mu, sd=0.0005)
  if(muStar>0){
    logLRatio_mu = llRatio_mu(mu=mu, muStar=muStar, theta=theta, sigma=sigma, 
                              tau=tau, dst=dst)
    
    if(logU[2] < logLRatio_mu){
      mu = muStar
      accept_count['mu'] = accept_count['mu'] + 1
    }
  }
  
  sigmaStar = rnorm(1, mean=sigma, sd=0.05)
  if(sigmaStar>0){
    logLRatio_sigma = llRatio_sigma(sigma=sigma, sigmaStar=sigmaStar, 
    								mu=mu, tau=tau, theta=theta, dst=dst)
    
    if(logU[3] < logLRatio_sigma){
      sigma = sigmaStar
      accept_count['sigma'] = accept_count['sigma'] + 1
    }
  }
  
  for(i in 1:infectionCount){
    iStar=i
    tauiStar = rnorm(1, mean=tau[i], sd=1)
    lowerBound = df$tauLowerBound[i]
    upperBound = df$tauUpperBound[i]
    
    if((tauiStar<upperBound) & (tauiStar>lowerBound)){
      logLRatio_tau = llRatio_tau(tau=tau, tauiStar=tauiStar, iStar=i, 
      							  mu=mu, theta=theta, sigma=sigma, dst=dst)
      
      if(logU[(i+3)] < logLRatio_tau){
        tau[i] = tauiStar
        accept_count[i] = accept_count[i] + 1
      }
    }
  }
  
  chain[r,] = c(tau, theta, sigma, mu)
}
)

chain = chain[complete.cases(chain),]
write.table(chain, file="~/prelim/data/mcmc_chain_R.txt")

accept_count[c('user', 'system', 'elapsed')] = summary(t)
write.table(accept_count, file="~/prelim/data/mcmc_performance_summary_R.txt")
