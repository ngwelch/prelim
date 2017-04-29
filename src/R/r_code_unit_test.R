library(RUnit)
source('~/prelim/src/R/likelihood_functions.R')

testdata = read.csv(file="/Users/nwelch/prelim/data/infectionDataTest.csv")
testdist = read.csv(file="/Users/nwelch/prelim/data/infectedDistancesTest.csv")

#priors
mu_a=0.7; mu_b=0.004
theta_a=0.8; theta_b=10
sigma_a=0.5; sigma_b=100

# initial/test values
u = mu_a*mu_b
th = theta_a*theta_b
s = sigma_a*sigma_b
t = c(14, 10, 6)
tl = max(t)

###############################################################################
# getLogLambda Tests
###############################################################################
gll = getLogLambda(tau=t, theta=th, mu=u, sigma=s, dst=testdist, Tlast=tl)
  
###############################################################################
# llRatio_mu Tests
###############################################################################
mu_increase = llRatio_mu(mu=u, muStar=0.01, theta=th, sigma=s, tau=t, Tlast=tl,
                         dst=testdist)
mu_decrease = llRatio_mu(mu=u, muStar=0.001, theta=th, sigma=s, tau=t, Tlast=tl,
                         dst=testdist)

###############################################################################
# llRatio_theta Tests
###############################################################################
theta_helper_1_increase = getSumLessT_theta(theta=th, thetaStar=20, mu=u, 
                                            tau=t, sigma=s, dst=testdist, 
                                            Tlast=tl)

theta_helper_1_decrease = getSumLessT_theta(theta=th, thetaStar=1, mu=u, 
                                            tau=t, sigma=s, dst=testdist, 
                                            Tlast=tl)

theta_helper_2_increase = getSumGreaterT_theta(theta=th, thetaStar=20, mu=u,
                                               tau=t, sigma=s, dst=testdist,
                                               Tlast=10)

theta_helper_2_decrease = getSumGreaterT_theta(theta=th, thetaStar=1, mu=u,
                                               tau=t, sigma=s, dst=testdist,
                                               Tlast=10)

theta_increase = llRatio_theta(theta=th, thetaStar=20, mu=u, tau=t, 
                               sigma=s, dst=testdist, Tlast=tl)
theta_decrease = llRatio_theta(theta=th, thetaStar=1, mu=u, tau=t, 
                               sigma=s, dst=testdist, Tlast=tl)
###############################################################################
# llRatio_sigma Tests
###############################################################################
sigma_helper_1_increase = getSumLessT_sigma(sigma=s, sigmaStar=100, mu=u, 
                                            tau=t, theta=th, dst=testdist, 
                                            Tlast=tl)
sigma_helper_1_decrease = getSumLessT_sigma(sigma=s, sigmaStar=1, mu=u, 
                                            tau=t, theta=th, dst=testdist, 
                                            Tlast=tl)
sigma_helper_2_increase = getSumGreaterT_sigma(sigma=s, sigmaStar=100, mu=u, 
                                               tau=t, theta=th, dst=testdist, 
                                               Tlast=10)
sigma_helper_2_decrease = getSumGreaterT_sigma(sigma=s, sigmaStar=1, mu=u, 
                                               tau=t, theta=th, dst=testdist, 
                                               Tlast=10)

sigma_increase = llRatio_sigma(sigma=s, sigmaStar=100, mu=u, tau=t, theta=th, 
                               dst=testdist, Tlast=tl)
sigma_decrease = llRatio_sigma(sigma=s, sigmaStar=1, mu=u, tau=t, theta=th, 
                               dst=testdist, Tlast=tl)
###############################################################################
# llRatio_tau Tests
###############################################################################
tau_helper_1_increase = getMinSum(iStar=1, tauiStar=t[1]+0.5, mu=u, tau=t, 
                                  theta=th, sigma=s, dst=testdist)
tau_helper_1_decrease = getMinSum(iStar=1, tauiStar=t[1]-0.5, mu=u, tau=t, 
                                  theta=th, sigma=s, dst=testdist)
  
tau_helper_2_increase = getMaxSum(iStar=2, tauiStar=t[2]+0.5, mu=u, tau=t, 
                                  theta=th, sigma=s, dst=testdist)
tau_helper_2_decrease = getMaxSum(iStar=2, tauiStar=t[2]-0.5, mu=u, tau=t, 
                                  theta=th, sigma=s, dst=testdist)

tau_helper_3_increase = getSumIntervalBelowTauStar(iStar=3, tauiStar=18, mu=u,
                                                   tau=t, theta=th, sigma=s, 
                                                   dst=testdist)
tau_helper_3_increase2 = getSumIntervalBelowTauStar(iStar=2, tauiStar=15, 
                                                    mu=u, tau=t, theta=th, 
                                                    sigma=s, dst=testdist)

tau_helper_3_decrease = getSumIntervalBelowTauStar(iStar=3, tauiStar=12, 
                                                   mu=u, tau=t, theta=th, 
                                                   sigma=s, dst=testdist)
tau_helper_3_decrease2 = getSumIntervalBelowTauStar(iStar=3, tauiStar=5, 
                                                   mu=u, tau=t, theta=th, 
                                                   sigma=s, dst=testdist)
tau_helper_3_decrease3 = getSumIntervalBelowTauStar(iStar=3, tauiStar=9, 
                                                    mu=u, tau=t, theta=th, 
                                                    sigma=s, dst=testdist)



tau_helper_4_decrease = getSumIntervalAboveTauStar(iStar=1, tauiStar=5, 
                                                   mu=u, tau=t, theta=th, 
                                                   sigma=s, dst=testdist)
tau_helper_4_decrease2 = getSumIntervalAboveTauStar(iStar=1, tauiStar=8, 
                                                   mu=u, tau=t, theta=th, 
                                                   sigma=s, dst=testdist)
tau_increase = llRatio_tau(tau=t, tauiStar=t[1]+0.5, iStar=1, mu=u, 
                           theta=th, sigma=s, dst=testdist, Tlast=tl)
tau_decrease = llRatio_tau(tau=t, tauiStar=t[1]-0.5, iStar=1, mu=u, 
                           theta=th, sigma=s, dst=testdist, Tlast=tl)
