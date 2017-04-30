using DataFrames;
using Distributions;

# priors
mu_a=0.7; 
mu_b=0.004;
theta_a=0.8; 
theta_b=10;
sigma_a=0.5; 
sigma_b=100;

##############################################################################
# log-likelihood
##############################################################################

function getLogLambda(tau, theta, mu, sigma, dst, Tlast=30)
    iLessTlast = find(x-> x<=Tlast, tau)
    gaussian = Normal(mu, sigma)
    
    sumLog = 0
    for i = iLessTlast
        sumTauJLessTauI = 0
        jLessI = find(x-> x < tau[i], tau)

        for j = jLessI
            sumTauJLessTauI = sumTauJLessTauI + pdf(gaussian, dst[i,j])
        end
        sumLog = sumLog + log(mu + theta*sumTauJLessTauI)
    end
    sumLog
end

##############################################################################
# mu update
##############################################################################

muPrior = Gamma(mu_a, mu_b)

function llRatio_mu(mu, muStar, theta, sigma, tau, dst; 
                        Tlast=30, mu_a=0.7, mu_b=0.004)
    tauLessTIndex = find(x-> x<=Tlast, tau)
    tauGreaterTIndex = find(x-> x>Tlast, tau)
    
    dmu = mu - muStar
    
    dllTerm1 = sum(tau[tauLessTIndex]*dmu)
    dllTerm2 = length(tauGreaterTIndex) * Tlast * dmu
    
    dllTerm3a = getLogLambda(tau, theta, muStar, sigma, dst, Tlast)
    dllTerm3b = getLogLambda(tau, theta, mu, sigma, dst, Tlast)
    dllTerm3 = dllTerm3a - dllTerm3b
    
    logPriorStar = log(pdf(muPrior, muStar))
    logPrior = log(pdf(muPrior, mu))
    
    dll = dllTerm1 + dllTerm2 + dllTerm3 + logPriorStar - logPrior
    dll
end

##############################################################################
# sigma update
##############################################################################

function getSumLessT_sigma(sigma, sigmaStar, mu, tau, theta, dst; 
                                  Tlast=30.)
    tauLessTIndex = find(x -> x<=Tlast, tau) 
    gaussian = Normal(mu, sigma)
    gaussianStar = Normal(mu, sigmaStar)
    
    sumLessT = 0
    for i = tauLessTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j=tauJLessTauIIndex
            tmp = (tau[i] - tau[j])*(pdf(gaussian, dst[i,j]) -
                                   pdf(gaussianStar, dst[i,j]))
            sum_ij = sum_ij + tmp 
        end
        sumLessT = sumLessT + theta*sum_ij
    end
    sumLessT
end


function getSumGreaterT_sigma(sigma, sigmaStar, mu, tau, theta, dst; 
                                     Tlast=30.)

    tauGreaterTIndex = find(x -> x>Tlast, tau)
    gaussian = Normal(mu, sigma)
    gaussianStar = Normal(mu, sigmaStar)
    
    sumGreaterT = 0
    for i=tauGreaterTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j=tauJLessTauIIndex
            sum_ij = sum_ij + 
                (Tlast - tau[j])*(pdf(gaussian, dst[i,j]) -
                                  pdf(gaussianStar, dst[i,j]))
        end
        sumGreaterT = sumGreaterT + theta*sum_ij
    end
    sumGreaterT
end


sigmaPrior = Gamma(sigma_a, sigma_b)

function llRatio_sigma(sigma, sigmaStar, mu, tau, theta, dst; 
                              Tlast=30., sigma_a=0.5, sigma_b=100.)
    
    dllTerm1 = getSumLessT_sigma(sigma, sigmaStar, mu, tau, theta, dst)

    dllTerm2 = getSumGreaterT_sigma(sigma, sigmaStar, mu, tau, theta, dst)
    
    dllTerm3 = getLogLambda(tau, theta, mu, sigmaStar, dst, Tlast) -
        getLogLambda(tau, theta, mu, sigma, dst, Tlast)
    
    priorStar = pdf(sigmaPrior, sigmaStar)
    prior = pdf(sigmaPrior, sigma)
    
    dll = dllTerm1 + dllTerm2 + dllTerm3 + priorStar - prior
    
    dll
end

##############################################################################
# theta update
##############################################################################

function getSumLessT_theta(theta, thetaStar, mu, tau, sigma, dst, Tlast=30.)
    gaussian = Normal(mu, sigma)
    
    tauLessTIndex = find(x -> x<=Tlast, tau)
    
    sumLessT = 0
    for i = tauLessTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j = tauJLessTauIIndex
            sum_ij = sum_ij + 
                (tau[i] - tau[j])*pdf(gaussian, dst[i,j])
        end
        sumLessT = sumLessT + (theta - thetaStar)*sum_ij
    end
    sumLessT
end


function getSumGreaterT_theta(theta, thetaStar, mu, tau, sigma, dst, 
                                     Tlast=30.)
    
    gaussian = Normal(mu, sigma)
    tauGreaterTIndex = find(x -> x>Tlast, tau)
    
    sumGreaterT = 0
    for i=tauGreaterTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j=tauJLessTauIIndex
            sum_ij = sum_ij + 
                (Tlast - tau[j])*pdf(gaussian, dst[i,j])
        end
        sumGreaterT = sumGreaterT + (theta - thetaStar)*sum_ij
    end
    sumGreaterT
end


thetaPrior = Gamma(theta_a, theta_b)

function llRatio_theta(theta, thetaStar, mu, tau, sigma, dst, 
                         Tlast=30, theta_a=0.8, theta_b=10)
    
    dllTerm1 = getSumLessT_theta(theta, thetaStar, mu, tau, sigma, dst, 
                                        Tlast)
    
    dllTerm2 = getSumGreaterT_theta(theta, thetaStar, mu, tau, sigma, dst, 
                                           Tlast)
    
    dllTerm3 = getLogLambda(tau, thetaStar, mu, sigma, dst, Tlast) - 
                    getLogLambda(tau, theta, mu, sigma, dst, Tlast)
      
    logPriorStar = log(pdf(thetaPrior, thetaStar))
    logPrior = log(pdf(thetaPrior, theta))
    
    dll = dllTerm1 + dllTerm2 + dllTerm3 + logPriorStar - logPrior
    dll
end

##############################################################################
# tau update
##############################################################################

function getMinSum(iStar, tauiStar, mu, tau, theta, sigma, dst)
    gaussian = Normal(mu, sigma)
    ti = tau[iStar]
    minTau = min(tauiStar, ti)

    jLessMin = find(x -> x<minTau, tau)
    minSum = 0
    for j = jLessMin
        minSum = minSum + pdf(gaussian, dst[iStar,j])
    end
    minSum = (ti - tauiStar)*theta*minSum
    minSum
end


function getMaxSum(iStar, tauiStar, mu, tau, theta, sigma, dst)
    gaussian = Normal(mu, sigma)
    ti = tau[iStar]
    maxTau = max(tauiStar, ti)
    
    jGreaterMax = find(x -> x>maxTau, tau)
    maxSum = 0
    for j = jGreaterMax
        maxSum = maxSum + pdf(gaussian, dst[iStar,j])
    end
    maxSum = (tauiStar - ti)*theta*maxSum
    maxSum
end


function getSumIntervalBelowTauStar(iStar, tauiStar, mu, tau, theta, sigma, 
                                           dst)
    
    gaussian = Normal(mu, sigma)
    ti = tau[iStar]
    
    jInBetween = find(x -> ((x>ti) & (x<tauiStar)), tau)
    sumBelow = 0
    for j = jInBetween
        sumBelow = sumBelow +
            (2*tau[j] - ti - tauiStar)*pdf(gaussian, dst[iStar,j])
    end
  sumBelow = theta*sumBelow
  sumBelow
end


function getSumIntervalAboveTauStar(iStar, tauiStar, mu, tau, theta, sigma, 
                                           dst)
    gaussian = Normal(mu, sigma)
    ti = tau[iStar]
    
    jInBetween = find(x -> ((x>tauiStar) & (x<ti)), tau)
    sumAbove = 0
    for j = jInBetween
            sumAbove = sumAbove +
                (ti + tauiStar - 2*tau[j])*pdf(gaussian, dst[iStar,j])
    end
    sumAbove = theta*sumAbove
    sumAbove
end


function llRatio_tau(tau, tauiStar, iStar, mu, theta, sigma, dst, Tlast=30.)
    N = length(tau)
    
    dllTerm1 = mu*(tau[iStar] - tauiStar)
    dllTerm2 = getMinSum(iStar, tauiStar, mu, tau, theta, sigma, dst)
    dllTerm3 = getMaxSum(iStar, tauiStar, mu, tau, theta, sigma, dst)
    dllTerm4 = getSumIntervalBelowTauStar(iStar, tauiStar, mu, tau, 
                                          theta, sigma, dst)
    dllTerm5 = getSumIntervalAboveTauStar(iStar, tauiStar, mu, tau, 
                                          theta, sigma, dst)
    
    tauStar = copy(tau)
    tauStar[iStar] = tauiStar
    dllTerm6 = getLogLambda(tauStar, theta, mu, sigma, dst, Tlast)
    dllTerm7 = getLogLambda(tau, theta, mu, sigma, dst, Tlast)
    
    dll = dllTerm1 + dllTerm2 + dllTerm3 + dllTerm4 + dllTerm5 +
        dllTerm6 - dllTerm7
    dll
end

