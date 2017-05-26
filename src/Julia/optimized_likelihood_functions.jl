#!/usr/bin/env julia

using DataFrames;
using Distributions;

##############################################################################
# theta * f(x_i - x_j | sigma)
##############################################################################

function getThetafxixj(theta, sigma, dxixj)
    dstSq = dxixj.*dxixj
    twoSigmaSq = 2*(sigma^2)
    thetafxixj = (theta/(twoSigmaSq*pi))*exp(-dstSq/twoSigmaSq)
    thetafxixj
end

##############################################################################
# log-intensity
##############################################################################

function getLogLambda(tau, mu, thetafxixj, Tlast=30.)
    tauLessTIndex = find(x-> x<=Tlast, tau)

    sumLog = 0
    for i = tauLessTIndex
        sumTauJLessTauI = 0
        tauJLessTauIIndex = find(x-> x < tau[i], tau)
        for j = tauJLessTauIIndex
            sumTauJLessTauI = sumTauJLessTauI + thetafxixj[i,j]
        end
        sumLog = sumLog + log(mu + sumTauJLessTauI)
    end
    sumLog
end

##############################################################################
# mu update
##############################################################################

function llRatio_mu(mu, muStar, tau, thetafxixj, 
                        Tlast=30.0, mu_shape=0.7, mu_rate=0.004)
    muPrior = Gamma(mu_shape, 1./mu_rate)

    tauLessTIndex = find(x-> x<=Tlast, tau)
    tauGreaterT = count(x-> x>Tlast, tau)

    dmu = mu - muStar

    dllTerm1 = sum(tau[tauLessTIndex]*dmu)
    dllTerm2 = tauGreaterT * Tlast * dmu

    dllTerm3a = getLogLambda(tau, muStar, thetafxixj, Tlast)
    dllTerm3b = getLogLambda(tau, mu, thetafxixj, Tlast)
    dllTerm3 = dllTerm3a - dllTerm3b

    logPriorStar = log(pdf(muPrior, muStar))
    logPrior = log(pdf(muPrior, mu))

    dll = dllTerm1 + dllTerm2 + dllTerm3 + logPriorStar - logPrior
    dll
end

##############################################################################
# sigma update
##############################################################################

function getSumLessT_sigma(dthetafxixj, mu, tau, Tlast=30.0)
    tauLessTIndex = find(x -> x<=Tlast, tau)

    sumLessT = 0
    for i = tauLessTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j=tauJLessTauIIndex
            #reversing i & j to make overall term positive
            tmp = (tau[j] - tau[i])*dthetafxixj[i,j]
            sum_ij = sum_ij + tmp
        end
        sumLessT = sumLessT + sum_ij
    end
    sumLessT
end


function getSumGreaterT_sigma(dthetafxixj, mu, tau, Tlast=30.)

    tauGreaterTIndex = find(x -> x>Tlast, tau)

    sumGreaterT = 0
    for i=tauGreaterTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j=tauJLessTauIIndex
            #reversing j & T to make overall term positive
            sum_ij = sum_ij +
                (tau[j]-Tlast)*dthetafxixj[i,j]
        end
        sumGreaterT = sumGreaterT + sum_ij
    end
    sumGreaterT
end


function llRatio_sigma(sigma, sigmaStar, thetafxixj, 
        mu, tau, Tlast=30.0, sigma_shape=0.5, sigma_scale=100.)
    
    thetafxixjStar = getThetafxixj(theta, sigmaStar, dxixj);
    dthetafxixj = thetafxixjStar - thetafxixj
    sigmaPrior = Gamma(sigma_shape, sigma_scale)

    dllTerm1 = getSumLessT_sigma(dthetafxixj, mu, tau, Tlast)

    dllTerm2 = getSumGreaterT_sigma(dthetafxixj, mu, tau, Tlast)

    dllTerm3 = getLogLambda(tau, mu, thetafxixjStar, Tlast) -
        getLogLambda(tau, mu, thetafxixj, Tlast)

    logPriorStar = log(pdf(sigmaPrior, sigmaStar))
    logPrior = log(pdf(sigmaPrior, sigma))

    dll = dllTerm1 + dllTerm2 + dllTerm3 + logPriorStar - logPrior

    dll
end

##############################################################################
# theta update
##############################################################################

function getSumLessT_theta(theta, thetaStar, 
        thetafxixj, thetafxixjStar, dthetafxixj,
        mu, tau, sigma, dxixj, Tlast=30.0)
    
    tauLessTIndex = find(x -> x<=Tlast, tau)

    sumLessT = 0
    for i = tauLessTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j = tauJLessTauIIndex
            sum_ij = sum_ij +
                (tau[i] - tau[j])*dthetafxixj[i,j]
        end
        sumLessT = sumLessT + sum_ij
    end
    sumLessT
end

function getSumGreaterT_theta(theta, thetaStar, 
        thetafxixj, thetafxixjStar, dthetafxixj,
        mu, tau, sigma, dxixj, Tlast=30.0)
    
    tauGreaterTIndex = find(x -> x>Tlast, tau)

    sumGreaterT = 0
    for i=tauGreaterTIndex
        sum_ij = 0
        tauJLessTauIIndex = find(x-> x<tau[i], tau)
        for j=tauJLessTauIIndex
            sum_ij = sum_ij + (Tlast - tau[j])*dthetafxixj[i,j]
        end
        sumGreaterT = sumGreaterT + sum_ij
    end
    sumGreaterT
end


function llRatio_theta(theta, thetaStar, mu, tau, sigma, dxixj, 
                         Tlast=30.0, theta_shape=0.8, theta_rate=10.0)

    thetafxixj = getThetafxixj(theta, sigma, dxixj)
    thetafxixjStar = getThetafxixj(thetaStar, sigma, dxixj)
    dthetafxixj = thetafxixj - thetafxixjStar
    
    thetaPrior = Gamma(theta_shape, 1./theta_rate)
    
    dllTerm1 = getSumLessT_theta(theta, thetaStar, 
        thetafxixj, thetafxixjStar, dthetafxixj, 
        mu, tau, sigma, dxixj, Tlast)

    dllTerm2 = getSumGreaterT_theta(theta, thetaStar, 
        thetafxixj, thetafxixjStar, dthetafxixj, 
        mu, tau, sigma, dxixj, Tlast)

    dllTerm3 = getLogLambda(tau, mu, thetafxixjStar, Tlast) -
                    getLogLambda(tau, mu, thetafxixj, Tlast)

    logPriorStar = log(pdf(thetaPrior, thetaStar))
    logPrior = log(pdf(thetaPrior, theta))

    dll = dllTerm1 + dllTerm2 + dllTerm3 + logPriorStar - logPrior
    dll
end

##############################################################################
# tau update
##############################################################################

function getMinSum(iStar, tauiStar, tau, thetafxixj)
    ti = tau[iStar]
    minTau = min(tauiStar, ti)

    jLessMin = find(x -> x<minTau, tau)
    minSum = 0
    for j = jLessMin
        minSum = minSum + thetafxixj[iStar,j]
    end
    minSum = (ti - tauiStar)*minSum
    minSum
end


function getMaxSum(iStar, tauiStar, tau, thetafxixj)
    ti = tau[iStar]
    maxTau = max(tauiStar, ti)

    jGreaterMax = find(x -> x>maxTau, tau)
    maxSum = 0
    for j = jGreaterMax
        maxSum = maxSum + thetafxixj[iStar,j]
    end
    maxSum = (tauiStar - ti)*maxSum
    maxSum
end


function getSumIntervalBelowTauStar(iStar, tauiStar, tau, thetafxixj)

    ti = tau[iStar]

    jInBetween = find(x -> ((x>ti) & (x<tauiStar)), tau)
    sumBelow = 0
    for j = jInBetween
        sumBelow = sumBelow + (2*tau[j] - ti - tauiStar)*thetafxixj[iStar,j]
    end
  sumBelow
end


function getSumIntervalAboveTauStar(iStar, tauiStar, tau, thetafxixj)
    ti = tau[iStar]

    jInBetween = find(x -> ((x>tauiStar) & (x<ti)), tau)
    sumAbove = 0
    for j = jInBetween
            sumAbove = sumAbove + (ti + tauiStar - 2*tau[j])*thetafxixj[iStar,j]
    end
    sumAbove
end


function llRatio_tau(tau, tauiStar, iStar, mu, thetafxixj, lastLogLambda, Tlast=30.0)
    
    N = length(tau)

    dllTerm1 = mu*(tau[iStar] - tauiStar)
    dllTerm2 = getMinSum(iStar, tauiStar, tau, thetafxixj)
    dllTerm3 = getMaxSum(iStar, tauiStar, tau, thetafxixj)
    dllTerm4 = getSumIntervalBelowTauStar(iStar, tauiStar, tau, thetafxixj)
    dllTerm5 = getSumIntervalAboveTauStar(iStar, tauiStar, tau, thetafxixj)

    tauStar = copy(tau)
    tauStar[iStar] = tauiStar
    dllTerm6 = getLogLambda(tauStar, mu, thetafxixj, Tlast) - lastLogLambda

    dll = dllTerm1 + dllTerm2 + dllTerm3 + dllTerm4 + dllTerm5 + dllTerm6
    dll
end

