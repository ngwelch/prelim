#!/usr/bin/env julia

using DataFrames;
using Distributions;

function run_mcmc_sampler(mu, theta, sigma, df, dst; 
                              trials=100, Tlast=30.0,
                              mu_a=0.7, mu_b=0.004,
                              sigma_a=0.5, sigma_b=100.0,
                              theta_a=0.8, theta_b=10.0)

    @everywhere include("/Users/nwelch/prelim/src/Julia/likelihood_functions.jl");

    #setup 
    infectionCount = sum(df[:t6]);
    accept = fill(0, 4);
    accept[4] = trials;

    zero = now()-now()
    mu_time = zero;
    sigma_time = zero; 
    tau_time = zero;
    theta_time = zero;

    firstInfection = indmin(df[:tau]);
    lastInfection = length(df[:tau]);
    infectionIndex = firstInfection:lastInfection;

    tau = convert(Array, copy(df[:tau]));
    tauLowerBound = df[:tauLowerBound];
    tauUpperBound = df[:tauUpperBound];

    chain = Array{Float64}(trials, (3+infectionCount));
    chain[1,1:infectionCount] = copy(tau[infectionIndex]);
    chain[1,(infectionCount+1):(infectionCount+3)] = copy([mu sigma theta]);

    U = Uniform();

    for r=2:trials
    
        tau[infectionIndex] = copy(chain[r-1, 1:infectionCount])
        mu = copy(chain[(r-1), (infectionCount+1)])
        sigma = copy(chain[(r-1), (infectionCount+2)])
        theta = copy(chain[(r-1), (infectionCount+3)])
    
    
        # sample U to generate the acceptance criteria
        u = rand(U, (infectionCount+3))
        logU = log(u)
    
    
        # shared function arguments
        tauLessTIndex = find(x-> x<=Tlast, tau)
        tauGreaterTIndex = find(x-> x>Tlast, tau)
    
    
        # update mu
        mu_start = now()
        J_mu = Normal(mu, 0.0005)
        muStar = rand(J_mu, 1)[1]
        if muStar>0.0
            logLRatio_mu = llRatio_mu(mu, muStar, theta, sigma, tau, dst)
                              
            if logU[1] < logLRatio_mu
                mu = copy(muStar)
                accept[1] = accept[1] + 1
            end
        end
        mu_time = mu_time + (now() - mu_start);
  
    
        # update sigma
        sigma_start = now()
        J_sigma = Normal(sigma, 0.05)
        sigmaStar = rand(J_sigma, 1)[1]
        if sigmaStar>0.0
        
            gaussian_s = Normal(mu, sigma)
            gaussianStar_s = Normal(mu, sigmaStar)
        
            logLRatio_sigma = llRatio_sigma(sigma, sigmaStar,
                gaussian_s, gaussianStar_s, mu, tau, theta, dst)
        
            if logU[2] < logLRatio_sigma  
                sigma = copy(sigmaStar)
                accept[2] = accept[2] + 1
            end
        end
        sigma_time = sigma_time + (now() - sigma_start);
    
    
        # update theta
        theta_start = now()
        J_theta = Normal(theta, 0.005)
        thetaStar = rand(J_theta, 1)[1]
        if thetaStar>0.
            gaussian_t = Normal(mu, sigma)
            logLRatio_theta = llRatio_theta(theta, thetaStar,
                mu, tau, sigma, gaussian_t, dst)
        
            if logU[3] < logLRatio_theta 
                theta = copy(thetaStar)
                accept[3] = accept[3] + 1
            end 
        end
        theta_time = theta_time + (now() - theta_start);
    
    
        #tau update
        tau_start = now()
        gaussian_tau = Normal(mu, sigma)
        curTau = copy(tau)
        tau = convert(SharedArray, tau)
        @sync @parallel for iStar=infectionIndex
        
            J_taui = Normal(tau[iStar], 1)
            tauiStar = rand(J_taui, 1)[1]

            lowerBound = tauLowerBound[iStar]
            upperBound = tauUpperBound[iStar]
            
            if (tauiStar<upperBound) & (tauiStar>lowerBound)
                logLRatio_tau = llRatio_tau(curTau, tauiStar, iStar, mu, theta,
                                                    sigma, gaussian_tau, dst)
                index = iStar-firstInfection+1
                if logU[(index+3)] < logLRatio_tau
                   tau[iStar] = copy(tauiStar)
                end
            end
        end
        tau_time = tau_time + (now() - tau_start);
        
        chain[r,1:infectionCount] = copy(tau[infectionIndex]);
        chain[r,(infectionCount+1):(infectionCount+3)] = copy([mu sigma theta]);
    end

    metrics = DataFrame(
        mu_accept=accept[1],
        sigma_accept=accept[2],
        theta_accept=accept[3],
        mu_time=mu_time,
        sigma_time=sigma_time,
        theta_time=theta_time,
        tau_time=tau_time);

    chainDF = DataFrame(chain);
    chainColNames = [Symbol("tau$i") for i in 1:(infectionCount)];
    chainColNames = append!(chainColNames, ["mu","sigma","theta"]);
    names!(chainDF, chainColNames);

    writetable("/Users/nwelch/prelim/data/mcmc_chain_metrics_julia.csv", 
        metrics, separator = ',', header = true);

    writetable("/Users/nwelch/prelim/data/mcmc_chain_julia.csv", 
        chainDF, separator = ',', header = true);

    metrics
end
