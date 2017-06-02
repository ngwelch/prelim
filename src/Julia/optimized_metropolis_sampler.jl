#!/usr/bin/env julia

using DataFrames;
using Distributions;

function run_mcmc_sampler(mu, theta, sigma, df, dxixj; 
                              trials=100, Tlast=30.0,
                              mu_shape=0.7, mu_rate=0.004,
                              theta_shape=0.8, theta_rate=10.0,
                              sigma_shape=0.5, sigma_scale=100.0,
                              v_mu=0.0005, v_sigma=0.05, v_theta=0.005
                              metricsFile="mcmc_chain_metrics_julia",
                              chainFile="mcmc_chain_julia")

    @everywhere include("/Users/nwelch/prelim/src/Julia/optimized_likelihood_functions.jl");

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
    tau = convert(SharedArray, tau)
    tauLowerBound = df[:tauLowerBound];
    tauUpperBound = df[:tauUpperBound];

    chain = Array{Float64}(trials, (3+infectionCount));
    chain[1,1:infectionCount] = copy(tau[infectionIndex]);
    chain[1,(infectionCount+1):(infectionCount+3)] = copy([mu sigma theta]);

    U = Uniform();

    for r=2:trials
    
        # sample U to generate the acceptance criteria
        u = rand(U, (infectionCount+3))
        logU = log(u)
   
    
        # update mu
        thetafxixj = getThetafxixj(theta, sigma, dxixj)
        mu_start = now()
        J_mu = Normal(mu, v_mu)
        muStar = rand(J_mu, 1)[1]
        if muStar>0.0
            logLRatio_mu = llRatio_mu(mu, muStar, tau, thetafxixj, 
                                          Tlast, mu_shape, mu_rate)
  
            if logU[1] < logLRatio_mu
                mu = copy(muStar)
                accept[1] = accept[1] + 1
            end
        end
        mu_time = mu_time + (now() - mu_start);
  
    
        # update sigma
        thetafxixj = getThetafxixj(theta, sigma, dxixj)
        sigma_start = now()
        J_sigma = Normal(sigma, v_sigma)
        sigmaStar = rand(J_sigma, 1)[1]
        if sigmaStar>0.0
        
            logLRatio_sigma = llRatio_sigma(sigma, sigmaStar, thetafxixj, 
                                                   mu, tau, Tlast, 
                                                   sigma_shape, sigma_scale)

            if logU[2] < logLRatio_sigma  
                sigma = copy(sigmaStar)
                accept[2] = accept[2] + 1
            end
        end
        sigma_time = sigma_time + (now() - sigma_start);
    
    
        # update theta
        thetafxixj = getThetafxixj(theta, sigma, dxixj)
        theta_start = now()
        J_theta = Normal(theta, v_theta)
        thetaStar = rand(J_theta, 1)[1]
        if thetaStar>0.
            logLRatio_theta = llRatio_theta(theta, thetaStar, mu, tau, sigma, 
                                                   dxixj, Tlast, 
                                                   theta_shape, theta_rate)

            if logU[3] < logLRatio_theta 
                theta = copy(thetaStar)
                accept[3] = accept[3] + 1
            end 
        end
        theta_time = theta_time + (now() - theta_start);
    
    
        #tau update
        tau_start = now()
        gaussian_tau = Normal(mu, sigma)
        lastLogLambda = getLogLambda(tau, mu, thetafxixj, Tlast)
        for iStar=infectionIndex
        #@sync @parallel for iStar=infectionIndex
        
            J_taui = Normal(tau[iStar], 1)
            tauiStar = rand(J_taui, 1)[1]

            lowerBound = tauLowerBound[iStar]
            upperBound = tauUpperBound[iStar]
            
            if (tauiStar<upperBound) & (tauiStar>lowerBound)
                logLRatio_tau = llRatio_tau(tau, tauiStar, iStar, mu, 
                                                thetafxixj, lastLogLambda, Tlast)
                index = iStar-firstInfection+1
                if logU[(index+3)] < logLRatio_tau
                    tau[iStar] = copy(tauiStar)
                    lastLogLambda = getLogLambda(tau, mu, thetafxixj, Tlast)
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

    writetable("/Users/nwelch/prelim/data/"*metricsFile*".csv", 
        metrics, separator = ',', header = true);

    writetable("/Users/nwelch/prelim/data/"*chainFile*".csv", 
        chainDF, separator = ',', header = true);

end
