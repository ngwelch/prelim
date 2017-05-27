#!/usr/bin/env julia


include("/Users/nwelch/prelim/src/Julia/optimized_metropolis_sampler.jl");

df = readtable("/Users/nwelch/prelim/data/plantDataTest.csv");
dxixj = convert(Matrix, readtable("/Users/nwelch/prelim/data/dxixjTest.csv"));

mu = 4.85996;
sigma = 1.0;
theta = 0.1;

run_mcmc_sampler(mu, theta, sigma, df, dxixj, trials=65000,
                     mu_shape=1.307877, mu_rate=0.06334962,
                     theta_shape=0.8, theta_rate=0.1,
                     sigma_shape=0.5, sigma_scale=100,
                     v_mu=0.005, v_sigma=0.55, v_theta=0.05)

#v_mu: 0.005    |-      |-      |-  |new th
#v_sg: 0.005    |0.01   |0.1    |1  |new th
#v_th: 0.05     |-      |-      |-  |new th

#previous values, sigma autocorrelation too high
#v_mu: 0.005    |-      |-      |-
#v_sg: 0.5      |0.005  |0.5    |0.75
#v_th: 0.05     |-      |-      |-
