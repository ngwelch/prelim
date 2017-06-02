#!/usr/bin/env julia

home = homedir()
include(home*"/prelim/src/Julia/optimized_metropolis_sampler.jl");

df = readtable(home*"/prelim/data/plantDataTest.csv");
dxixj = convert(Matrix, readtable(home*"/prelim/data/dxixjTest.csv"));

mu = 0.003;
sigma = 1.0;
theta = 0.1;

run_mcmc_sampler(mu, theta, sigma, df, dxixj, trials=20000,
                     mu_shape=1.307877, mu_rate=1.900488,   #2 plants/wk
                     theta_shape=0.85, theta_rate=0.145,    #94% CI
                     sigma_shape=0.75, sigma_scale=8,       #0.1-50m
                     v_mu=0.001, v_sigma=0.05, v_theta=0.05,
                     metricsFile="mcmc_chain_metrics_julia_20k",
                     chainFile="mcmc_chain_julia_20k")

#v_mu: 0.005    |-      |-      |-  |new th
#v_sg: 0.005    |0.01   |0.1    |1  |new th
#v_th: 0.05     |-      |-      |-  |new th

#v_mu: 0.005    |-      |-      |-      |new mu_rate
#v_sg: 0.5      |0.005  |0.5    |0.75   |new th shape & rate
#v_th: 0.05     |-      |-      |-      |new shape & scale

#v_mu: 0.005    |0.007  |0.005  |0.005
#v_sg: 0.55     |0.575  |0.05   |
#v_th:          |       |       |


