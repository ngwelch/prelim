#!/usr/bin/env julia

home = homedir()
include(home*"/prelim/src/Julia/optimized_metropolis_sampler.jl");

df = readtable(home*"/prelim/data/plantDataTest.csv");
dxixj = convert(Matrix, readtable(home*"/prelim/data/dxixjTest.csv"));

mu = 0.003;
sigma = 1.2;
theta = 0.12;

run_mcmc_sampler(mu, theta, sigma, df, dxixj, trials=10000,
                     mu_shape=0.971366, mu_rate=0.6881532,  #2 plants/wk
                     theta_shape=0.85, theta_rate=0.145,    #94% CI
                     sigma_shape=0.75, sigma_scale=8,       #0.1-50m
                     v_mu=0.003, v_sigma=0.9, v_theta=0.05,
                     metricsFile="10k_25pct_aws_metrics",
                     chainFile="10k_25pct_aws_chain")

#v_mu: 0.005	|0.0005	|0.001 	56%	|0.002	40%
#v_sg: 0.05	|0.1	|0.2	73%	|0.4	49%
#v_th: 0.01	|0.05	|	28%	|	30%

#v_mu: 0.004	20%	| 	25%	|	26%
#v_sg: 0.8	36%	|0.1	88%	|0.5	48%
#v_th:		30%	|	29%	|	31%

