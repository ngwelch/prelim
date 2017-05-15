#!/usr/bin/env julia

mu = 0.003;
sigma = 5.0;
theta = 0.105;

include("/Users/nwelch/prelim/src/Julia/metropolis_sampler.jl");

df = readtable("/Users/nwelch/prelim/data/plantDataTest.csv");
dst = readtable("/Users/nwelch/prelim/data/plantDistanceTest.csv");

run_mcmc_sampler(mu, theta, sigma, df, dst, trials=10000)
