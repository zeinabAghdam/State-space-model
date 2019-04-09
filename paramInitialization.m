function [  ] = paramInitialization( input_args )

global groupList
global tstart
global tend
global stepsize
global m
global kappaInterval
global muInterval
global dataVariableName
global muTransitionKappa
global kappaTransitionSigma
global transitionModel

groupList = [6, 8, 8];
tstart = 44;
tend = 65;
stepsize = 3;
m  = 30;
kappaInterval = logspace(0,0.9,m);
muInterval = linspace(-pi, pi, m);
dataVariableName = 'phase_tinn_af_2_nlm';
muTransitionKappa = 100;
kappaTransitionSigma = 0.05;
transitionModel.mu.pdf = @(mu_tm1, kappa_tm1) (@(mu)(0.0005 + circ_vmpdf(mu, mu_tm1, muTransitionKappa)));  
transitionModel.kappa.pdf = @(mu_tm1, kappa_tm1) (@(kappa)(normpdf(kappa, kappa_tm1, kappaTransitionSigma)));


end

