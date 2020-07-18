# State Space Model

The files provided detect the changes in the signal using a Bayesian Change Point Model. The states ```mu``` and ```kappa``` are the mean and standard deviation of a von Mises distribution that have been discretized equally. 

The parameters ```muTransitionKappa``` and ```kappaTransitionSigma``` in ```main2.m``` file determine how fast the transition between the states occur.  

