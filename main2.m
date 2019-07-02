clear all
clc
close all
%%
global kappaInterval
global muInterval

% m: number of state-space discretization
% (0): You need to change this according to your data and problem
timeCounter = 1;
m  = 25;
ix = 1;
subjectNumber = 1;
subjectList = 3;

% load data 
load('example.mat')
% no specific time range for data, set all to 1
tstart = 1; tend = 1; stepsize =1;
%%
% The interval of the 
% (3): This is the discretization of your state space model 
% You also need to set the parameters here based on your data 
% As we were working with the phase data, I used a von mises pdf, therefore
% the parameters are the mean and the concentration parameters (kappa).  
kappaInterval = logspace(0,0.8,m);
muInterval = linspace(-pi, pi, m);
%%
% Setting up the state matrix (mean vs kappa).
for i = 1:m
    for j=1:m      
        StateSpace{ix} = [muInterval(i), kappaInterval(j)];
        ix = ix + 1;
    end
end
numStates = length(StateSpace);
%%
% This is very important how to set the mean and kappa transitions.
% (4): These two parameters are based on your data and the assumptions you make on it.
% Do the changes in your data supposed to occur gradually or rather more sudden. These two 
% determine how fast do you want the changes in your data to occur. They
% are the transition probabilities.
muTransitionKappa = 10;
kappaTransitionSigma = 0.05;
bias_term = 0.0005;

transitionModel.mu.pdf = @(mu_tm1, kappa_tm1) (@(mu)(bias_term + circ_vmpdf(mu, mu_tm1, muTransitionKappa)));  
transitionModel.kappa.pdf = @(mu_tm1, kappa_tm1) (@(kappa)(normpdf(kappa, kappa_tm1, kappaTransitionSigma)));

mTransition = zeros(numStates, numStates);
mat_All_mJointDist = cell(length(tstart:stepsize:tend),1);

for jx = 1:numStates 
    previousState.mu    = StateSpace{jx}(1);
    previousState.kappa = StateSpace{jx}(2);
    for ix = 1:numStates
        targetState.mu    = StateSpace{ix}(1);
        targetState.kappa = StateSpace{ix}(2);
        % The transition probability distributions
        % p(mu_t, kappa_t | mu_t-1, kappa_t-1) = p(mu_t|mu_t-1) * p(kappa_t|kappa_t-1)
        funcMu    = transitionModel.mu.pdf(previousState.mu, previousState.kappa);
        funcKappa = transitionModel.kappa.pdf(previousState.mu, previousState.kappa);
        mTransition(jx, ix) = funcMu(targetState.mu) * funcKappa(targetState.kappa);
    end
    mTransition(jx, :) = mTransition(jx, :) / sum(mTransition(jx, :));
end
%%
figure, 

vData = data; 
TMAX = length(vData);
% Initialization of mJointDist for t=1
mJointDist = zeros(numStates, TMAX);
for ix = 1:numStates
    muState = StateSpace{ix}(1);
    kappaState = StateSpace{ix}(2);

    % how to initialize the probabilities of the first data point
    mJointDist(ix,1) = circ_vmpdf(vData(1),muState,kappaState);
end

mJointDist(:, 1) = mJointDist(:, 1) ./ sum(mJointDist(:,1));

for t=2:(TMAX-1)
    t
    for ix = 1:numStates  
        % p(d_t | mu_t , kappa_t ) ~ vonMises(d_t , mu_t, kappa_t)
        targetState.mu    = StateSpace{ix}(1);
        targetState.kappa = StateSpace{ix}(2);
        predProb = circ_vmpdf(vData(t), targetState.mu , targetState.kappa);

        %vStateProb = zeros(numStates,1);
        %for jx = 1:numStates
        %    transitionProb = mTransition(jx, ix);
        %    vStateProb(jx) = mJointDist(jx, t-1) * transitionProb;
        %end
        % Optimized version of the loop above:
        vStateProb = mTransition(:, ix) .* mJointDist(:, t-1);

        mJointDist(ix,t) = predProb * sum(vStateProb);
        % p(d_t+1|mu_t, kappa_t)
        nextPredProb  = circ_vmpdf(vData(t+1), targetState.mu , targetState.kappa);
        vPredProb(ix) = nextPredProb * mJointDist(ix,t);
    end
    mJointDist(:, t) = mJointDist(:, t) ./ sum(mJointDist(:,t));
    vLikelihood(t)   = sum(vPredProb);
end


vSummedStateLikelihoods = sum(mJointDist, 2);
mSummedStateLikelihoods = zeros(m, m);
ix = 1;
for i = 1:m
    for j=1:m
        mSummedStateLikelihoods(i, j) = vSummedStateLikelihoods(ix);
        ix = ix + 1;
    end
end
subplot(1,3,timeCounter),
contour(mSummedStateLikelihoods);
mat_All_mJointDist{timeCounter} = mJointDist;
timeCounter = timeCounter + 1;

%%
% These two sections down here are showing the likelihood of the state
% transitions. 
%%
%for every trial plot the state space area.
% here we are showing for every 50th trial but you can change this as you
% wish. It is showing the state space for every sample ( which is our
% trial here) for specific time timeInd. Here we had discretized our state
% space into m=25 states, thereby, it shows it as a 25x25 matrix. 
timeCounter = 1;
mJointDist = mat_All_mJointDist{timeCounter};
plotCounter = 1;
figure
for trial =1:50:TMAX
    matStatesPerTrial = zeros(m,m);
    ix = 1;
    for i = 1:m
        for j=1:m
            matStatesPerTrial(i,j) = mJointDist(ix,trial);
            ix = ix+1;
        end
    end    
     subplot(3,7,plotCounter),
     contour(matStatesPerTrial);
     plotCounter = plotCounter +1;
    M(trial) = getframe;
end
timeCounter = timeCounter +1;
%%
% Plotting the marginal likelihood over the concentration or the mean
% parameter. Uncomment the line at 197 for checking the changes in the mean
% parameter. 
% It plots the likelihood of being at any of the concentration or tha mean states. 

mKappaPlot = zeros(m, TMAX);
mMeanPlot = zeros(m,TMAX);
timeCounter = 1;
mJointDist = mat_All_mJointDist{timeCounter};
for trial =1:1:TMAX
    ix = 1;
    for i = 1:m
        for j=1:m
            mKappaPlot(j,trial) = mKappaPlot(j,trial) + mJointDist(ix,trial);
            mMeanPlot(i,trial) = mMeanPlot(i,trial)+ mJointDist(ix,trial); 
            ix = ix+1;
        end
    end
end

figure,
subplot(2,1,1), 
plot(vData,'.'), 
subplot(2,1,2), 
imagesc(1:1:length(mKappaPlot), kappaInterval, mKappaPlot )
xlabel('samples')
ylabel('kappa values')
    