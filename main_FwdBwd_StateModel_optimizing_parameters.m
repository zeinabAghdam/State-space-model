% The forward-backward model 

% One problem of the forward method as described in the paper ( see Readme) 
% is that it can heavily depend on the initialization parameter at t=1. 
% Depending on how the state likelihoods for the first observation is
% initialized, state likelihoods at the subsequent times can be far from
% the ideal. 

%  A similar effect is observed if there is a sudden change in the
%  underlying states at a later time in the data or if the initial data
%  points are outliers (noise). Because the forward method does not look
%  ahead of the current time, it adapts to changes in the data only with a
%  delay. 

% We add a backward pass to compensate for these
% effects, and to remove any influence of the direction of time on the 
% the results

% We also mitigate the initialization problem by using the result of the forward pass at tN to initialize
% the backward pass, and vice-verse. Unless the data distribution drastically changes close to the end
% of the data tN, the forward method usually converges sufficiently at that point, thus providing a very
% stable initialization for the backward pass. After the backward pass is done, we run the forward pass
% for a second time, using the result at t1 from the backward pass for its initialization. The results
% of the second forward pass together with the results of the backward pass are used in analysis. The
% results of the first, randomly initialized forward pass are discarded in the final result to remove the
% influence of the original random initialization
%%
clc,
clear all,
close all,

global vSubjectNumbers
global vCategory
global intensityList
global vDataVariables
global directoryAlloc
global m
global tstart
global tend
global stepsize
global muInterval
global kappaInterval
global transitionModel
global muTransitionKappa
global kappaTransitionSigma
global numStates
global timeLength
global vecTransitionKappa
global vecTransitionSigma
global interVariance
global intraVariance

%% Load data 
load('example.mat')
%% Parameter initialization
% m: number of state discretization
% tstart:stepsize:tend : the time interval to be analyzed
% All the parameters in this section have be adjusted according to your
% data and problem 

m =  20;
stepsize = 10;

kappaInterval = logspace(0,1.6,m);
muInterval = linspace(-pi, pi, m);

%parameterSize = 8;
%vecTransitionKappa = linspace(200,500,parameterSize);
%vecTransitionSigma = linspace(0.1,0.6,parameterSize);
%% Initialize the model parameters.
mat_all_mJointDist = cell(length(tstart:stepsize:tend),1);
mat_all_summedLikelihoods = cell(length(tstart:stepsize:tend),1);
%%
% We don't know the optimum parameters, kappaTransition and the
% meanTransition. There are various methods, and you can define your own
% objective function for doing so, however, one way was based on maximizing
% the fisher's criteria, in the sense that the parameters are adjusted
% according to differences betweeen the groups of data that we're
% analyzing.

kappa = 10;
sigma = 0.6071;
sigma = 0.2;

[numStates, mTransition,StateSpace] = modelStateMatrixInitialization(kappa, sigma);

vData = data;
TMAX = length(vData);
%vData = vData(1:TMAX);
% vData = locallyCenterData(vData,100);

mJointDistFw = compute_mjointDist( StateSpace, mTransition, vData );
mJointDistBw = compute_mjointDist( StateSpace, mTransition, vData(end:-1:1), mJointDistFw(:, end-1));
mJointDistFw = compute_mjointDist( StateSpace, mTransition, vData, mJointDistBw(:, end-1) );
mJointDist = mJointDistFw .* mJointDistBw(:, end:-1:1);

%normalizing the mJointDist matrix
mJointDist = normalize_mJointDist(mJointDist, TMAX);
vSummedStateLikelihoods = sum(mJointDist, 2);
mSummedStateLikelihoods = compute_mSummedStateLikelihoods(vSummedStateLikelihoods);   

%%
% Plotting the marginal likelihood over the concentration or the mean
% parameter. Uncomment the line at 197 for checking the changes in the mean
% parameter. 
% It plots the likelihood of being at any of the concentration or tha mean states. 
mKappaPlot = zeros(m, TMAX);
mMeanPlot = zeros(m,TMAX);
timeCounter = 1;
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
imagesc(1:1:length(mKappaPlot), kappaInterval, mKappaPlot ), colormap hot
xlabel('samples')
ylabel('kappa values')