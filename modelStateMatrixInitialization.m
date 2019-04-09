function [numStates, mTransition, StateSpace] = modelStateMatrixInitialization( muTransitionKappa, kappaTransitionSigma )

% compute the state space matrix.
% compute the transition matrix.

StateSpace = computeStateSpace();
numStates = length(StateSpace);
mTransition = computeTransitionModel(StateSpace, numStates, muTransitionKappa, kappaTransitionSigma);

end

