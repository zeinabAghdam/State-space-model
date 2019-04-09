function [mJointDist, vLikelihood] = compute_mjointDist( StateSpace, mTransition, vData, vInitDist )
    numStates = length(StateSpace);
    TMAX = length(vData);
    % Initialization of mJointDist for t=1
    mJointDist = zeros(numStates, TMAX);
    if nargin == 4
        mJointDist(:, 1) = vInitDist;
    else
        for ix = 1:numStates
            muState = StateSpace{ix}(1);
            kappaState = StateSpace{ix}(2);
            mJointDist(ix,1) = circ_vmpdf(vData(1),muState,kappaState);
            % This sum is probably wrong, but maybe it works
            % mJointDist(ix,1) = sum(circ_vmpdf(vData(1:50),muState,kappaState));
        end
        mJointDist(:, 1) = mJointDist(:, 1) ./ sum(mJointDist(:,1));
    end
    
    vLikelihood = zeros(TMAX, 1);
    vPredProb = zeros(numStates, 1);
    for t=2:(TMAX-1)
        
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
end

