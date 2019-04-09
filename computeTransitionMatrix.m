function [StateSpace, mTransition] = computeTransitionMatrix()

    global transitionModel
    global kappaInterval
    global muInterval
    global m 
    ix = 1;
    for i = 1:m
        for j=1:m      
            StateSpace{ix} = [muInterval(i), kappaInterval(j)];
            ix = ix + 1;
        end
    end
    % StateSpace = @(i) (muIndex = floor(i / m)+1; kappaIndex = mod(i, m)+1; [muInterval(muIndex), kappaInderval(kappaIndex)]);
    numStates = length(StateSpace);
       
    mTransition = zeros(numStates, numStates);

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
end