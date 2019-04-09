function [ mTransition ] = computeTransitionModel( StateSpace, numStates, muTransitionKappa, kappaTransitionSigma )

% compute the transition model matrix
global transitionModel

muTransitionConst = 0.0005;
kappaTransitionConst = 0.0;

transitionModel.mu.pdf = @(mu_tm1, kappa_tm1) (@(mu)(muTransitionConst + circ_vmpdf(mu, mu_tm1, muTransitionKappa)));  
transitionModel.kappa.pdf = @(mu_tm1, kappa_tm1) (@(kappa)(kappaTransitionConst +normpdf(kappa, kappa_tm1, kappaTransitionSigma)));

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
            mTransition(jx, ix) = funcMu(StateSpace{ix}(1)) * funcKappa(StateSpace{ix}(2));
        end
        mTransition(jx, :) = mTransition(jx, :) / sum(mTransition(jx, :));
    end

end

