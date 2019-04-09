function [vecParameters] = extract_sigmaKappa(specificConfig)
    numConfig = 1;
    parameterSize = 10;
    vecTransitionKappa = linspace(100,300,parameterSize);
    vecTransitionSigma = linspace(0.005,0.6,parameterSize);
    ix = 1;
    
    for numKappa = 1:parameterSize
        kappa = vecTransitionKappa(numKappa);
        for numSigma = 1:parameterSize
            sigma = vecTransitionSigma(numSigma);
            if specificConfig == numConfig 
                vecParameters(ix,:) = [kappa sigma];
                ix = ix+1;
            end
            %[numConfig kappa sigma]
            numConfig = numConfig +1;

        end
    end
end
