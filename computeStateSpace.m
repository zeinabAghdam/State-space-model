function [ StateSpace ] = computeStateSpace()
% ------------------------------------------
% construct the state space matrix.
% ------------------------------------------
global m
global muInterval
global kappaInterval
ix = 1;
for i = 1:m
    for j=1:m      
        StateSpace{ix} = [muInterval(i), kappaInterval(j)];
        ix = ix + 1;
    end
end

end

