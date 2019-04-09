function [ mSummedStateLikelihoods ] = compute_mSummedStateLikelihoods( vSummedStateLikelihoods )
    global m
    mSummedStateLikelihoods = zeros(m, m);
    ix = 1;
    for i = 1:m
        for j=1:m
            mSummedStateLikelihoods(i, j) = vSummedStateLikelihoods(ix);
            ix = ix + 1;
        end
    end
end

