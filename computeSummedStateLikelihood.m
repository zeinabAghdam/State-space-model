function [ mSummedStateLikelihoods ] = computeSummedStateLikelihood( mJointDist, n )
    global grid_size
    vSummedStateLikelihoods = sum(mJointDist, n);
    mSummedStateLikelihoods = zeros(grid_size, grid_size);
    ix = 1;
    if n==2
        for i = 1:grid_size
            for j=1:grid_size
                mSummedStateLikelihoods(i, j) = vSummedStateLikelihoods(ix);
                ix = ix + 1;
            end
        end
    else
      ix = 1;
      for i = 1:grid_size
            for j=1:grid_size
                mSummedStateLikelihoods(j, i) = vSummedStateLikelihoods(ix);
                ix = ix + 1;
            end
      end
    end
end

