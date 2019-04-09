function [mJointDist] = normalize_mJointDist(mJointDist, TMAX)
% normalizing the mJointDist
    parfor t=1:TMAX
        if sum(mJointDist(:, t)) > 0
            mJointDist(:, t) = mJointDist(:, t) / sum(mJointDist(:, t));
        end
    end
end
